"""Regression tests for main.py's mtime-invalidated intermediate-file cache
(load_inputs()'s caching contract -- see its docstring, and CLAUDE.md's
"intermediate/<mol>_<normal|emit>_data.txt" note).

Run from ``Github/scoring-functions/``:
    py -m pytest tests/test_intermediate_cache.py
    py tests/test_intermediate_cache.py       (standalone; no pytest needed)

Exercises against the real 'water' normal-mode data already in this repo
(data/logs/water.log + data/gjf/water.gjf), matching the convention the rest
of tests/ already uses (real molecules, not synthetic fixtures). Each test
snapshots and restores the intermediate file's bytes/mtime and the log's
mtime, so running this file has no lasting effect on repo state.
"""
import os
import sys
from unittest import mock

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
DATA_DIR = os.path.join(ROOT, "data")

import main                                                          # noqa: E402
from src.parser import GaussianParser                                # noqa: E402

MOL, MODE = "water", "normal"


def _inter_path():
    dirs = main.resolve_dirs(DATA_DIR)
    return main.intermediate_path(dirs, MOL, MODE)


def _log_path():
    return os.path.join(DATA_DIR, "logs", f"{MOL}.log")


class _StateGuard:
    """Snapshot + restore the intermediate file's bytes/mtime and the log's
    mtime around a test, so this file's tests don't leave the repo dirty."""

    def __enter__(self):
        self.inter_path = _inter_path()
        self.log_path = _log_path()
        self.had_inter = os.path.exists(self.inter_path)
        if self.had_inter:
            with open(self.inter_path, "rb") as f:
                self.inter_bytes = f.read()
            self.inter_mtime = os.path.getmtime(self.inter_path)
        self.log_mtime = os.path.getmtime(self.log_path)
        return self

    def __exit__(self, *exc):
        if self.had_inter:
            with open(self.inter_path, "wb") as f:
                f.write(self.inter_bytes)
            os.utime(self.inter_path, (self.inter_mtime, self.inter_mtime))
        elif os.path.exists(self.inter_path):
            os.remove(self.inter_path)
        os.utime(self.log_path, (self.log_mtime, self.log_mtime))


def test_first_call_writes_mode_suffixed_intermediate():
    """A cache miss parses from source and writes exactly
    data/intermediate/<mol>_<mode_type>_data.txt (not the old unsuffixed
    '<mol>_data.txt', which collided between normal/EMIT for the same mol)."""
    with _StateGuard():
        if os.path.exists(_inter_path()):
            os.remove(_inter_path())
        main.load_inputs(MOL, MODE, DATA_DIR)
        assert os.path.exists(_inter_path())
        assert _inter_path().endswith(f"{MOL}_{MODE}_data.txt")


def test_cache_hit_skips_reparse():
    """Second call with nothing changed must NOT re-parse the log -- it
    should come straight from the cached intermediate instead."""
    with _StateGuard():
        if os.path.exists(_inter_path()):
            os.remove(_inter_path())
        raw1, _ = main.load_inputs(MOL, MODE, DATA_DIR)

        with mock.patch.object(
            GaussianParser, "parse",
            side_effect=AssertionError("cache hit must not re-parse the log"),
        ):
            raw2, _ = main.load_inputs(MOL, MODE, DATA_DIR)

        assert raw2["atoms"] == raw1["atoms"]
        assert len(raw2["modes"]) == len(raw1["modes"])
        assert raw2["bonds"] == raw1["bonds"]


def test_stale_source_forces_reparse():
    """Touching the log to be newer than the cached intermediate must trigger
    a real re-parse, not a silent (stale) cache hit."""
    with _StateGuard():
        main.load_inputs(MOL, MODE, DATA_DIR)  # ensure a cache exists
        inter_path = _inter_path()
        future = os.path.getmtime(inter_path) + 5
        os.utime(_log_path(), (future, future))

        call_count = {"n": 0}
        real_parse = GaussianParser.parse

        def _counting_parse(self, *a, **kw):
            call_count["n"] += 1
            return real_parse(self, *a, **kw)

        with mock.patch.object(GaussianParser, "parse", _counting_parse):
            main.load_inputs(MOL, MODE, DATA_DIR)

        assert call_count["n"] == 1, "stale cache should have forced exactly one re-parse"
        assert os.path.getmtime(inter_path) >= os.path.getmtime(_log_path())


def test_use_cache_false_forces_reparse_even_if_fresh():
    """use_cache=False bypasses an otherwise-fresh cache (still refreshes it)."""
    with _StateGuard():
        main.load_inputs(MOL, MODE, DATA_DIR)  # fresh cache now exists

        call_count = {"n": 0}
        real_parse = GaussianParser.parse

        def _counting_parse(self, *a, **kw):
            call_count["n"] += 1
            return real_parse(self, *a, **kw)

        with mock.patch.object(GaussianParser, "parse", _counting_parse):
            main.load_inputs(MOL, MODE, DATA_DIR, use_cache=False)

        assert call_count["n"] == 1, "use_cache=False should force a re-parse"


if __name__ == "__main__":
    # Standalone runner so the suite works even without pytest installed.
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    passed = 0
    for fn in tests:
        try:
            fn()
            print(f"PASS  {fn.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"FAIL  {fn.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            print(f"ERROR {fn.__name__}: {type(e).__name__}: {e}")
    print(f"\n{passed}/{len(tests)} passed")
    sys.exit(0 if passed == len(tests) else 1)
