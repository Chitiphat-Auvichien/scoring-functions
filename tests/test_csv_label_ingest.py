"""Regression tests for src/csv_label_ingest.py -- the 2026-07-05 default
(sole source since the 2026-07-07 xlsx-ingest retirement) label/citation-key
source: data/characterised_modes.csv, data/ref-label_citation.csv (see the
module docstring). The 6 molecules these CSVs do not cover (C2H2/C2H4/C2H6/
H2O2/iso-C4H10/n-C4H10) are outside data/mol_list_method.csv's 72-molecule
roster entirely -- genuinely out of scope, not a gap needing a fallback (the
old xlsx-fallback mechanism and its tests are removed).

**2026-07-08 (data_score.csv retirement):** ``data/data_score.csv`` is no
longer read by this module at all -- ``load_label_csvs()`` returns only
``"characterised_modes"``/``"citations"``, and ``build_label_lookup()`` /
``get_label()`` no longer supply an ``ideal`` value (that is now sourced
structurally from ``data/mol_list_method.csv``'s ``mol_type`` column by
``src.library_ingest.attach_ideal_tags()``, entirely outside this module --
see ``tests/test_library_ingest.py``). ``get_label()`` now returns a
2-tuple, ``(ref_label, ref_key)``.

Run from ``Github/scoring-functions/``:
    py -m pytest tests/                          (with pytest)
    py tests/test_csv_label_ingest.py           (standalone; no pytest needed)
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import pandas as pd                                                # noqa: E402

from src.csv_label_ingest import (                                 # noqa: E402
    load_label_csvs, build_label_lookup, get_label, build_citation_table,
    ALLOWED_REF_LABELS, LABEL_CSV_FILES,
)

DATA_DIR = os.path.join(ROOT, "data")


def test_label_csv_files_no_longer_includes_data_score():
    """2026-07-08: data_score.csv is retired from this module's read path
    entirely."""
    assert set(LABEL_CSV_FILES.keys()) == {"characterised_modes", "citations"}
    assert "data_score" not in LABEL_CSV_FILES


def test_load_label_csvs_reads_both_with_expected_columns():
    tables = load_label_csvs(DATA_DIR)
    assert set(tables.keys()) == {"characterised_modes", "citations"}
    assert {"molecule", "mode", "type", "ref"}.issubset(tables["characterised_modes"].columns)
    assert {"code", "doi"}.issubset(tables["citations"].columns)


def test_allowed_ref_labels_includes_the_new_sb_class():
    assert set(ALLOWED_REF_LABELS) == {"bend", "stretch", "SB"}


def test_build_label_lookup_accepts_sb_and_rejects_garbage():
    tables = {
        "characterised_modes": pd.DataFrame([
            {"molecule": "MOLA", "mode": 1, "type": "bend", "ref": None},
            {"molecule": "MOLA", "mode": 2, "type": "SB", "ref": "Bar1999"},
            {"molecule": "MOLA", "mode": 3, "type": "garbage", "ref": None},
        ]),
    }
    lookup = build_label_lookup(tables)
    assert get_label(lookup, "MOLA", 1) == ("bend", None)
    assert get_label(lookup, "MOLA", 2) == ("SB", "Bar1999")
    assert get_label(lookup, "MOLA", 3)[0] is None  # unrecognized type -> None, not silently kept


def test_get_label_missing_key_returns_none_pair():
    lookup = build_label_lookup({
        "characterised_modes": pd.DataFrame(columns=["molecule", "mode", "type", "ref"]),
    })
    assert get_label(lookup, "NOPE", 1) == (None, None)


def test_build_label_lookup_does_not_return_an_ideal_key():
    """'ideal' is no longer part of this module's output at all -- confirms
    the retirement is complete, not just untested."""
    tables = {
        "characterised_modes": pd.DataFrame([
            {"molecule": "MOLB", "mode": 1, "type": "stretch", "ref": None},
        ]),
    }
    lookup = build_label_lookup(tables)
    assert set(lookup[("MOLB", 1)].keys()) == {"ref_label", "ref_key"}


def test_build_citation_table_skips_blank_codes_and_reads_real_file():
    tables = load_label_csvs(DATA_DIR)
    citations = build_citation_table(tables)
    assert "" not in citations
    assert "Shi1972" in citations
    assert citations["Shi1972"]["doi"] is not None or citations["Shi1972"]["author"] is not None


if __name__ == "__main__":
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
