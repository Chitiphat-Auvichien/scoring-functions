"""Regression tests for src/csv_label_ingest.py -- the 2026-07-05 default
label/citation-key source (data/data_score.csv, data/characterised_modes.csv,
data/ref-label_citation.csv), superseding the xlsx workbook's data_score/
characterised modes sheets for ref_label/ideal/citation content across ALL
molecules (with an xlsx fallback for the handful of molecules the new CSVs
don't cover yet -- see the module docstring).

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
    ALLOWED_REF_LABELS,
)

DATA_DIR = os.path.join(ROOT, "data")


def test_load_label_csvs_reads_all_three_with_expected_columns():
    tables = load_label_csvs(DATA_DIR)
    assert set(tables.keys()) == {"data_score", "characterised_modes", "citations"}
    assert {"molecule", "mode", "type", "ideal"}.issubset(tables["data_score"].columns)
    assert {"molecule", "mode", "ref"}.issubset(tables["characterised_modes"].columns)
    assert {"code", "doi"}.issubset(tables["citations"].columns)


def test_allowed_ref_labels_includes_the_new_sb_class():
    assert set(ALLOWED_REF_LABELS) == {"bend", "stretch", "SB"}


def test_build_label_lookup_accepts_sb_and_rejects_garbage():
    tables = {
        "data_score": pd.DataFrame([
            {"molecule": "MOLA", "mode": 1, "type": "bend", "ideal": "no"},
            {"molecule": "MOLA", "mode": 2, "type": "SB", "ideal": "no"},
            {"molecule": "MOLA", "mode": 3, "type": "garbage", "ideal": "no"},
        ]),
        "characterised_modes": pd.DataFrame(columns=["molecule", "mode", "ref"]),
    }
    lookup = build_label_lookup(tables)
    assert get_label(lookup, "MOLA", 1)[0] == "bend"
    assert get_label(lookup, "MOLA", 2)[0] == "SB"
    assert get_label(lookup, "MOLA", 3)[0] is None  # unrecognized type -> None, not silently kept


def test_build_label_lookup_new_csv_supersedes_fallback_when_molecule_is_covered():
    """A molecule present in the new CSVs always wins over the xlsx
    fallback, even if the fallback disagrees -- the 'full switch' contract."""
    new_tables = {
        "data_score": pd.DataFrame([
            {"molecule": "MOLA", "mode": 1, "type": "stretch", "ideal": "yes"},
        ]),
        "characterised_modes": pd.DataFrame([
            {"molecule": "MOLA", "mode": 1, "ref": "NewKey2026"},
        ]),
    }
    fallback_ds = pd.DataFrame([
        {"molecule": "MOLA", "mode": 1, "type": "bend", "ideal": "no"},  # disagrees -- must lose
    ])
    lookup = build_label_lookup(new_tables, fallback_ds=fallback_ds)
    ref_label, ideal, ref_key = get_label(lookup, "MOLA", 1)
    assert ref_label == "stretch"
    assert ideal == "yes"
    assert ref_key == "NewKey2026"


def test_build_label_lookup_falls_back_for_molecule_new_csvs_do_not_cover():
    """A molecule ABSENT from the new CSVs entirely (e.g. the 6 real
    molecules this applies to today: C2H2/C2H4/C2H6/H2O2/iso-C4H10/
    n-C4H10) must reproduce the xlsx fallback's own values exactly, with no
    ref_key (the new CSVs supply citation keys; the xlsx fallback never
    did)."""
    new_tables = {
        "data_score": pd.DataFrame([
            {"molecule": "MOLA", "mode": 1, "type": "stretch", "ideal": "yes"},
        ]),
        "characterised_modes": pd.DataFrame(columns=["molecule", "mode", "ref"]),
    }
    fallback_ds = pd.DataFrame([
        {"molecule": "MOLA", "mode": 1, "type": "bend", "ideal": "no"},
        {"molecule": "MOLB", "mode": 1, "type": "bend", "ideal": "no"},  # not in new CSVs at all
    ])
    lookup = build_label_lookup(new_tables, fallback_ds=fallback_ds)

    # MOLA: covered by new CSVs -> new values win.
    assert get_label(lookup, "MOLA", 1) == ("stretch", "yes", None)
    # MOLB: uncovered -> xlsx fallback values, ref_key always None.
    assert get_label(lookup, "MOLB", 1) == ("bend", "no", None)


def test_build_label_lookup_fallback_never_recognizes_sb():
    """The xlsx fallback path predates the SB class entirely -- a fallback
    row with type=='SB' (hypothetically) would map to None, since the old
    xlsx sheet's own filter (bend/stretch only) is preserved verbatim for
    the fallback so uncovered molecules reproduce EXACTLY their pre-
    2026-07-05 values, never a value the old pipeline could not have
    produced."""
    new_tables = {
        "data_score": pd.DataFrame(columns=["molecule", "mode", "type", "ideal"]),
        "characterised_modes": pd.DataFrame(columns=["molecule", "mode", "ref"]),
    }
    fallback_ds = pd.DataFrame([
        {"molecule": "MOLC", "mode": 1, "type": "SB", "ideal": "no"},
    ])
    lookup = build_label_lookup(new_tables, fallback_ds=fallback_ds)
    assert get_label(lookup, "MOLC", 1)[0] is None


def test_get_label_missing_key_returns_none_triple():
    lookup = build_label_lookup({
        "data_score": pd.DataFrame(columns=["molecule", "mode", "type", "ideal"]),
        "characterised_modes": pd.DataFrame(columns=["molecule", "mode", "ref"]),
    })
    assert get_label(lookup, "NOPE", 1) == (None, None, None)


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
