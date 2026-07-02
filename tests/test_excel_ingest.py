"""Regression tests for src/excel_ingest.py (Phase 3 library ingest).

Run from ``Github/scoring-functions/``:
    py -m pytest tests/                (with pytest)
    py tests/test_excel_ingest.py     (standalone; no pytest needed)

Design note: opening the ~10 MB workbook with openpyxl takes on the order of
a minute, so these tests deliberately do NOT re-run the full
``build_library_scores()``/xlsx-reading pipeline (that has already been run
and independently spot-checked this session -- see
IMPLEMENTATION_PLAN.md's Phase-3 Changelog entry for the verification
numbers). Instead they check the already-committed, already-generated
``data/results/library_scores.csv`` artifact (fast, plain pandas), plus a
couple of pure-logic unit tests of ``resolve_log_basename`` that need no I/O
at all. This mirrors the repo's existing pattern of treating
``data/results/*.csv`` as checked-in goldens.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import pandas as pd                                                # noqa: E402

from src.excel_ingest import (                                     # noqa: E402
    resolve_log_basename, EXCLUDED_MOLECULES, attach_geometry_classification,
)
from src.classifier import is_clean_external                       # noqa: E402

LIB_CSV = os.path.join(ROOT, "data", "results", "library_scores.csv")


def _load():
    return pd.read_csv(LIB_CSV)


def test_resolve_log_basename_direct_and_mapped():
    """Unit test of the Excel-name -> on-disk-log-basename resolution (no I/O
    beyond os.path.exists checks against the real repo tree)."""
    data_dir = os.path.join(ROOT, "data")
    # Direct-match molecules (same name in Excel and on disk).
    assert resolve_log_basename("SnH4", data_dir) == "SnH4"
    assert resolve_log_basename("H2S", data_dir) == "H2S"
    # Explicitly name-mapped molecules.
    assert resolve_log_basename("H2O", data_dir) == "water"
    assert resolve_log_basename("C6H6", data_dir) == "benzene"
    assert resolve_log_basename("CO2", data_dir) == "co2_mp2_3-21g"
    assert resolve_log_basename("Cl2O", data_dir) == "ocl2"
    assert resolve_log_basename("OF2", data_dir) == "of2"
    assert resolve_log_basename("Br2O", data_dir) == "br2o"
    assert resolve_log_basename("SeBr2", data_dir) == "SeBr2-cc"
    # An Excel-only molecule with no .log/.gjf pair in this repo.
    assert resolve_log_basename("BH3", data_dir) is None
    assert resolve_log_basename("NotAMolecule", data_dir) is None


def test_library_excludes_gramicidin_fragment():
    """Gly5 (Decision 5, gramicidin fragment) must not appear at all."""
    df = _load()
    assert not (df["molecule"] == "Gly5").any()
    assert not df["molecule"].isin(EXCLUDED_MOLECULES).any()


def test_library_covers_tab_ideal_and_tab_nonideal_molecules():
    """All 11 tab:ideal entries and the tab:nonideal bent-AB2 + two-center
    exemplars are present (checked exhaustively against the .tex this
    session; see excel_ingest.py's module docstring for the 5-bromide gap
    that is the only real omission)."""
    df = _load()
    present = set(df["molecule"].unique())
    tab_ideal = {"SnO2", "TeH2", "InH3", "SbH3", "IH3", "SnH4", "XeH4",
                 "TeH4", "SbH5", "XeOH4", "TeH6"}
    assert tab_ideal <= present, tab_ideal - present

    tab_nonideal_bent_ab2 = {"OF2", "Cl2O", "Br2O", "H2S", "SF2", "SCl2",
                              "SBr2", "H2Se", "SeF2", "SeCl2", "SeBr2", "H2O2"}
    assert tab_nonideal_bent_ab2 <= present, tab_nonideal_bent_ab2 - present


def test_geometry_backed_molecules_have_external_rows():
    """The 25 molecules with a real .log/.gjf pair get n_T+n_R appended
    'external' rows (never present in the raw Excel sheet)."""
    df = _load()
    ext = df[df["kind"] == "external"]
    counts = ext.groupby("molecule").size()
    # Linear molecules (n_R=2) get 5 external rows; everything else gets 6.
    linear = {"CO2", "CS2", "CSe2", "CTe2"}
    for mol, n in counts.items():
        expected = 5 if mol in linear else 6
        assert n == expected, f"{mol}: expected {expected} external rows, got {n}"
    # 22 hydride-library molecules with logs + water/benzene/CO2 = 25.
    assert len(counts) == 25, sorted(counts.index)
    assert {"SnH4", "H2S", "C6H6", "H2O", "CO2"} <= set(counts.index)


def test_external_rows_classify_clean_translation_rotation():
    """Ideal T/R references (Eckart-Sayvetz, exact for normal modes) must
    classify clean (bare "Tx".."Rz", no trailing "*") for every geometry-backed
    molecule -- this is the completeness guarantee, not an anecdotal check."""
    df = _load()
    ext = df[df["kind"] == "external"]
    bad = ext[~ext["predicted_label"].apply(is_clean_external)]
    assert len(bad) == 0, bad[["molecule", "mode_index", "predicted_label"]]


def test_water_o_series_internal_rows_correctly_unmerged():
    """H2O/OF2/Cl2O/Br2O's Excel rows do not match this repo's own logs
    (discovered this session -- a real data-provenance mismatch, not a bug);
    their internal rows must be left has_geometry=False rather than
    force-merged, while their external (T/R) rows -- independent of the
    vibrational-frequency mismatch -- are still attached normally."""
    df = _load()
    for mol in ("H2O", "OF2", "Cl2O", "Br2O"):
        internal = df[(df["molecule"] == mol) & (df["kind"] == "internal")]
        assert len(internal) > 0
        assert not internal["has_geometry"].any(), mol
        external = df[(df["molecule"] == mol) & (df["kind"] == "external")]
        assert external["has_geometry"].all(), mol
        assert external["predicted_label"].apply(is_clean_external).all(), mol


def test_bond_scores_sum_to_v_stretch():
    """Per-bond s_AB (parsed out of the semicolon-joined string) sums back to
    V_Stretch for every internal row that has bond detail (eq:bondscore)."""
    df = _load()
    checked = 0
    for _, row in df[df["kind"] == "internal"].iterrows():
        if not isinstance(row["s_AB"], str) or not row["s_AB"]:
            continue
        total = sum(float(part.split(":")[1]) for part in row["s_AB"].split(";"))
        assert abs(total - row["V_Stretch"]) < 1e-3, (row["molecule"], row["mode_index"])
        checked += 1
    assert checked > 100, f"only checked {checked} rows -- unexpectedly few"


def test_unresolved_engine_mode_is_reported_not_silently_skipped():
    """formula-auditor finding (2026-07-02): a "Vib i" row that resolves to
    NO engine mode at all (m is None) must count as a mismatch -- not be
    silently `continue`d past -- so the molecule still trips the
    all-or-nothing gate and shows up in the skip report. Synthesize a
    molecule ('H2O' -> data/logs/water.log, a real geometry-backed molecule)
    with one bogus, unresolvable mode_index (999; water only has 3
    vibrational modes) alongside a mismatch that SHOULD have triggered a
    warning even before this fix, to isolate the new code path."""
    bogus_row = {
        "molecule": "H2O", "mode_index": 999, "kind": "internal",
        "freq": 12345.0, "ref_label": "bend", "ideal": "no",
        "V_Stretch": 0.5, "delta_b_mean": None, "s_AB": "", "rel_db": "",
        "has_geometry": False, "predicted_label": None,
        "predicted_annotation": None,
        "Tx": None, "Ty": None, "Tz": None, "Rx": None, "Ry": None, "Rz": None,
    }
    df_in = pd.DataFrame([bogus_row])
    data_dir = os.path.join(ROOT, "data")
    df_out, skip_report = attach_geometry_classification(df_in, data_dir)

    assert len(skip_report) == 1
    entry = skip_report[0]
    assert entry["molecule"] == "H2O"
    assert entry["n_mismatched"] == 1
    mode_index, engine_freq, excel_freq = entry["example"]
    assert mode_index == 999
    assert engine_freq is None          # the m-is-None sentinel, not skipped
    assert excel_freq == 12345.0

    # The bogus internal row must NOT have been half-merged.
    internal_out = df_out[df_out["kind"] == "internal"]
    assert len(internal_out) == 1
    assert not bool(internal_out.iloc[0]["has_geometry"])
    assert internal_out.iloc[0]["predicted_label"] is None

    # External (T/R) rows are independent of this gate and are still
    # appended normally for this geometry-backed molecule.
    external_out = df_out[df_out["kind"] == "external"]
    assert len(external_out) == 6
    assert external_out["has_geometry"].all()


def test_ideal_stretch_bend_populations_do_not_overlap():
    """Precondition src/calibrate.py's threshold derivation relies on: the
    ideal-molecule stretch/bend V_Stretch populations are cleanly separated
    (verified this session: gap width ~0.73)."""
    df = _load()
    ideal = df[(df["kind"] == "internal") & (df["ideal"] == "yes")]
    stretch_min = ideal.loc[ideal["ref_label"] == "stretch", "V_Stretch"].min()
    bend_max = ideal.loc[ideal["ref_label"] == "bend", "V_Stretch"].max()
    assert bend_max < stretch_min
    assert stretch_min > 0.9
    assert bend_max < 0.2


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
