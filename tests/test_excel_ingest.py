"""Regression tests for src/excel_ingest.py (Phase 3 library ingest).

Run from ``Github/scoring-functions/``:
    py -m pytest tests/                (with pytest)
    py tests/test_excel_ingest.py     (standalone; no pytest needed)

**Updated 2026-07-04 for the dual-source architecture** (see excel_ingest.py's
module docstring + IMPLEMENTATION_PLAN.md's RESUME HERE): ``source="excel"``
was restored as the default (interim decision, reproducing the population/
values the manuscript's currently-typeset figures were built from) alongside
the still-fully-working ``source="gaussian"`` disk-driven path from
2026-07-03. The checked-in ``data/results/library_scores.csv`` golden this
file reads is now EXCEL-sourced (~69 molecules / 659 rows, ``has_geometry``
False for the ~45 Excel-only molecules, True for the ~25 that also have an
on-disk ``.log``/``.gjf`` pair) -- NOT the disk-driven ~25-molecule population
this file's tests were written against on 2026-07-03. Tests that check
excel-sourced-golden-specific facts read ``_load()`` (the checked-in CSV)
directly; tests that check source="gaussian"-specific behavioral guarantees
(every row geometry-backed, molecule count == disk roster) instead build a
small in-memory ``source="gaussian"`` DataFrame via ``_gaussian_df()``
(memoized per test session) so they stay valid regardless of which source
currently produced the checked-in golden.

Design note: opening the ~10 MB workbook with openpyxl takes on the order of a
minute, so these tests avoid re-running the full ``source="excel"``
``build_library_scores()`` pipeline (already run and spot-checked this
session -- see IMPLEMENTATION_PLAN.md's RESUME HERE for the verification
numbers); the cheaper ``source="gaussian"`` build (one Excel sheet + 25 small
Gaussian logs, no ``data_mode&bond`` read) is used instead for the handful of
tests that need a live disk-driven DataFrame. Otherwise these tests check:
(a) the already-committed, already-generated ``data/results/library_scores.csv``
artifact (fast, plain pandas) as a checked-in golden, (b) pure-logic/no-I/O
unit tests of ``resolve_log_basename``/``resolve_excel_molecule_name``,
(c) the disk listing ``discover_geometry_molecules`` (fast, just os.listdir),
(d) a direct, fast call to ``score_geometry_molecule`` on water (9 modes,
sub-second) to validate the real-engine row-building logic independent of the
CSV golden, and (e) a fully synthetic, no-I/O unit test of
``attach_excel_labels``'s frequency-gating logic (the one piece of ingest
logic most likely to silently regress).
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import pandas as pd                                                # noqa: E402

from src.excel_ingest import (                                     # noqa: E402
    resolve_log_basename, resolve_excel_molecule_name,
    discover_geometry_molecules, score_geometry_molecule,
    attach_excel_labels, build_library_scores,
    EXCLUDED_MOLECULES, EXCEL_TO_LOG, _EXTERNAL_SLOTS,
)
from src.classifier import is_clean_external                       # noqa: E402

LIB_CSV = os.path.join(ROOT, "data", "results", "library_scores.csv")
DATA_DIR = os.path.join(ROOT, "data")


def _load():
    return pd.read_csv(LIB_CSV)


_gaussian_df_cache = {}


def _gaussian_df():
    """Memoized source="gaussian" build (disk-driven, ~25 molecules) --
    independent of whatever source produced the checked-in library_scores.csv
    golden. Computed once per test session (module-level dict cache), not
    once per test, to keep the suite fast."""
    if "df" not in _gaussian_df_cache:
        _gaussian_df_cache["df"] = build_library_scores(data_dir=DATA_DIR, source="gaussian")
    return _gaussian_df_cache["df"]


def test_resolve_log_basename_direct_and_mapped():
    """Unit test of the Excel-name -> on-disk-log-basename resolution (no I/O
    beyond os.path.exists checks against the real repo tree). Unaffected by
    the disk-driven rearchitecture -- this function is unchanged."""
    # Direct-match molecules (same name in Excel and on disk).
    assert resolve_log_basename("SnH4", DATA_DIR) == "SnH4"
    assert resolve_log_basename("H2S", DATA_DIR) == "H2S"
    # Explicitly name-mapped molecules.
    assert resolve_log_basename("H2O", DATA_DIR) == "water"
    assert resolve_log_basename("C6H6", DATA_DIR) == "benzene"
    assert resolve_log_basename("CO2", DATA_DIR) == "co2_mp2_3-21g"
    assert resolve_log_basename("Cl2O", DATA_DIR) == "ocl2"
    assert resolve_log_basename("OF2", DATA_DIR) == "of2"
    assert resolve_log_basename("Br2O", DATA_DIR) == "br2o"
    assert resolve_log_basename("SeBr2", DATA_DIR) == "SeBr2-cc"
    # An Excel-only molecule with no .log/.gjf pair in this repo.
    assert resolve_log_basename("BH3", DATA_DIR) is None
    assert resolve_log_basename("NotAMolecule", DATA_DIR) is None


def test_resolve_excel_molecule_name_mapped_direct_and_absent():
    """The reverse direction (on-disk basename -> Excel molecule name) that
    build_library_scores() actually uses to label rows. Mapped names take
    priority over a direct match; a basename with no Excel counterpart at all
    (not in _LOG_TO_EXCEL, not itself an Excel molecule name) returns None,
    which build_library_scores() then falls back to using the basename itself
    -- not an error (module docstring point 5)."""
    excel_molecules = {"SnH4", "H2S", "TeH6"}
    # Explicitly mapped (EXCEL_TO_LOG reversed).
    assert resolve_excel_molecule_name("water", excel_molecules) == "H2O"
    assert resolve_excel_molecule_name("benzene", excel_molecules) == "C6H6"
    assert resolve_excel_molecule_name("co2_mp2_3-21g", excel_molecules) == "CO2"
    assert resolve_excel_molecule_name("ocl2", excel_molecules) == "Cl2O"
    assert resolve_excel_molecule_name("of2", excel_molecules) == "OF2"
    assert resolve_excel_molecule_name("br2o", excel_molecules) == "Br2O"
    assert resolve_excel_molecule_name("SeBr2-cc", excel_molecules) == "SeBr2"
    # Direct match (basename itself is an Excel molecule name).
    assert resolve_excel_molecule_name("SnH4", excel_molecules) == "SnH4"
    assert resolve_excel_molecule_name("TeH6", excel_molecules) == "TeH6"
    # No Excel counterpart at all.
    assert resolve_excel_molecule_name("totally_unmapped_xyz", excel_molecules) is None


def test_every_excel_to_log_target_is_present_in_reverse_map():
    """EXCEL_TO_LOG's forward map must be exactly invertible (no duplicate
    on-disk targets silently shadowing each other) -- a defensive check on the
    fixed name-mapping table itself."""
    targets = list(EXCEL_TO_LOG.values())
    assert len(targets) == len(set(targets)), "EXCEL_TO_LOG has duplicate on-disk targets"


def test_discover_geometry_molecules_is_disk_driven_and_scales():
    """The disk listing that drives build_library_scores()'s whole molecule
    loop -- must reflect exactly what's in data/logs/ + data/gjf/ right now,
    with no hardcoded roster. Currently 25 molecules (grows automatically as
    more .log/.gjf pairs are dropped in; no code change needed -- see module
    docstring)."""
    bases = discover_geometry_molecules(DATA_DIR)
    assert bases == sorted(bases)  # sorted, deterministic order
    assert len(bases) == len(set(bases))  # no duplicates
    for expected in ("water", "benzene", "co2_mp2_3-21g", "SnH4", "H2S"):
        assert expected in bases, expected
    assert len(bases) == 25


def test_library_row_count_and_molecule_count_match_disk_roster():
    """source="gaussian"'s molecule count must equal discover_geometry_
    molecules()'s count exactly (one row-group per on-disk molecule, no
    Excel-only molecules sneaking in and no on-disk molecule silently
    dropped). Uses the in-memory source="gaussian" build (_gaussian_df), NOT
    the checked-in golden -- that golden is source="excel" by default as of
    2026-07-04 and legitimately has ~69 molecules (see
    test_excel_sourced_default_has_full_population below)."""
    df = _gaussian_df()
    bases = discover_geometry_molecules(DATA_DIR)
    assert df["molecule"].nunique() == len(bases) == 25


def test_library_excludes_gramicidin_fragment():
    """Gly5 (Decision 5, gramicidin fragment) must not appear at all."""
    df = _load()
    assert not (df["molecule"] == "Gly5").any()
    assert not df["molecule"].isin(EXCLUDED_MOLECULES).any()


def test_all_rows_are_geometry_backed():
    """Under source="gaussian", EVERY row comes from a real .log/.gjf pair --
    has_geometry is unconditionally True. Uses _gaussian_df(), not the
    checked-in golden (which is source="excel" by default as of 2026-07-04
    and has has_geometry==False for its ~45 Excel-only molecules -- see
    test_excel_sourced_default_has_full_population)."""
    df = _gaussian_df()
    assert df["has_geometry"].all()
    assert len(df) > 0


def test_geometry_backed_molecules_have_expected_external_row_counts():
    """Every molecule gets n_T+n_R appended 'external' rows: 5 for linear
    molecules (n_R=2), 6 otherwise."""
    df = _load()
    ext = df[df["kind"] == "external"]
    counts = ext.groupby("molecule").size()
    linear = {"CO2", "CS2", "CSe2", "CTe2"}
    for mol, n in counts.items():
        expected = 5 if mol in linear else 6
        assert n == expected, f"{mol}: expected {expected} external rows, got {n}"
    assert len(counts) == 25, sorted(counts.index)
    assert {"SnH4", "H2S", "C6H6", "H2O", "CO2"} <= set(counts.index)


def test_external_rows_classify_clean_translation_rotation():
    """Ideal T/R references (Eckart-Sayvetz, exact for normal modes) must
    classify clean (bare "Tx".."Rz", no trailing "*") for every geometry-backed
    molecule -- the completeness guarantee, not an anecdotal check."""
    df = _load()
    ext = df[df["kind"] == "external"]
    bad = ext[~ext["predicted_label"].apply(is_clean_external)]
    assert len(bad) == 0, bad[["molecule", "mode_index", "predicted_label"]]
    # And ref_label is the unconditional structural truth, never null, for
    # every external row (independent of any Excel join).
    assert ext["ref_label"].isin(("translation", "rotation")).all()


def test_frequency_mismatched_molecules_leave_internal_label_null_only_gaussian():
    """(source="gaussian") H2O/OF2/Cl2O/Br2O's Excel rows do not match this
    repo's own logs (discovered in a prior session -- a real data-provenance
    mismatch, not a bug). Under source="gaussian" their internal rows are
    still fully SCORED (has_geometry True, V_Stretch/predicted_label
    populated) -- only ref_label/ideal are left null because the frequency-
    gated label join is skipped. External (T/R) rows are completely
    unaffected. Uses _gaussian_df(), not the checked-in (source="excel")
    golden -- see test_excel_sourced_frequency_mismatched_molecules_keep_
    ref_label_without_geometry_overlay for that source's differently-shaped
    (but equally correct) contract for the same 4 molecules."""
    df = _gaussian_df()
    for mol in ("H2O", "OF2", "Cl2O", "Br2O"):
        internal = df[(df["molecule"] == mol) & (df["kind"] == "internal")]
        assert len(internal) > 0, mol
        assert internal["ref_label"].isna().all(), mol
        assert internal["ideal"].isna().all(), mol
        # Scores are NOT affected by the label-join skip.
        assert internal["has_geometry"].all(), mol
        assert internal["V_Stretch"].notna().all(), mol
        assert internal["predicted_label"].notna().all(), mol

        external = df[(df["molecule"] == mol) & (df["kind"] == "external")]
        assert external["has_geometry"].all(), mol
        assert external["ref_label"].isin(("translation", "rotation")).all(), mol
        assert external["predicted_label"].apply(is_clean_external).all(), mol


def test_excel_sourced_default_has_full_population():
    """The checked-in library_scores.csv golden is source="excel" by default
    (2026-07-04 interim decision -- see excel_ingest.py's module docstring):
    the full ~69-molecule hydride library, not just the ~25 with on-disk
    geometry. Matches the pre-2026-07-03 (commit 149fc62) population this
    repo's manuscript figures were originally built from."""
    df = _load()
    assert df["molecule"].nunique() == 69
    assert len(df) == 659
    # Both geometry-backed (on-disk .log/.gjf) and Excel-only rows coexist.
    assert df["has_geometry"].any()
    assert not df["has_geometry"].all()


def test_excel_sourced_frequency_mismatched_molecules_keep_ref_label_without_geometry_overlay():
    """(source="excel") H2O/OF2/Cl2O/Br2O's internal rows get their
    V_Stretch/ref_label/ideal DIRECTLY from Excel (never gated by the
    frequency check -- that check only gates the geometry OVERLAY of
    predicted_label/Tx..Rz), so ref_label stays populated even though
    has_geometry is False for these 4 (their on-disk .log frequencies don't
    match Excel's, so the overlay is correctly skipped). This is the
    source="excel" mirror of test_frequency_mismatched_molecules_leave_
    internal_label_null_only_gaussian -- same 4 molecules, differently-shaped
    but equally correct contract under the other source."""
    df = _load()
    for mol in ("H2O", "OF2", "Cl2O", "Br2O"):
        internal = df[(df["molecule"] == mol) & (df["kind"] == "internal")]
        assert len(internal) > 0, mol
        assert internal["ref_label"].isin(("bend", "stretch")).all(), mol
        assert internal["V_Stretch"].notna().all(), mol
        # Geometry overlay skipped (frequency mismatch) -- no engine
        # predicted_label/Tx..Rz for these internal rows.
        assert not internal["has_geometry"].any(), mol
        assert internal["predicted_label"].isna().all(), mol

        # External (ideal T/R) rows are unaffected -- geometry-only, no
        # frequency dependency at all.
        external = df[(df["molecule"] == mol) & (df["kind"] == "external")]
        assert external["has_geometry"].all(), mol
        assert external["ref_label"].isin(("translation", "rotation")).all(), mol
        assert external["predicted_label"].apply(is_clean_external).all(), mol


def test_bond_scores_sum_to_v_stretch():
    """Per-bond s_AB (parsed out of the semicolon-joined string) sums back to
    V_Stretch for every internal row that has bond detail (eq:bondscore).
    Every internal row now carries bond detail regardless of label (the new
    module always calls score_bonds() directly, not classify_all_modes()'s
    filtered subset)."""
    df = _load()
    checked = 0
    for _, row in df[df["kind"] == "internal"].iterrows():
        assert isinstance(row["s_AB"], str) and row["s_AB"], \
            (row["molecule"], row["mode_index"])
        total = sum(float(part.split(":")[1]) for part in row["s_AB"].split(";"))
        assert abs(total - row["V_Stretch"]) < 1e-3, (row["molecule"], row["mode_index"])
        checked += 1
    assert checked > 100, f"only checked {checked} rows -- unexpectedly few"


def test_rel_db_string_present_and_parseable_for_every_internal_row():
    """rel_db is the diagnostic companion column (fig:bondscores' x-axis);
    every internal row must carry a parseable, same-bond-count string."""
    df = _load()
    for _, row in df[df["kind"] == "internal"].iterrows():
        assert isinstance(row["rel_db"], str) and row["rel_db"], \
            (row["molecule"], row["mode_index"])
        s_ab_bonds = row["s_AB"].split(";")
        rel_db_bonds = row["rel_db"].split(";")
        assert len(s_ab_bonds) == len(rel_db_bonds), (row["molecule"], row["mode_index"])
        for part in rel_db_bonds:
            float(part.split(":")[1])  # must parse cleanly


def test_score_geometry_molecule_water_direct():
    """Direct, fast (no xlsx I/O) call to score_geometry_molecule() on water
    -- validates the real-engine row-building logic itself, independent of
    the checked-in CSV golden. Water: 3N=9 -> 3 T + 3 R + 3 internal rows."""
    rows = score_geometry_molecule("water", DATA_DIR)
    assert len(rows) == 9
    external = [r for r in rows if r["kind"] == "external"]
    internal = [r for r in rows if r["kind"] == "internal"]
    assert len(external) == 6  # water is non-linear: n_T=3, n_R=3
    assert len(internal) == 3
    assert {r["mode_index"] for r in external} == set(_EXTERNAL_SLOTS)
    assert sorted(r["mode_index"] for r in internal) == [1, 2, 3]

    for r in external:
        assert is_clean_external(r["predicted_label"])
        assert r["ref_label"] in ("translation", "rotation")
        assert r["ideal"] is None  # label-only join not yet attached
        assert r["has_geometry"] is True

    for r in internal:
        assert r["ref_label"] is None and r["ideal"] is None  # attached later
        assert r["has_geometry"] is True
        assert 0.0 <= r["V_Stretch"] <= 1.0
        assert isinstance(r["s_AB"], str) and r["s_AB"]
        assert r["delta_b_mean"] is not None and r["delta_b_mean"] >= 0.0
        # Per-bond s_AB sums to V_Stretch for this mode too.
        total = sum(float(part.split(":")[1]) for part in r["s_AB"].split(";"))
        assert abs(total - r["V_Stretch"]) < 1e-3


def test_attach_excel_labels_gates_on_frequency_and_leaves_unlisted_molecules_alone():
    """Fully synthetic (no xlsx I/O) unit test of attach_excel_labels()'s
    core contract: a matching, frequency-consistent internal row gets
    ref_label/ideal filled in; a mismatched (or Excel-absent) mode disqualifies
    the WHOLE molecule's internal-row label join (all-or-nothing gate); a
    molecule with no Excel counterpart at all is untouched (not an error);
    external rows are never touched."""
    tables = {"data_score": pd.DataFrame([
        {"molecule": "TESTMOL", "mode": 1, "freq": 1000.0, "type": "stretch", "ideal": "yes"},
        {"molecule": "TESTMOL", "mode": 2, "freq": 500.0, "type": "bend", "ideal": "no"},
        {"molecule": "TESTMOL2", "mode": 1, "freq": 1000.0, "type": "stretch", "ideal": "yes"},
    ])}

    df = pd.DataFrame([
        # TESTMOL: both internal modes match within tolerance -> fully joined.
        {"molecule": "TESTMOL", "mode_index": 1, "kind": "internal",
         "freq": 1000.0001, "ref_label": None, "ideal": None},
        {"molecule": "TESTMOL", "mode_index": 2, "kind": "internal",
         "freq": 500.0, "ref_label": None, "ideal": None},
        {"molecule": "TESTMOL", "mode_index": "Tx", "kind": "external",
         "freq": 0.0, "ref_label": "translation", "ideal": None},
        # TESTMOL2: mode 1 is off by 1.0 cm-1, outside atol=0.05/rtol=1e-4 ->
        # whole molecule's internal join skipped.
        {"molecule": "TESTMOL2", "mode_index": 1, "kind": "internal",
         "freq": 999.0, "ref_label": None, "ideal": None},
        # NOTINEXCEL: no Excel counterpart at all -> untouched, no warning.
        {"molecule": "NOTINEXCEL", "mode_index": 1, "kind": "internal",
         "freq": 42.0, "ref_label": None, "ideal": None},
    ])

    out, skip_report = attach_excel_labels(df, tables)

    tm = out[out["molecule"] == "TESTMOL"]
    assert tm.loc[tm["mode_index"] == 1, "ref_label"].iloc[0] == "stretch"
    assert tm.loc[tm["mode_index"] == 1, "ideal"].iloc[0] == "yes"
    assert tm.loc[tm["mode_index"] == 2, "ref_label"].iloc[0] == "bend"
    assert tm.loc[tm["mode_index"] == 2, "ideal"].iloc[0] == "no"
    # External row untouched (already correct before the call).
    assert tm.loc[tm["kind"] == "external", "ref_label"].iloc[0] == "translation"

    tm2 = out[out["molecule"] == "TESTMOL2"]
    assert tm2["ref_label"].isna().all()
    assert tm2["ideal"].isna().all()

    ne = out[out["molecule"] == "NOTINEXCEL"]
    assert ne["ref_label"].isna().all()

    assert len(skip_report) == 1
    assert skip_report[0]["molecule"] == "TESTMOL2"
    assert skip_report[0]["n_mismatched"] == 1
    mode_index, engine_freq, excel_freq = skip_report[0]["example"]
    assert mode_index == 1
    assert engine_freq == 999.0
    assert excel_freq == 1000.0


def test_attach_excel_labels_missing_excel_row_counts_as_mismatch():
    """A mode_index with no corresponding Excel row at all (not just a
    numeric mismatch) must also disqualify the molecule's join and be
    reported with excel_freq=None -- not silently skipped past."""
    tables = {"data_score": pd.DataFrame([
        {"molecule": "TESTMOL3", "mode": 1, "freq": 1000.0, "type": "stretch", "ideal": "yes"},
    ])}
    df = pd.DataFrame([
        {"molecule": "TESTMOL3", "mode_index": 1, "kind": "internal",
         "freq": 1000.0, "ref_label": None, "ideal": None},
        {"molecule": "TESTMOL3", "mode_index": 2, "kind": "internal",
         "freq": 2000.0, "ref_label": None, "ideal": None},  # no Excel row for mode 2
    ])
    out, skip_report = attach_excel_labels(df, tables)
    assert out["ref_label"].isna().all()  # all-or-nothing: mode 1's match doesn't survive
    assert len(skip_report) == 1
    mode_index, engine_freq, excel_freq = skip_report[0]["example"]
    assert excel_freq is None


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
