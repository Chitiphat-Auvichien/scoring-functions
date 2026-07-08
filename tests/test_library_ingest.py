"""Regression tests for src/library_ingest.py (formerly src/excel_ingest.py).

Run from ``Github/scoring-functions/``:
    py -m pytest tests/                (with pytest)
    py tests/test_library_ingest.py   (standalone; no pytest needed)

**Rewritten 2026-07-07 for the single, roster-driven pipeline** (final flip
to Gaussian-direct -- see IMPLEMENTATION_PLAN.md's RESUME HERE and
src/library_ingest.py's module docstring). The old dual-source
("source=excel"/"source=gaussian") architecture and its xlsx-workbook ingest
path are gone entirely; every one of ``data/mol_list_method.csv``'s 72
roster molecules now has a verified on-disk ``.log``/``.gjf`` pair, and
``library_scores.csv`` is built by iterating that roster directly. Tests
that used to check ``source="excel"``-specific behavior (the full-library
overlay, ``EXCEL_TO_LOG``/``resolve_excel_molecule_name``) are deleted
outright rather than adapted -- there is no equivalent concept left to test.

Design note: most tests here read the already-committed, already-regenerated
``data/results/library_scores.csv`` golden (fast, plain pandas) rather than
re-running the full ~70-molecule ``build_library_scores()`` (known-slow --
real Gaussian parsing + full Algorithm 1 classification for every roster
molecule). The handful of tests that need a live, uncached call either use
the fast single-molecule ``score_geometry_molecule()`` directly (water-sized,
sub-second) or monkeypatch ``load_mol_roster``/``discover_geometry_
molecules`` to restrict ``build_library_scores()`` to a tiny synthetic
population before calling it for real, so the fail-fast/warn-and-continue
control flow itself is exercised without the cost of a full run.
"""
import os
import sys
import tempfile
import warnings
from unittest.mock import patch

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import pandas as pd                                                # noqa: E402

from src.library_ingest import (                                   # noqa: E402
    load_mol_roster, check_roster_disk_consistency,
    resolve_log_basename, discover_geometry_molecules,
    score_geometry_molecule, attach_labels, build_library_scores,
    _EXTERNAL_SLOTS, SCHEMA_COLUMNS,
)
from src.classifier import is_clean_external                       # noqa: E402
from src.csv_label_ingest import build_label_lookup                # noqa: E402

LIB_CSV = os.path.join(ROOT, "data", "results", "library_scores.csv")
DATA_DIR = os.path.join(ROOT, "data")


def _load():
    return pd.read_csv(LIB_CSV)


# ---------------------------------------------------------------------------
# Roster (data/mol_list_method.csv) -- the authoritative (molecule, basename)
# registry.
# ---------------------------------------------------------------------------

def test_load_mol_roster_reads_all_72_with_basename_column():
    roster = load_mol_roster(DATA_DIR)
    assert {"molecule", "basename"}.issubset(roster.columns)
    assert len(roster) == 72
    assert roster["molecule"].is_unique
    assert roster["basename"].notna().all()


def test_load_mol_roster_raises_on_missing_required_column():
    with tempfile.TemporaryDirectory() as tmp_dir:
        with open(os.path.join(tmp_dir, "mol_list_method.csv"), "w") as f:
            f.write("molecule,shape\nH2O,bend\n")
        try:
            load_mol_roster(tmp_dir)
            assert False, "expected ValueError for a roster missing 'basename'"
        except ValueError as e:
            assert "basename" in str(e)


def test_resolve_log_basename_matches_roster_for_known_molecules():
    """resolve_log_basename() is now a pure roster lookup (no on-disk
    existence check) -- values must match data/mol_list_method.csv exactly."""
    assert resolve_log_basename("SbH3", DATA_DIR) == "SbH3_MP2_3-21G"
    assert resolve_log_basename("BrF3", DATA_DIR) == "brf3"
    assert resolve_log_basename("H2O", DATA_DIR) == "H2O-MP2-321G"
    assert resolve_log_basename("C6H6", DATA_DIR) == "C6H6_MP2_3-21G_D6h"
    assert resolve_log_basename("SeCl4", DATA_DIR) == "secl4_C2v"


def test_resolve_log_basename_returns_none_for_unknown_molecule():
    assert resolve_log_basename("NotAMolecule", DATA_DIR) is None


def test_check_roster_disk_consistency_real_roster_has_no_missing_or_orphaned():
    """The real, finalized roster: every basename resolves to an on-disk
    .log+.gjf pair, and no on-disk pair is unaccounted for (Phase A's own
    100%-coverage verification, re-checked here as a live regression test)."""
    roster = load_mol_roster(DATA_DIR)
    missing, orphaned = check_roster_disk_consistency(roster, DATA_DIR)
    assert missing == []
    assert orphaned == []


def test_check_roster_disk_consistency_detects_missing_row():
    roster = load_mol_roster(DATA_DIR)
    fake_row = pd.DataFrame([{"molecule": "FAKE", "basename": "not_a_real_basename_xyz"}])
    combined = pd.concat([roster, fake_row], ignore_index=True)
    missing, _ = check_roster_disk_consistency(combined, DATA_DIR)
    assert missing == [("FAKE", "not_a_real_basename_xyz")]


def test_check_roster_disk_consistency_detects_orphaned_basename():
    """Dropping a real roster row leaves its on-disk basename unreferenced --
    exactly the gramicidin-style 'file present, no roster row' case."""
    roster = load_mol_roster(DATA_DIR)
    dropped_basename = roster.loc[roster["molecule"] == "TeH2", "basename"].iloc[0]
    reduced = roster[roster["molecule"] != "TeH2"]
    _, orphaned = check_roster_disk_consistency(reduced, DATA_DIR)
    assert dropped_basename in orphaned


def test_discover_geometry_molecules_matches_roster_exactly():
    """The safety-net directory scan and the roster agree exactly (72
    basenames each) now that Phase A's coverage is complete."""
    bases = discover_geometry_molecules(DATA_DIR)
    assert bases == sorted(bases)
    assert len(bases) == len(set(bases))
    roster = load_mol_roster(DATA_DIR)
    assert set(bases) == set(roster["basename"])


# ---------------------------------------------------------------------------
# build_library_scores() control flow: fail-fast on a missing roster file,
# warn-and-continue on an orphaned disk file -- exercised on tiny synthetic
# populations so these stay fast.
# ---------------------------------------------------------------------------

def test_build_library_scores_raises_filenotfounderror_on_missing_basename():
    fake_roster = pd.DataFrame([{"molecule": "FAKE", "basename": "not_a_real_basename_xyz"}])
    with patch("src.library_ingest.load_mol_roster", return_value=fake_roster):
        try:
            build_library_scores(DATA_DIR)
            assert False, "expected FileNotFoundError"
        except FileNotFoundError as e:
            assert "FAKE" in str(e)
            assert "not_a_real_basename_xyz" in str(e)


def test_build_library_scores_warns_on_orphaned_disk_basename_and_still_builds():
    """Restrict the roster to just water (fast: 9 modes) while leaving the
    real 72-basename disk listing plus one extra fake basename -- missing
    stays empty (every roster row -- just water here -- resolves fine) but
    the fake basename (and every other real-but-unreferenced-by-this-
    reduced-roster basename) is reported as orphaned and warned about,
    non-fatally; the build still completes for the one roster row given."""
    full_roster = load_mol_roster(DATA_DIR)
    water_roster = full_roster[full_roster["molecule"] == "H2O"].reset_index(drop=True)
    real_bases = discover_geometry_molecules(DATA_DIR)
    fake_bases = real_bases + ["totally_fake_orphan_basename_xyz"]

    with patch("src.library_ingest.load_mol_roster", return_value=water_roster), \
         patch("src.library_ingest.discover_geometry_molecules", return_value=fake_bases):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            df = build_library_scores(DATA_DIR)

    orphan_msgs = [str(w.message) for w in caught
                   if "totally_fake_orphan_basename_xyz" in str(w.message)]
    assert len(orphan_msgs) == 1
    assert set(df["molecule"].unique()) == {"H2O"}
    assert df["has_geometry"].all()


# ---------------------------------------------------------------------------
# Checked-in golden (data/results/library_scores.csv).
# ---------------------------------------------------------------------------

def test_library_excludes_gramicidin_fragment():
    """Gly5 (Decision 5, gramicidin fragment) must not appear -- it is
    simply absent from mol_list_method.csv's 72-row roster now, not an
    explicit ingest-time exclusion list."""
    df = _load()
    assert not (df["molecule"] == "Gly5").any()


def test_full_population_has_geometry_for_every_row():
    """Every row in the checked-in golden is geometry-backed (has_geometry
    unconditionally True) -- the whole point of the roster-driven flip.
    Population may be 71 or 72 molecules depending on whether every roster
    row's files happen to parse/score cleanly (a molecule with files present
    but unusable connectivity is warned-and-skipped, not fabricated)."""
    df = _load()
    roster = load_mol_roster(DATA_DIR)
    assert df["has_geometry"].all()
    assert len(df) > 0
    assert set(df["molecule"].unique()) <= set(roster["molecule"])
    assert df["molecule"].nunique() >= len(roster) - 1  # at most 1 legitimate skip


def test_geometry_backed_molecules_have_expected_external_row_counts():
    """Every molecule gets n_T+n_R appended 'external' rows: 5 for linear
    molecules (n_R=2), 6 otherwise."""
    df = _load()
    ext = df[df["kind"] == "external"]
    counts = ext.groupby("molecule").size()
    linear = {"CO2", "CS2", "CSe2", "CTe2", "SnO2"}
    for mol, n in counts.items():
        expected = 5 if mol in linear else 6
        assert n == expected, f"{mol}: expected {expected} external rows, got {n}"
    assert {"SnH4", "H2S", "C6H6", "H2O", "CO2"} <= set(counts.index)


def test_external_rows_classify_clean_translation_rotation():
    """Ideal T/R references (Eckart-Sayvetz, exact for normal modes) must
    classify clean (bare "Tx".."Rz", no trailing "*") for every geometry-backed
    molecule -- the completeness guarantee, not an anecdotal check."""
    df = _load()
    ext = df[df["kind"] == "external"]
    bad = ext[~ext["predicted_label"].apply(is_clean_external)]
    assert len(bad) == 0, bad[["molecule", "mode_index", "predicted_label"]]
    assert ext["ref_label"].isin(("translation", "rotation")).all()


def test_bond_scores_sum_to_v_stretch():
    """Per-bond s_AB (parsed out of the semicolon-joined string) sums back to
    V_Stretch for every internal row (eq:bondscore). Every internal row
    carries bond detail regardless of label (score_geometry_molecule always
    calls score_bonds() directly, not classify_all_modes()'s filtered
    subset)."""
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


def test_ideal_stretch_bend_populations_do_not_overlap():
    """Precondition src/calibrate.py's threshold derivation relies on: the
    ideal-molecule stretch/bend V_Stretch populations are cleanly separated."""
    df = _load()
    ideal = df[(df["kind"] == "internal") & (df["ideal"] == "yes")]
    stretch_min = ideal.loc[ideal["ref_label"] == "stretch", "V_Stretch"].min()
    bend_max = ideal.loc[ideal["ref_label"] == "bend", "V_Stretch"].max()
    assert bend_max < stretch_min


# ---------------------------------------------------------------------------
# score_geometry_molecule() -- direct, fast, no-golden-CSV-needed check of
# the real-engine row-building logic itself.
# ---------------------------------------------------------------------------

def test_score_geometry_molecule_water_direct():
    """Water: 3N=9 -> 3 T + 3 R + 3 internal rows. Uses the roster's own
    on-disk basename for H2O, not the literal molecule name."""
    base = resolve_log_basename("H2O", DATA_DIR)
    rows = score_geometry_molecule(base, DATA_DIR)
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
        # construct_T/construct_R are synthetic -- never carry mu/k/irrep.
        assert r["reduced_mass"] is None
        assert r["force_constant"] is None
        assert r["irrep"] is None

    for r in internal:
        assert r["ref_label"] is None and r["ideal"] is None  # attached later
        assert r["has_geometry"] is True
        assert 0.0 <= r["V_Stretch"] <= 1.0
        assert isinstance(r["s_AB"], str) and r["s_AB"]
        assert r["delta_b_mean"] is not None and r["delta_b_mean"] >= 0.0
        total = sum(float(part.split(":")[1]) for part in r["s_AB"].split(";"))
        assert abs(total - r["V_Stretch"]) < 1e-3
        # 2026-07-08 Gaussian-direct parser rework: real internal modes carry
        # engine-parsed mu/k/irrep.
        assert r["reduced_mass"] is not None and r["reduced_mass"] > 0.0
        assert r["force_constant"] is not None and r["force_constant"] > 0.0
        assert r["irrep"] is not None and isinstance(r["irrep"], str)


def test_schema_columns_includes_mu_k_irrep():
    """reduced_mass/force_constant/irrep (2026-07-08) are appended to the end
    of the locked schema -- purely additive, existing columns untouched."""
    assert SCHEMA_COLUMNS[-3:] == ["reduced_mass", "force_constant", "irrep"]
    assert SCHEMA_COLUMNS[:-3] == [
        "molecule", "mode_index", "kind", "freq", "ref_label", "ideal",
        "V_Stretch", "delta_b_mean", "s_AB", "rel_db", "has_geometry",
        "predicted_label", "predicted_annotation",
        "Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "ref_key",
    ]


# ---------------------------------------------------------------------------
# attach_labels() -- fully synthetic (no CSV I/O) unit tests of the
# frequency-gating logic.
# ---------------------------------------------------------------------------

def test_attach_labels_gates_on_frequency_and_leaves_unlisted_molecules_alone():
    """A matching, frequency-consistent internal row gets ref_label/ideal
    filled in; a mismatched (or data_score.csv-absent) mode disqualifies the
    WHOLE molecule's internal-row label join (all-or-nothing gate); a
    molecule with no data_score.csv counterpart at all is untouched (not an
    error); external rows are never touched."""
    csv_tables = {"data_score": pd.DataFrame([
        {"molecule": "TESTMOL", "mode": 1, "freq": 1000.0, "type": "stretch", "ideal": "yes"},
        {"molecule": "TESTMOL", "mode": 2, "freq": 500.0, "type": "bend", "ideal": "no"},
        {"molecule": "TESTMOL2", "mode": 1, "freq": 1000.0, "type": "stretch", "ideal": "yes"},
    ]), "characterised_modes": pd.DataFrame(columns=["molecule", "mode", "ref"])}

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
        # NOTINCSV: no data_score.csv counterpart at all -> untouched, no warning.
        {"molecule": "NOTINCSV", "mode_index": 1, "kind": "internal",
         "freq": 42.0, "ref_label": None, "ideal": None},
    ])

    label_lookup = build_label_lookup(csv_tables)
    out, skip_report = attach_labels(df, csv_tables, label_lookup)

    tm = out[out["molecule"] == "TESTMOL"]
    assert tm.loc[tm["mode_index"] == 1, "ref_label"].iloc[0] == "stretch"
    assert tm.loc[tm["mode_index"] == 1, "ideal"].iloc[0] == "yes"
    assert tm.loc[tm["mode_index"] == 2, "ref_label"].iloc[0] == "bend"
    assert tm.loc[tm["mode_index"] == 2, "ideal"].iloc[0] == "no"
    assert tm.loc[tm["kind"] == "external", "ref_label"].iloc[0] == "translation"

    tm2 = out[out["molecule"] == "TESTMOL2"]
    assert tm2["ref_label"].isna().all()
    assert tm2["ideal"].isna().all()

    ne = out[out["molecule"] == "NOTINCSV"]
    assert ne["ref_label"].isna().all()

    assert len(skip_report) == 1
    assert skip_report[0]["molecule"] == "TESTMOL2"
    assert skip_report[0]["n_mismatched"] == 1
    mode_index, engine_freq, ds_freq = skip_report[0]["example"]
    assert mode_index == 1
    assert engine_freq == 999.0
    assert ds_freq == 1000.0


def test_attach_labels_missing_data_score_row_counts_as_mismatch():
    """A mode_index with no corresponding data_score.csv row at all (not
    just a numeric mismatch) must also disqualify the molecule's join and be
    reported with ds_freq=None -- not silently skipped past."""
    csv_tables = {"data_score": pd.DataFrame([
        {"molecule": "TESTMOL3", "mode": 1, "freq": 1000.0, "type": "stretch", "ideal": "yes"},
    ]), "characterised_modes": pd.DataFrame(columns=["molecule", "mode", "ref"])}
    df = pd.DataFrame([
        {"molecule": "TESTMOL3", "mode_index": 1, "kind": "internal",
         "freq": 1000.0, "ref_label": None, "ideal": None},
        {"molecule": "TESTMOL3", "mode_index": 2, "kind": "internal",
         "freq": 2000.0, "ref_label": None, "ideal": None},  # no data_score row for mode 2
    ])
    label_lookup = build_label_lookup(csv_tables)
    out, skip_report = attach_labels(df, csv_tables, label_lookup)
    assert out["ref_label"].isna().all()  # all-or-nothing: mode 1's match doesn't survive
    assert len(skip_report) == 1
    mode_index, engine_freq, ds_freq = skip_report[0]["example"]
    assert ds_freq is None


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
