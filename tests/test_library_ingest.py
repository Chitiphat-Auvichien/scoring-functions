"""Regression tests for src/library_ingest.py (formerly src/excel_ingest.py).

Run from ``Github/scoring-functions/``:
    py -m pytest tests/                (with pytest)
    py tests/test_library_ingest.py   (standalone; no pytest needed)

**Rewritten 2026-07-07 for the single, roster-driven pipeline** (final flip
to Gaussian-direct -- see IMPLEMENTATION_PLAN.md's RESUME HERE and
src/library_ingest.py's module docstring). The old dual-source
("source=excel"/"source=gaussian") architecture and its xlsx-workbook ingest
path are gone entirely; every one of ``data/mol_list_method.csv``'s roster
molecules now has a verified on-disk ``.log``/``.gjf`` pair, and
``library_scores.csv`` is built by iterating that roster directly. Tests
that used to check ``source="excel"``-specific behavior (the full-library
overlay, ``EXCEL_TO_LOG``/``resolve_excel_molecule_name``) are deleted
outright rather than adapted -- there is no equivalent concept left to test.

**2026-07-09: OH4/OF4 excluded from the roster** (72 -> 70 molecules).
Neither is a genuine stationary point at this project's MP2/3-21G level
(imaginary/negative frequencies -- e.g. OH4 mode 1 was -303.26 cm^-1), so
their "normal modes" are not physically meaningful vibrations of a real
minimum and cannot be validly compared to the TeH4 ideal see-saw template.
Every roster-count assertion in this file reflects 70, not 72.

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
    score_geometry_molecule, attach_labels, attach_ideal_tags,
    build_library_scores, multi_centre_molecules,
    out_of_calibration_scope_molecules, _central_atom_index,
    _basename_to_molecule_map, regenerate_characterised_modes,
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

def test_load_mol_roster_reads_all_77_molecules():
    """77 rows: 11 'ideal' + 56 'non-ideal' + 1 'multi-centre' (C6H6) + 9
    'test' (a held-out transferability-test set: CH4, C4H4, C10H16, PCl5,
    C3H6, B3N3H6, CHCl3, CH3CN, C3O3H6). 2026-08-11: the redundant
    `basename` column was dropped from mol_list_method.csv entirely --
    `molecule` now doubles as the on-disk basename (every on-disk file was
    already renamed to match `molecule` in an earlier commit). (Test name
    kept in the same style for history/grep-ability; the docstring is the
    source of truth for the current count.)"""
    roster = load_mol_roster(DATA_DIR)
    assert "molecule" in roster.columns
    assert "basename" not in roster.columns
    assert len(roster) == 77
    assert roster["molecule"].is_unique
    assert roster["molecule"].notna().all()
    counts = roster["mol_type"].value_counts()
    assert counts["ideal"] == 11
    assert counts["non-ideal"] == 56
    assert counts["multi-centre"] == 1
    assert counts["test"] == 9


def test_load_mol_roster_raises_on_missing_required_column():
    with tempfile.TemporaryDirectory() as tmp_dir:
        with open(os.path.join(tmp_dir, "mol_list_method.csv"), "w") as f:
            f.write("shape,foo\nbend,1\n")
        try:
            load_mol_roster(tmp_dir)
            assert False, "expected ValueError for a roster missing 'molecule'"
        except ValueError as e:
            assert "molecule" in str(e)


def test_resolve_log_basename_matches_roster_for_known_molecules():
    """resolve_log_basename() is now a pure roster-membership lookup (no
    on-disk existence check) that returns the molecule name itself (the
    `basename` column no longer exists) -- values must match
    data/mol_list_method.csv exactly."""
    assert resolve_log_basename("SbH3", DATA_DIR) == "SbH3"
    assert resolve_log_basename("BrF3", DATA_DIR) == "BrF3"
    assert resolve_log_basename("H2O", DATA_DIR) == "H2O"
    assert resolve_log_basename("C6H6", DATA_DIR) == "C6H6"
    assert resolve_log_basename("SeCl4", DATA_DIR) == "SeCl4"
    assert resolve_log_basename("CH4", DATA_DIR) == "CH4"  # test-category molecule


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
    fake_row = pd.DataFrame([{"molecule": "FAKE"}])
    combined = pd.concat([roster, fake_row], ignore_index=True)
    missing, _ = check_roster_disk_consistency(combined, DATA_DIR)
    assert missing == [("FAKE", "FAKE")]


def test_check_roster_disk_consistency_detects_orphaned_basename():
    """Dropping a real roster row leaves its on-disk basename (== molecule
    name) unreferenced -- exactly the gramicidin-style 'file present, no
    roster row' case."""
    roster = load_mol_roster(DATA_DIR)
    dropped_basename = roster.loc[roster["molecule"] == "TeH2", "molecule"].iloc[0]
    reduced = roster[roster["molecule"] != "TeH2"]
    _, orphaned = check_roster_disk_consistency(reduced, DATA_DIR)
    assert dropped_basename in orphaned


def test_discover_geometry_molecules_matches_roster_exactly():
    """The safety-net directory scan and the roster agree exactly (77
    basenames each) now that Phase A's coverage is complete."""
    bases = discover_geometry_molecules(DATA_DIR)
    assert bases == sorted(bases)
    assert len(bases) == len(set(bases))
    roster = load_mol_roster(DATA_DIR)
    assert set(bases) == set(roster["molecule"])


# ---------------------------------------------------------------------------
# build_library_scores() control flow: fail-fast on a missing roster file,
# warn-and-continue on an orphaned disk file -- exercised on tiny synthetic
# populations so these stay fast.
# ---------------------------------------------------------------------------

def test_build_library_scores_raises_filenotfounderror_on_missing_basename():
    fake_roster = pd.DataFrame([{"molecule": "not_a_real_basename_xyz"}])
    with patch("src.library_ingest.load_mol_roster", return_value=fake_roster):
        try:
            build_library_scores(DATA_DIR)
            assert False, "expected FileNotFoundError"
        except FileNotFoundError as e:
            assert "not_a_real_basename_xyz" in str(e)


def test_build_library_scores_warns_on_orphaned_disk_basename_and_still_builds():
    """Restrict the roster to just water (fast: 9 modes) while leaving the
    real 77-basename disk listing plus one extra fake basename -- missing
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
    simply absent from mol_list_method.csv's 77-row roster now, not an
    explicit ingest-time exclusion list."""
    df = _load()
    assert not (df["molecule"] == "Gly5").any()


def test_full_population_has_geometry_for_every_row():
    """Every row in the checked-in golden is geometry-backed (has_geometry
    unconditionally True) -- the whole point of the roster-driven flip.

    The checked-in data/results/library_scores.csv golden predates the
    2026-08-11 roster expansion (9 'test'-category transferability
    molecules added) and was NOT regenerated as part of that change
    (regenerating it is a separate, deliberate follow-up -- see
    IMPLEMENTATION_PLAN.md), so it is compared here against the roster's
    non-test subset (ideal/non-ideal/multi-centre, 68 rows) rather than the
    full 77-row roster. Population may be 67 or 68 molecules depending on
    whether every roster row's files happen to parse/score cleanly (a
    molecule with files present but unusable connectivity is
    warned-and-skipped, not fabricated)."""
    df = _load()
    roster = load_mol_roster(DATA_DIR)
    non_test_roster = roster[roster["mol_type"] != "test"]
    assert df["has_geometry"].all()
    assert len(df) > 0
    assert set(df["molecule"].unique()) <= set(non_test_roster["molecule"])
    assert df["molecule"].nunique() >= len(non_test_roster) - 1  # at most 1 legitimate skip


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
    """Per-bond |s_AB| (parsed out of the semicolon-joined string) sums back to
    V_Stretch for every internal row (eq:bondscore). s_AB is signed (stretch
    vs. compress); V_Stretch itself sums magnitudes. Every internal row
    carries bond detail regardless of label (score_geometry_molecule always
    calls score_bonds() directly, not classify_all_modes()'s filtered
    subset)."""
    df = _load()
    checked = 0
    for _, row in df[df["kind"] == "internal"].iterrows():
        assert isinstance(row["s_AB"], str) and row["s_AB"], \
            (row["molecule"], row["mode_index"])
        total = sum(abs(float(part.split(":")[1])) for part in row["s_AB"].split(";"))
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
        assert r["d_CA"] is None  # d_CA is an internal-mode-only quantity

    for r in internal:
        assert r["ref_label"] is None and r["ideal"] is None  # attached later
        assert r["has_geometry"] is True
        assert 0.0 <= r["V_Stretch"] <= 1.0
        assert isinstance(r["s_AB"], str) and r["s_AB"]
        assert r["delta_b_mean"] is not None and r["delta_b_mean"] >= 0.0
        # s_AB is signed (stretch vs. compress); V_Stretch sums magnitudes.
        total = sum(abs(float(part.split(":")[1])) for part in r["s_AB"].split(";"))
        assert abs(total - r["V_Stretch"]) < 1e-3
        # 2026-07-08 Gaussian-direct parser rework: real internal modes carry
        # engine-parsed mu/k/irrep.
        assert r["reduced_mass"] is not None and r["reduced_mass"] > 0.0
        assert r["force_constant"] is not None and r["force_constant"] > 0.0
        assert r["irrep"] is not None and isinstance(r["irrep"], str)
        # mol_type defaulted to None (not passed) -> d_CA never attempted.
        assert r["d_CA"] is None


# ---------------------------------------------------------------------------
# d_CA (central/hub-atom displacement amplitude) -- _central_atom_index()
# and score_geometry_molecule(mol_type=...)'s gating of it.
# ---------------------------------------------------------------------------

def test_central_atom_index_finds_the_unique_hub():

    # Water: O(0) bonded to H(1), H(2) -- O has degree 2 == n_atoms-1.
    assert _central_atom_index([(0, 1), (0, 2)], 3) == 0


def test_central_atom_index_returns_none_for_no_hub_or_multiple_candidates():

    # Benzene-like ring: every atom has degree 2, none has degree n-1=5.
    ring = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
    assert _central_atom_index(ring, 6) is None
    # Two atoms both happen to have degree n_atoms-1 (e.g. two bridging
    # centres) -- ambiguous, must return None, not guess.
    assert _central_atom_index([(0, 1), (0, 2), (1, 2)], 3) is None


def test_score_geometry_molecule_computes_d_ca_for_non_ideal_water():
    """Water tagged mol_type='non-ideal' (its real roster tag) gets a real
    d_CA: the oxygen (the unique central atom, bonded to both hydrogens) is
    atom index 0 by construction (GaussianParser preserves file order and
    O is listed first in H2O's geometry)."""

    base = resolve_log_basename("H2O", DATA_DIR)
    rows = score_geometry_molecule(base, DATA_DIR, mol_type="non-ideal")
    internal = [r for r in rows if r["kind"] == "internal"]
    external = [r for r in rows if r["kind"] == "external"]
    for r in internal:
        assert r["d_CA"] is not None and r["d_CA"] >= 0.0
    for r in external:
        assert r["d_CA"] is None


def test_score_geometry_molecule_skips_d_ca_for_multi_centre_or_unset_mol_type():
    base = resolve_log_basename("H2O", DATA_DIR)
    for mol_type in (None, "multi-centre", "something-unrecognized"):
        rows = score_geometry_molecule(base, DATA_DIR, mol_type=mol_type)
        assert all(r["d_CA"] is None for r in rows if r["kind"] == "internal")


def test_multi_centre_molecules_tags_c6h6_and_nothing_else_in_the_real_roster():
    """multi_centre_molecules() must reproduce the real mol_list_method.csv's
    mol_type=='multi-centre' rows exactly -- as of the 77-molecule roster
    (2026-08-11: 'test' category added, basename column dropped), that is
    C6H6 alone (verified, not assumed). No longer drives
    src.calibrate.SINGLE_CENTRE_ONLY_EXCLUDE directly (that's
    out_of_calibration_scope_molecules() now); still used for figure
    filtering."""

    result = multi_centre_molecules(DATA_DIR)
    assert result == frozenset({"C6H6"})


def test_out_of_calibration_scope_molecules_is_c6h6_plus_the_9_test_molecules():
    """out_of_calibration_scope_molecules() (drives
    src.calibrate.SINGLE_CENTRE_ONLY_EXCLUDE) is the complement of the
    ideal/non-ideal inclusion filter: mol_type not in {'ideal', 'non-ideal'}
    -- currently C6H6 (multi-centre) plus the 9 held-out 'test'
    transferability molecules, 10 names total."""
    result = out_of_calibration_scope_molecules(DATA_DIR)
    expected = {
        "C6H6", "CH4", "C4H4", "C10H16", "PCl5", "C3H6", "B3N3H6",
        "CHCl3", "CH3CN", "C3O3H6",
    }
    assert result == frozenset(expected)
    assert len(result) == 10


def test_schema_columns_includes_mu_k_irrep_and_d_ca():
    """reduced_mass/force_constant/irrep (2026-07-08 parser rework) and d_CA
    (2026-07-08 data_score.csv retirement) are appended to the end of the
    locked schema -- purely additive, existing columns untouched."""
    assert SCHEMA_COLUMNS[-4:] == ["reduced_mass", "force_constant", "irrep", "d_CA"]
    assert SCHEMA_COLUMNS[:-4] == [
        "molecule", "mode_index", "kind", "freq", "ref_label", "ideal",
        "V_Stretch", "delta_b_mean", "s_AB", "rel_db", "has_geometry",
        "predicted_label", "predicted_annotation",
        "Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "ref_key",
    ]


# ---------------------------------------------------------------------------
# attach_labels() -- fully synthetic (no CSV I/O) unit tests of the
# frequency-gating logic. 2026-07-08: repointed off data_score.csv onto
# characterised_modes.csv; 'ideal' is no longer set by attach_labels() at
# all (see attach_ideal_tags() tests below) -- get_label() now returns a
# (ref_label, ref_key) 2-tuple.
# ---------------------------------------------------------------------------

def test_attach_labels_gates_on_frequency_and_leaves_unlisted_molecules_alone():
    """A matching, frequency-consistent internal row gets ref_label filled
    in; a mismatched (or characterised_modes.csv-absent) mode disqualifies
    the WHOLE molecule's internal-row label join (all-or-nothing gate); a
    molecule with no characterised_modes.csv counterpart at all is untouched
    (not an error); external rows are never touched. 'ideal' is untouched by
    this function entirely (sourced separately by attach_ideal_tags())."""
    csv_tables = {"characterised_modes": pd.DataFrame([
        {"molecule": "TESTMOL", "mode": 1, "freq": 1000.0, "type": "stretch", "ref": "Foo1970"},
        {"molecule": "TESTMOL", "mode": 2, "freq": 500.0, "type": "bend", "ref": None},
        {"molecule": "TESTMOL2", "mode": 1, "freq": 1000.0, "type": "stretch", "ref": None},
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
        # NOTINCSV: no characterised_modes.csv counterpart at all -> untouched, no warning.
        {"molecule": "NOTINCSV", "mode_index": 1, "kind": "internal",
         "freq": 42.0, "ref_label": None, "ideal": None},
    ])

    label_lookup = build_label_lookup(csv_tables)
    out, skip_report = attach_labels(df, csv_tables, label_lookup)

    tm = out[out["molecule"] == "TESTMOL"]
    assert tm.loc[tm["mode_index"] == 1, "ref_label"].iloc[0] == "stretch"
    assert tm.loc[tm["mode_index"] == 1, "ref_key"].iloc[0] == "Foo1970"
    assert tm.loc[tm["mode_index"] == 2, "ref_label"].iloc[0] == "bend"
    # attach_labels() never touches 'ideal' -- stays whatever it was passed in as.
    assert tm["ideal"].isna().all()
    assert tm.loc[tm["kind"] == "external", "ref_label"].iloc[0] == "translation"

    tm2 = out[out["molecule"] == "TESTMOL2"]
    assert tm2["ref_label"].isna().all()

    ne = out[out["molecule"] == "NOTINCSV"]
    assert ne["ref_label"].isna().all()

    assert len(skip_report) == 1
    assert skip_report[0]["molecule"] == "TESTMOL2"
    assert skip_report[0]["n_mismatched"] == 1
    mode_index, engine_freq, cm_freq = skip_report[0]["example"]
    assert mode_index == 1
    assert engine_freq == 999.0
    assert cm_freq == 1000.0


def test_attach_labels_missing_characterised_modes_row_counts_as_mismatch():
    """A mode_index with no corresponding characterised_modes.csv row at all
    (not just a numeric mismatch) must also disqualify the molecule's join
    and be reported with cm_freq=None -- not silently skipped past."""
    csv_tables = {"characterised_modes": pd.DataFrame([
        {"molecule": "TESTMOL3", "mode": 1, "freq": 1000.0, "type": "stretch", "ref": None},
    ])}
    df = pd.DataFrame([
        {"molecule": "TESTMOL3", "mode_index": 1, "kind": "internal",
         "freq": 1000.0, "ref_label": None, "ideal": None},
        {"molecule": "TESTMOL3", "mode_index": 2, "kind": "internal",
         "freq": 2000.0, "ref_label": None, "ideal": None},  # no characterised_modes row for mode 2
    ])
    label_lookup = build_label_lookup(csv_tables)
    out, skip_report = attach_labels(df, csv_tables, label_lookup)
    assert out["ref_label"].isna().all()  # all-or-nothing: mode 1's match doesn't survive
    assert len(skip_report) == 1
    mode_index, engine_freq, cm_freq = skip_report[0]["example"]
    assert cm_freq is None


# ---------------------------------------------------------------------------
# attach_ideal_tags() -- 'ideal' is now sourced structurally from
# mol_list_method.csv's per-molecule mol_type column, unconditional (not
# gated by attach_labels()'s frequency check).
# ---------------------------------------------------------------------------

def test_attach_ideal_tags_maps_mol_type_onto_internal_rows_only():

    roster = pd.DataFrame([
        {"molecule": "IDEALMOL", "mol_type": "ideal"},
        {"molecule": "NONIDEALMOL", "mol_type": "non-ideal"},
        {"molecule": "MULTIMOL", "mol_type": "multi-centre"},
    ])
    df = pd.DataFrame([
        {"molecule": "IDEALMOL", "kind": "internal", "ideal": None},
        {"molecule": "IDEALMOL", "kind": "external", "ideal": None},
        {"molecule": "NONIDEALMOL", "kind": "internal", "ideal": None},
        {"molecule": "MULTIMOL", "kind": "internal", "ideal": None},
    ])
    out = attach_ideal_tags(df, roster)
    assert out.loc[(out["molecule"] == "IDEALMOL") & (out["kind"] == "internal"), "ideal"].iloc[0] == "yes"
    assert out.loc[(out["molecule"] == "IDEALMOL") & (out["kind"] == "external"), "ideal"].isna().all()
    assert out.loc[out["molecule"] == "NONIDEALMOL", "ideal"].iloc[0] == "no"
    assert out.loc[out["molecule"] == "MULTIMOL", "ideal"].isna().all()


def test_attach_ideal_tags_is_unconditional_even_without_a_ref_label():
    """Unlike the old data_score.csv-sourced 'ideal', attach_ideal_tags()
    does not require a successful ref_label join -- it is a structural,
    per-molecule property."""

    roster = pd.DataFrame([{"molecule": "NOLABELMOL", "mol_type": "ideal"}])
    df = pd.DataFrame([
        {"molecule": "NOLABELMOL", "kind": "internal", "ref_label": None, "ideal": None},
    ])
    out = attach_ideal_tags(df, roster)
    assert out["ideal"].iloc[0] == "yes"
    assert out["ref_label"].isna().all()  # unaffected -- separate concern


# ---------------------------------------------------------------------------
# resync_reference_metadata() (2026-07-08, Phase 2 of the "Gaussian-direct
# intermediate file" plan) -- fixes attach_labels()'s stale-freq gate at the
# source by resyncing characterised_modes.csv's freq/k(/mu) from the SAME
# on-disk log the engine scores. (Also resynced data_score.csv until that
# file was deleted from disk 2026-07-09 -- see the function's docstring.)
# ---------------------------------------------------------------------------

def test_resync_reference_metadata_real_roster_has_no_mismatches_or_missing_logs():
    """Live, read-only (write=False) check against the real, already-fixed
    data/ tree.

    **2026-08-11** (mol_list_method.csv format migration): the roster grew
    from 68 to 77 molecules (the 9 held-out 'test'-category transferability
    molecules -- CH4, C4H4, C10H16, PCl5, C3H6, B3N3H6, CHCl3, CH3CN,
    C3O3H6 -- were already present on disk with matching
    characterised_modes.csv rows, so this stays a clean 0-mismatch resync).
    resync_reference_metadata() iterates the roster directly (not the
    checked-in library_scores.csv golden, which was NOT regenerated as part
    of that migration), so it sees all 77 rows here."""
    from src.library_ingest import resync_reference_metadata

    report = resync_reference_metadata(DATA_DIR, write=False)
    assert report["skipped_mode_count_mismatch"] == []
    assert report["skipped_no_log"] == []
    assert report["skipped_no_rows"] == []
    assert len(report["resynced"]) == 77


def test_resync_reference_metadata_is_a_true_dry_run_when_write_false():
    """write=False must not touch characterised_modes.csv on disk -- a live
    safety check, not just a docstring promise (this test would fail
    loudly, corrupting the real tracked CSV, if resync_reference_metadata
    ever stopped honoring write=False). data_score.csv is no longer part of
    this function's scope (deleted from disk 2026-07-09)."""
    from src.library_ingest import resync_reference_metadata

    cm_path = os.path.join(DATA_DIR, "characterised_modes.csv")
    with open(cm_path, "rb") as f:
        cm_before = f.read()

    resync_reference_metadata(DATA_DIR, write=False)

    with open(cm_path, "rb") as f:
        assert f.read() == cm_before


def _write_csv_utf8sig(path, rows, columns):
    pd.DataFrame(rows, columns=columns).to_csv(path, index=False, encoding="utf-8-sig")


def test_resync_reference_metadata_synthetic_fixture():
    """Fully isolated (tempdir, no real CSVs touched): a stale-freq/k row
    gets corrected from the on-disk log, `irrep` is left untouched even
    though it disagrees with the engine's raw ASCII token (A1 vs 'A₁' --
    see src/library_ingest.py's module comment on why irrep resync is
    deliberately out of scope), a mode-count-mismatch molecule is skipped
    and reported (not guessed), and a molecule with zero rows in
    characterised_modes.csv is reported under skipped_no_rows, not silently
    ignored or fabricated.

    **2026-07-09 (OH4/OF4 exclusion session):** `resync_reference_metadata`
    no longer reads/writes `data_score.csv` (deleted from disk this session,
    an already-retired file) -- this fixture was rewritten to put all the
    stale/mismatch fixture rows in `characterised_modes.csv` directly
    instead of `data_score.csv`.

    **2026-08-11:** the `basename` column was dropped (molecule name ==
    on-disk basename now), so this fixture's water log/gjf pair is copied
    under each synthetic molecule's own name (STALEMOL.log/.com,
    MISMATCHMOL.log/.com) instead of being shared via a basename column."""
    import shutil
    from src.library_ingest import resync_reference_metadata

    with tempfile.TemporaryDirectory() as tmp_dir:
        os.makedirs(os.path.join(tmp_dir, "logs"))
        os.makedirs(os.path.join(tmp_dir, "gjf"))
        for mol in ("STALEMOL", "MISMATCHMOL"):
            shutil.copy(os.path.join(DATA_DIR, "logs", "H2O.log"),
                        os.path.join(tmp_dir, "logs", f"{mol}.log"))
            shutil.copy(os.path.join(DATA_DIR, "gjf", "H2O.com"),
                        os.path.join(tmp_dir, "gjf", f"{mol}.com"))

        # Roster: STALEMOL (real water log/gjf, stale CSV values to fix),
        # MISMATCHMOL (own copy of the same log, but the CSV claims a 4th
        # mode that doesn't exist -- 3 engine modes vs. 4 claimed),
        # NOROWSMOL (no CSV row at all, and no log copy needed -- not an
        # error, just nothing to resync).
        roster = pd.DataFrame([
            {"molecule": "STALEMOL"},
            {"molecule": "MISMATCHMOL"},
            {"molecule": "NOROWSMOL"},
        ])
        roster.to_csv(os.path.join(tmp_dir, "mol_list_method.csv"), index=False)

        cm_cols = ["molecule", "mode", "freq", "μ", "k", "irrep"]
        cm_rows = [
            # STALEMOL: mode 1 deliberately stale (engine says 1722.457), mode
            # 2/3 already exact -- proves per-mode-only-when-changed reporting.
            {"molecule": "STALEMOL", "mode": 1, "freq": 1700.0, "μ": 1.0,
             "k": 1.0, "irrep": "A₁"},
            {"molecule": "STALEMOL", "mode": 2, "freq": 3501.5073, "μ": 1.0,
             "k": 7.4906, "irrep": "A₁"},
            {"molecule": "STALEMOL", "mode": 3, "freq": 3660.7973, "μ": 1.0,
             "k": 8.5478, "irrep": "B₂"},
            # MISMATCHMOL: 4 rows claimed, engine (same water log) has only 3.
            {"molecule": "MISMATCHMOL", "mode": 1, "freq": 1700.0, "μ": 1.0,
             "k": 1.0, "irrep": "A₁"},
            {"molecule": "MISMATCHMOL", "mode": 2, "freq": 3501.0, "μ": 1.0,
             "k": 7.5, "irrep": "A₁"},
            {"molecule": "MISMATCHMOL", "mode": 3, "freq": 3660.0, "μ": 1.0,
             "k": 8.5, "irrep": "B₂"},
            {"molecule": "MISMATCHMOL", "mode": 4, "freq": 9999.0, "μ": 1.0,
             "k": 9.9, "irrep": "X"},
        ]
        _write_csv_utf8sig(os.path.join(tmp_dir, "characterised_modes.csv"), cm_rows, cm_cols)

        report = resync_reference_metadata(tmp_dir, write=True)

        assert report["skipped_no_rows"] == ["NOROWSMOL"]
        assert len(report["skipped_mode_count_mismatch"]) == 1
        assert report["skipped_mode_count_mismatch"][0]["molecule"] == "MISMATCHMOL"

        resynced_by_mol = {e["molecule"]: e for e in report["resynced"]}
        assert set(resynced_by_mol) == {"STALEMOL"}
        stale_changes = resynced_by_mol["STALEMOL"]["changes"]
        # Only mode 1's freq actually changed (modes 2/3 were already exact).
        freq_changes = [c for c in stale_changes if c["field"] == "freq"]
        assert len(freq_changes) == 1
        assert freq_changes[0]["mode"] == 1
        assert freq_changes[0]["old"] == "1700.0"
        assert freq_changes[0]["new"] == "1722.457"

        cm_after = pd.read_csv(os.path.join(tmp_dir, "characterised_modes.csv"), dtype=str)
        stale_after = cm_after[cm_after["molecule"] == "STALEMOL"].set_index("mode")
        assert stale_after.loc["1", "freq"] == "1722.457"
        assert stale_after.loc["2", "freq"] == "3501.5073"
        assert stale_after.loc["3", "freq"] == "3660.7973"
        # irrep is deliberately NEVER touched, even for the corrected mode.
        assert stale_after.loc["1", "irrep"] == "A₁"

        # MISMATCHMOL's row was left completely untouched (still the stale,
        # never-corrected values) -- a real structural problem, not guessed.
        cm_mismatch_after = cm_after[cm_after["molecule"] == "MISMATCHMOL"].set_index("mode")
        assert cm_mismatch_after.loc["2", "freq"] == "3501.0"


# ---------------------------------------------------------------------------
# regenerate_characterised_modes() (2026-07-08) -- direct disk scan
# (NOT restricted to mol_list_method.csv's roster), preserving manual
# literature columns and reporting (never silently dropping) any molecule
# that loses its on-disk log/gjf pair.
# ---------------------------------------------------------------------------

def test_basename_to_molecule_map_translates_known_roster_basenames():
    m = _basename_to_molecule_map(DATA_DIR)
    assert m["SbH3"] == "SbH3"
    assert m["C6H6"] == "C6H6"
    assert m["H2O"] == "H2O"
    assert m["CH4"] == "CH4"  # test-category molecule


def test_regenerate_characterised_modes_dry_run_against_real_tree():
    """Live, read-only (write=False) check against the real data/ tree,
    which is IDEMPOTENT under regenerate_characterised_modes(): the file
    already covers all on-disk basenames, so a fresh dry run reports old ==
    new molecule counts, nothing dropped, nothing added.

    **2026-08-11** (mol_list_method.csv format migration): the roster's
    'test' category (9 held-out transferability molecules) was added, and
    their .log/.gjf pairs plus characterised_modes.csv rows were already on
    disk before this session -- this disk-scan-based function (unrestricted
    by the roster) already covered them, so the count moves from 68 to 77
    but idempotency is unaffected."""
    report = regenerate_characterised_modes(DATA_DIR, write=False)
    assert report["n_disk_basenames"] == 77
    assert report["n_old_molecules"] == 77
    assert report["n_new_molecules"] == 77
    assert report["dropped_molecules"] == []
    assert report["added_molecules"] == []
    assert report["parse_failures"] == []
    assert report["n_rows_written"] > 0
    # Every row already existed at its (molecule, mode) key (idempotent
    # regeneration), so every row's manual columns are "preserved".
    assert report["n_rows_with_preserved_manual_labels"] == report["n_rows_written"]


def test_regenerate_characterised_modes_is_a_true_dry_run_when_write_false():
    cm_path = os.path.join(DATA_DIR, "characterised_modes.csv")
    with open(cm_path, "rb") as f:
        before = f.read()
    regenerate_characterised_modes(DATA_DIR, write=False)
    with open(cm_path, "rb") as f:
        assert f.read() == before


def test_regenerate_characterised_modes_synthetic_fixture_preserves_manual_columns():
    """Fully isolated (tempdir, no real CSVs touched): an existing
    (molecule, mode) row's manual columns (shape/type/sym/description/
    νₖ/ref/Note) survive a regeneration untouched, its irrep is PRESERVED
    (not overwritten with the raw engine token, even though it differs --
    see regenerate_characterised_modes()'s docstring for why), a genuinely
    new row gets blank manual columns + the raw engine irrep, and a
    molecule with no on-disk match anymore is reported dropped."""
    import shutil

    with tempfile.TemporaryDirectory() as tmp_dir:
        os.makedirs(os.path.join(tmp_dir, "logs"))
        os.makedirs(os.path.join(tmp_dir, "gjf"))
        shutil.copy(os.path.join(DATA_DIR, "logs", "H2O.log"),
                    os.path.join(tmp_dir, "logs", "H2O.log"))
        shutil.copy(os.path.join(DATA_DIR, "gjf", "H2O.com"),
                    os.path.join(tmp_dir, "gjf", "H2O.com"))
        # No roster file at all -- regenerate_characterised_modes() must not
        # need one (disk-scan only); _basename_to_molecule_map() falls back
        # to the raw basename when the roster can't be read.

        old_cols = ["shape", "molecule", "mode", "freq", "μ", "k", "irrep",
                    "type", "sym", "description", "νₖ", "ref", "Note"]
        old_rows = [
            # Existing row for the water log's mode 1 -- manual columns +
            # a deliberately "wrong" (but author-curated) irrep that must
            # survive untouched.
            {"shape": "bend", "molecule": "H2O", "mode": "1",
             "freq": "0", "μ": "0", "k": "0", "irrep": "HAND_RESOLVED",
             "type": "bend", "sym": "-", "description": "scissor", "νₖ": "ν2",
             "ref": "Foo1970", "Note": ""},
            # A molecule with no on-disk match anymore -- must be dropped
            # and reported.
            {"shape": "-", "molecule": "GONEMOL", "mode": "1", "freq": "0",
             "μ": "0", "k": "0", "irrep": "A", "type": "stretch", "sym": "-",
             "description": "", "νₖ": "", "ref": "", "Note": ""},
        ]
        _write_csv_utf8sig(os.path.join(tmp_dir, "characterised_modes.csv"), old_rows, old_cols)

        report = regenerate_characterised_modes(tmp_dir, write=True)

        assert report["dropped_molecules"] == ["GONEMOL"]
        assert report["added_molecules"] == []  # H2O already had a row
        assert report["n_rows_written"] == 3  # water: 3 internal modes
        assert report["n_rows_with_preserved_manual_labels"] == 1  # only mode 1

        new_df = pd.read_csv(os.path.join(tmp_dir, "characterised_modes.csv"),
                              dtype=str, keep_default_na=False)
        assert set(new_df["molecule"]) == {"H2O"}
        row1 = new_df[new_df["mode"] == "1"].iloc[0]
        assert row1["type"] == "bend" and row1["ref"] == "Foo1970"
        assert row1["irrep"] == "HAND_RESOLVED"  # preserved, not overwritten
        assert float(row1["freq"]) != 0.0  # freq WAS overwritten from the real log

        row2 = new_df[new_df["mode"] == "2"].iloc[0]
        assert row2["type"] == "" and row2["ref"] == ""  # genuinely new -> blank manual cols
        assert row2["irrep"] != ""  # but DOES get the raw engine irrep (nothing to lose)


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
