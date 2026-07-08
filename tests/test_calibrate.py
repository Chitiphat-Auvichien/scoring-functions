"""Regression tests for src/calibrate.py (Phase 3 tau calibration) and the
label-level validation it unlocks.

Run from ``Github/scoring-functions/``:
    py -m pytest tests/                (with pytest)
    py tests/test_calibrate.py        (standalone; no pytest needed)

Design note (mirrors tests/test_library_ingest.py): the plateau sweep itself is
cheap (small molecules, ~190 grid points in well under a minute) but does
depend on data/logs/ parsing for the geometry-backed library + benzene EMIT; the two derived
artifacts (``data/results/thresholds.json``,
``data/results/tau_sensitivity_sweep.csv``) are already committed goldens, so
most tests here check those directly rather than re-running the sweep. One
test (`test_sweep_tau_tr_reproduces_frozen_sweep_on_a_reduced_grid`) does
re-run `sweep_tau_tr` on a smaller grid as an integration check.

Threshold policy (see tests/test_classifier.py's docstring for the other
half of this decision): this file is the dedicated home for validating
Thresholds.calibrated() -- the Phase-3 CALIBRATED values loaded from
thresholds.json -- re-confirming the exact same behavioral targets
test_classifier.py pins under the PROVISIONAL Thresholds() defaults. Task 4's
explicit ask ("re-verify benzene EMIT 34/35/36 under calibrated thresholds")
is satisfied here.
"""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import pandas as pd                                                 # noqa: E402

from main import load_inputs, build_scorer_and_final                # noqa: E402
from src.classifier import (                                        # noqa: E402
    classify_all_modes, Thresholds,
    is_external_label, is_mixed_external,
    STRETCHING, BENDING, MIXED_STRETCH_BEND,
)
from src.calibrate import (                                         # noqa: E402
    find_plateau, freeze_tau_tr, derive_stretch_bend_thresholds,
    confusion_matrix_stats, filter_single_centre_library,
    SINGLE_CENTRE_ONLY_EXCLUDE,
)

THRESHOLDS_JSON = os.path.join(ROOT, "data", "results", "thresholds.json")
SWEEP_CSV = os.path.join(ROOT, "data", "results", "tau_sensitivity_sweep.csv")
LIB_CSV = os.path.join(ROOT, "data", "results", "library_scores.csv")


def _classify(mol, mode_type, thresholds):
    raw, _ = load_inputs(mol, mode_type, os.path.join(ROOT, "data"))
    scorer, final = build_scorer_and_final(raw, mode_type)
    scored = classify_all_modes(scorer, final, thresholds)
    return {m["name"]: m for m in scored}


def test_thresholds_json_is_well_formed():
    with open(THRESHOLDS_JSON) as f:
        data = json.load(f)
    assert 0.0 < data["tau_B"] < data["tau_S"] < 1.0
    assert 0.0 < data["tau_TR"] < 1.0
    lo, hi = data["plateau_tau_TR_range"]
    assert lo <= data["tau_TR"] <= hi


def test_calibrated_thresholds_differ_from_but_are_close_to_provisional():
    """The whole point of calibration: derive real numbers from data, not
    just reproduce the provisional placeholders -- but they should land in
    the same ballpark (the manuscript's own "~0.9 / ~0.2" language)."""
    calibrated = Thresholds.calibrated()
    provisional = Thresholds()
    assert calibrated.tau_TR == provisional.tau_TR == 0.95  # unaffected by this library
    assert calibrated.tau_S != provisional.tau_S
    assert calibrated.tau_B != provisional.tau_B
    assert abs(calibrated.tau_S - 0.9) < 0.02
    assert abs(calibrated.tau_B - 0.2) < 0.03


# Historical reference only (NOT test-enforced as SINGLE_CENTRE_ONLY_EXCLUDE's
# value anymore -- see below): the 2026-07-05 scope decision originally named
# these 7 molecules as two-, four-, or six-centre topologies that don't fit
# the "hydride library" single-centre AB_n scope.
_HISTORICAL_SINGLE_CENTRE_EXCLUDE_7 = frozenset({
    "C2H2", "C2H4", "C2H6", "H2O2", "C6H6", "iso-C4H10", "n-C4H10",
})


def test_single_centre_only_exclude_matches_scope_decision():
    """**2026-07-08: repointed to be roster-derived**
    (src.library_ingest.multi_centre_molecules()) rather than the original
    hardcoded 7-name frozenset -- see src/calibrate.py's own comment above
    SINGLE_CENTRE_ONLY_EXCLUDE for the full rationale. Six of the historical
    7 names are outside the finalized 72-molecule roster entirely (removed
    2026-07-07) and were always no-ops when applied to the real
    library_scores.csv; only C6H6 was ever actually present to be dropped,
    and mol_list_method.csv already tags it mol_type=='multi-centre'
    (verified, not assumed). H2O is deliberately NOT excluded (it is
    single-centre itself and stays in the non-ideal population as a
    familiar illustrative molecule)."""
    assert SINGLE_CENTRE_ONLY_EXCLUDE == {"C6H6"}
    assert "H2O" not in SINGLE_CENTRE_ONLY_EXCLUDE

    # Byte-identical FILTERING behavior check: the roster-derived set and
    # the historical hardcoded 7-name set must drop exactly the same
    # molecules from the real, current library_scores.csv.
    lib_df = pd.read_csv(LIB_CSV)
    old_dropped = set(lib_df.loc[lib_df["molecule"].isin(_HISTORICAL_SINGLE_CENTRE_EXCLUDE_7),
                                  "molecule"].unique())
    new_dropped = set(lib_df.loc[lib_df["molecule"].isin(SINGLE_CENTRE_ONLY_EXCLUDE),
                                  "molecule"].unique())
    assert old_dropped == new_dropped == {"C6H6"}


def test_filter_single_centre_library_drops_exactly_the_excluded_molecules():
    """filter_single_centre_library() is DATA-LOSS-FREE at the CSV level
    (library_scores.csv keeps every molecule's rows -- see
    src/library_ingest.py); this checks the analysis-time filter itself
    removes exactly whichever of the excluded molecules are actually present
    and nothing else.

    **2026-07-07 update (roster-driven pipeline):** six of
    SINGLE_CENTRE_ONLY_EXCLUDE's seven molecules (C2H2, C2H4, C2H6, H2O2,
    iso-C4H10, n-C4H10) are two-, four-, or six-centre topologies that are
    also absent from data/mol_list_method.csv's finalized 72-molecule
    roster entirely -- they no longer appear in library_scores.csv at all
    (not merely filtered out), so filter_single_centre_library() is a no-op
    for them. Only C6H6 (benzene, six-centre, still in the roster) is
    actually present to be dropped by this filter now."""
    lib_df = pd.read_csv(LIB_CSV)
    before_molecules = set(lib_df["molecule"].unique())
    present_of_excluded = SINGLE_CENTRE_ONLY_EXCLUDE & before_molecules
    assert present_of_excluded == {"C6H6"}

    filtered = filter_single_centre_library(lib_df)
    after_molecules = set(filtered["molecule"].unique())

    assert before_molecules - after_molecules == {"C6H6"}
    assert lib_df["molecule"].nunique() - 1 == filtered["molecule"].nunique()
    # filter_single_centre_library is idempotent (re-filtering drops nothing more).
    assert set(filter_single_centre_library(filtered)["molecule"].unique()) == after_molecules


def test_confusion_matrix_stats_applies_single_centre_filter_internally():
    """confusion_matrix_stats() must apply the scope filter itself,
    regardless of whether the caller already pre-filtered -- this is what
    makes plot_confusion_matrix's rigorous_df/nonideal_df tiers AND the
    direct-call regression tests above both automatically correct without
    each needing their own copy of the filter."""
    lib_df = pd.read_csv(LIB_CSV)
    calibrated = Thresholds.calibrated()

    # Calling with the raw (unfiltered) library and with an already-filtered
    # one must give identical results -- the filter is idempotent and
    # applied unconditionally inside the function. Compared field-by-field
    # (not dict equality) since some fields are legitimately NaN (e.g.
    # translation/rotation's recall_nonideal -- no external row is ever
    # ideal=='no') and NaN != NaN under plain ``==``.
    res_raw = confusion_matrix_stats(lib_df, calibrated)
    res_prefiltered = confusion_matrix_stats(filter_single_centre_library(lib_df), calibrated)
    for cat in ("stretch", "bend", "translation", "rotation"):
        entry_raw = res_raw["per_category"][cat]
        entry_pre = res_prefiltered["per_category"][cat]
        assert entry_raw.keys() == entry_pre.keys()
        for key in entry_raw:
            a, b = entry_raw[key], entry_pre[key]
            if isinstance(a, float) and a != a:  # NaN
                assert isinstance(b, float) and b != b
            else:
                assert a == b, (cat, key, a, b)


def test_derive_stretch_bend_thresholds_matches_frozen_values():
    lib_df = pd.read_csv(LIB_CSV)
    tau_S, tau_B, stats = derive_stretch_bend_thresholds(lib_df)
    with open(THRESHOLDS_JSON) as f:
        frozen = json.load(f)
    assert abs(tau_S - frozen["tau_S"]) < 1e-9
    assert abs(tau_B - frozen["tau_B"]) < 1e-9
    assert stats["gap_width"] > 0.5  # a wide, clean, zero-overlap gap


def test_derive_stretch_bend_thresholds_rejects_overlapping_populations():
    """Defensive-coding check: if the ideal stretch/bend populations ever
    overlapped, the derivation must fail loudly, not silently round past it."""
    bad_df = pd.DataFrame([
        {"kind": "internal", "ideal": "yes", "ref_label": "stretch", "V_Stretch": 0.5},
        {"kind": "internal", "ideal": "yes", "ref_label": "bend", "V_Stretch": 0.6},
    ])
    try:
        derive_stretch_bend_thresholds(bad_df)
        assert False, "expected ValueError for overlapping populations"
    except ValueError:
        pass


def test_tau_sensitivity_sweep_plateau_and_accuracy():
    sweep = pd.read_csv(SWEEP_CSV)
    # Accuracy (library normal-mode T/R ground truth) is 100% everywhere --
    # exact completeness for normal modes holds across the whole tau_TR grid.
    assert (sweep["accuracy"] == 1.0).all()
    # Exactly one grid step shows a label change (benzene EMIT 19's Rz
    # assignment, |score|=0.333, crossing the grid at tau_TR~0.335 -- the
    # sole tau_TR-sensitive transition anywhere in the evaluated data; every
    # other external mode's classification is gated by tau_B, not tau_TR,
    # across this whole range).
    nonzero = sweep[sweep["label_change_fraction"] > 0]
    assert len(nonzero) == 1, nonzero
    lo, hi = find_plateau(sweep)
    tau_TR, plateau = freeze_tau_tr(sweep)
    assert (lo, hi) == plateau
    assert tau_TR == 0.95
    assert lo <= 0.95 <= hi
    assert hi - lo > 0.5  # a wide plateau, not a razor's edge


def test_benzene_emit_34_35_36_under_calibrated_thresholds():
    """Task-4 re-check: EMIT 34/35 flagged, EMIT 36 the documented blind
    spot -- re-verified under Thresholds.calibrated() (tau_S 0.9->0.90368,
    tau_B 0.2->0.17327), not just the provisional 0.95/0.9/0.2 defaults.
    V_Stretch for 34/35 (0.667/0.577) is far above the calibrated tau_B
    (0.17327) either way, so this is robust to the shift, as predicted.
    (2026-07-07: "benzene" -> the finalized roster basename
    C6H6_MP2_3-21G_D6h -- content-preserving rename, numbers unchanged.)"""
    calibrated = Thresholds.calibrated()
    t = _classify("C6H6_MP2_3-21G_D6h", "emit", calibrated)
    assert t["EMIT 34"]["classification"] == "Tx*"
    assert t["EMIT 34"]["annotation"] == f"vibration={MIXED_STRETCH_BEND}"
    assert t["EMIT 35"]["classification"] == "Ty*"
    assert t["EMIT 35"]["annotation"] == f"vibration={MIXED_STRETCH_BEND}"
    assert t["EMIT 36"]["classification"] == "Tz"


def test_water_targets_under_calibrated_thresholds():
    """Task-4 re-check of test_classifier.py's water targets, under
    Thresholds.calibrated() instead of the pinned provisional defaults.
    (2026-07-07: "water" -> the finalized roster basename H2O-MP2-321G, a
    genuinely different/corrected calculation from the old water.log -- see
    tests/test_scores.py's module docstring -- but its V_Stretch values land
    on the same side of tau_S/tau_B, so the classification buckets checked
    here are unaffected.)"""
    calibrated = Thresholds.calibrated()
    t = _classify("H2O-MP2-321G", "normal", calibrated)
    for lbl in ("Tx", "Ty", "Tz", "Rx", "Ry", "Rz"):
        assert t[lbl]["classification"] == lbl
    assert t["Vib 1"]["classification"] == BENDING
    assert t["Vib 2"]["classification"] == STRETCHING
    assert t["Vib 3"]["classification"] == STRETCHING


def test_degenerate_emit_block_sanity_check():
    """Sanity check only (Decision 8 explicitly retracts any block-consuming
    MECHANISM; this test confirms no such mechanism is needed, it does not
    add one). Benzene EMIT 1-9 share one eigenvalue (a genuinely degenerate
    9-fold block). Plain one-to-one Hungarian assignment claims AT MOST 2 of
    them for external slots (Rx, Ry) purely on individual |score|, and the
    other 7 fall through to Step 4 with their OWN individually-computed
    s[V_S] -- which correctly differ from each other (they are 9 distinct
    eigenvectors, not 9 copies of one mode), so DIFFERENT vib_label outcomes
    within the block is the CORRECT behavior, not an inconsistency. What must
    hold: (a) the assignment is fully deterministic (same input -> byte-
    identical output, run twice), (b) exactly n_T+n_R=6 modes total across
    the whole 36-mode pool receive an external-slot classification, with no
    over- or under-assignment caused by the degeneracy.
    """
    calibrated = Thresholds.calibrated()
    raw, _ = load_inputs("C6H6_MP2_3-21G_D6h", "emit", os.path.join(ROOT, "data"))
    scorer, final = build_scorer_and_final(raw, "emit")

    scored_a = classify_all_modes(scorer, final, calibrated)
    scored_b = classify_all_modes(scorer, final, calibrated)
    labels_a = [m["classification"] for m in scored_a]
    labels_b = [m["classification"] for m in scored_b]
    assert labels_a == labels_b  # deterministic, run-to-run stable

    n_slotted = sum(1 for lbl in labels_a if is_external_label(lbl))
    assert n_slotted == 6  # n_T(3) + n_R(3) for benzene (non-linear)

    # The degenerate block (EMIT 1-9, eigenvalue -23.7514 cm-1) genuinely
    # contains 2 of the 6 slot-winners (EMIT 6 -> Rx, EMIT 9 -> Ry) and 7
    # modes that fall through to Step 4 with their own vib_label -- these 7
    # are NOT all the same label, which is expected (different eigenvectors,
    # different s[V_S]), not a bug.
    by_name = {m["name"]: m for m in scored_a}
    degenerate_block = [f"EMIT {i}" for i in range(1, 10)]
    freqs = [by_name[n]["frequency"] for n in degenerate_block]
    # Numerically degenerate (same eigenvalue to within the EMIT file's print
    # precision, ~1e-6 relative), not necessarily bit-identical.
    assert max(freqs) - min(freqs) < 1e-4 * abs(freqs[0]), freqs
    block_labels = {n: by_name[n]["classification"] for n in degenerate_block}
    assert is_mixed_external(block_labels["EMIT 6"]) and block_labels["EMIT 6"].startswith("Rx")
    assert is_mixed_external(block_labels["EMIT 9"]) and block_labels["EMIT 9"].startswith("Ry")
    n_external_in_block = sum(1 for lbl in block_labels.values() if is_external_label(lbl))
    assert n_external_in_block == 2  # not 0, not 9 -- exactly the 2 slot-winners


def test_confusion_matrix_precision_perfect_recall_explained_by_mixed_bucket():
    """Clean-category confusion matrix over the full library (fig:confusion's
    numbers), RESTRICTED to the single-centre AB_n hydride-library scope
    (confusion_matrix_stats() applies src.calibrate.SINGLE_CENTRE_ONLY_EXCLUDE
    internally -- see that constant's docstring). Recall is 1.0 for
    translation/rotation (exact completeness); stretch/bend recall is NOT
    quite >=0.95 -- and that is not a defect: most stretch/bend 'misses'
    land in the MIXED bucket, not the opposite clean category.

    **Numbers re-derived 2026-07-08 (T-shaped/see-saw literature back-fill
    session)**: the author manually back-filled `shape`/`type`/`sym`/`ref`
    ground truth in `data/characterised_modes.csv` for 24 molecules that
    previously had NO row at all (the 10 non-ideal T-shaped AB3 families, the
    11 non-ideal see-saw AB4 families, plus `BBr3`/`OCl2`, and -- new
    information this session -- `SnO2`, closing the exact ideal-tier gap the
    PRIOR 2026-07-08 (`data_score.csv` retirement) session's docstring here
    described). `--library`/`--calibrate` were rerun: `thresholds.json`'s
    `tau_S`/`tau_B` are UNCHANGED bit-for-bit (0.9036817451504533 /
    0.17326891344050538) -- SnO2 contributes 2 stretch modes at V_Stretch=1.0
    and 2 bend modes at V_Stretch=0.0, comfortably inside the existing
    boundary, never an extremum (`ideal_stretch_n` 39->41, `ideal_bend_n`
    48->50 in `thresholds.json`, min/max unchanged). `attach_labels()`'s
    frequency-agreement gate had ZERO skip-report entries for all 24 newly
    labeled molecules (172/172 internal rows joined `ref_label` cleanly).

    **Precision is NO LONGER exactly 1.0 for stretch/bend -- a genuinely new
    finding, not a regression in the code:** of the 24 newly-labeled
    molecules' 172 modes, exactly 4 (all in the see-saw AB4 family, and only
    in `OH4`/`OF4` specifically -- no other see-saw or T-shaped molecule
    shows this) land on the OPPOSITE clean category from their literature
    label: `OH4` mode 5 (lit. bend, `V_Stretch`=0.9675, predicted clean
    stretch) and mode 7 (lit. stretch, `V_Stretch`=0.1283, predicted clean
    bend); `OF4` mode 4 (lit. bend, `V_Stretch`=0.9839, predicted clean
    stretch) and mode 6 (lit. stretch, `V_Stretch`=0.0468, predicted clean
    bend) -- all 4 are confidently on the wrong side of tau_S/tau_B, not
    borderline/mixed misses. This is the population's first-ever opposite-
    category crossing (0 cases in every previously-tested non-ideal
    molecule); both are the lightest-ligand members of the see-saw family
    (O central atom with H/F ligands), consistent with a genuine A1-symmetry
    stretch/bend coupling the literature's normal-mode-eigenvector label
    doesn't reflect but the geometric V_Stretch score picks up -- flagged for
    `lead-author`/`expert-reviewer-jcc` to discuss as a concrete limitation
    example, not fixed here (framework classifies; it does not adjudicate
    which of two coupled, same-irrep descriptions is "more correct").
    Precision: stretch 0.98450 (tp=127, n_pred=129, 2 false positives -- the
    2 `OH4`/`OF4` lit.-bend-predicted-stretch modes), bend 0.99010 (tp=200,
    n_pred=202, 2 false positives -- the 2 lit.-stretch-predicted-bend
    modes); translation/rotation stay exactly 1.0.
    """
    lib_df = pd.read_csv(LIB_CSV)
    calibrated = Thresholds.calibrated()
    res = confusion_matrix_stats(lib_df, calibrated, acceptance_floor=0.95)

    for cat in ("translation", "rotation"):
        assert res["per_category"][cat]["precision"] == 1.0, cat

    assert abs(res["per_category"]["stretch"]["precision"] - 0.98450) < 1e-4
    assert abs(res["per_category"]["bend"]["precision"] - 0.99010) < 1e-4

    assert res["per_category"]["translation"]["recall"] == 1.0
    assert res["per_category"]["rotation"]["recall"] == 1.0

    stretch = res["per_category"]["stretch"]
    assert abs(stretch["recall"] - 0.58796) < 1e-3
    assert abs(stretch["mixed_fraction"] - 0.40278) < 1e-3

    bend = res["per_category"]["bend"]
    assert abs(bend["recall"] - 0.89286) < 1e-3
    assert abs(bend["mixed_fraction"] - 0.09821) < 1e-3

    # The floor is NOT met overall, because stretch recall sits well under
    # 0.95 (bend also now falls short) -- reported honestly, not forced to pass.
    assert res["floor_met"] is False


def test_confusion_matrix_ideal_nonideal_recall_split():
    """Task-C addition (2026-07-02): confusion_matrix_stats() now also
    reports recall_ideal/recall_nonideal per category (additive; the pooled
    precision/recall keys above are unchanged). Mirrors
    src/figures.py::plot_confusion_matrix's ad hoc ideal/non-ideal split
    (commit 424a666) as a single formal, tested statistic instead of only
    living inside plotting code.

    recall_ideal['stretch']==1.0 and recall_ideal['bend']==1.0 must hold
    EXACTLY, not approximately -- tau_S/tau_B are derived
    (derive_stretch_bend_thresholds) as literally the min/max of this same
    ideal=='yes' population, so by construction no ideal-tier stretch/bend
    row can land on the wrong side of its own defining boundary. This test
    verifies that construction argument computationally rather than assuming
    it.

    **Numbers re-derived 2026-07-08 (T-shaped/see-saw literature back-fill
    session)** -- see the sibling test's docstring for the full mechanism.
    `n_ref_ideal` for stretch/bend is RESTORED to 41/50 (was 39/48 after the
    prior 2026-07-08 `data_score.csv`-retirement session, which had exposed
    -- but not fixed -- `SnO2` missing its `characterised_modes.csv` row;
    this session's back-fill closes that exact gap, +2 stretch/+2 bend, back
    to the original pre-any-gap size). `recall_ideal` remains EXACTLY 1.0 for
    both, unaffected -- by construction (tau_S/tau_B are literally this
    population's own min/max, so no ideal-tier row can land on the wrong
    side of its own defining boundary; ALSO independently confirmed this
    session: SnO2's 2 stretch modes score V_Stretch=1.0 and its 2 bend modes
    score V_Stretch=0.0, nowhere near tau_S/tau_B, so they were never
    boundary-defining -- `thresholds.json` is bit-identical before/after).
    The non-ideal tier grew substantially this session (23 new non-ideal
    molecules' 168 modes joined `ref_label` for the first time): bend
    recall_nonideal 0.86207 (n_ref_nonideal 174, was 85), stretch
    recall_nonideal 0.49143 (n_ref_nonideal 175, was 96) -- both recall
    figures move because the new T-shaped/see-saw population's mixed/
    opposite-category rate differs from the pre-existing non-ideal pool's
    average, not because of any change to previously-scored molecules'
    numbers (score-neutrality for all pre-existing rows -- unaffected by a
    label-only back-fill).
    """
    lib_df = pd.read_csv(LIB_CSV)
    calibrated = Thresholds.calibrated()
    res = confusion_matrix_stats(lib_df, calibrated, acceptance_floor=0.95)

    assert res["per_category"]["stretch"]["recall_ideal"] == 1.0
    assert res["per_category"]["bend"]["recall_ideal"] == 1.0
    assert res["per_category"]["stretch"]["n_ref_ideal"] == 41
    assert res["per_category"]["bend"]["n_ref_ideal"] == 50

    # Non-ideal tier -- grew this session (see docstring above).
    assert abs(res["per_category"]["bend"]["recall_nonideal"] - 0.86207) < 1e-3
    assert abs(res["per_category"]["stretch"]["recall_nonideal"] - 0.49143) < 1e-3
    assert res["per_category"]["bend"]["n_ref_nonideal"] == 174
    assert res["per_category"]["stretch"]["n_ref_nonideal"] == 175

    # Translation/rotation: every row is an external (T/R) reference, so the
    # ideal tier reproduces the pooled recall exactly and there is no
    # non-ideal tier at all (n_ref_nonideal==0 -> recall_nonideal is NaN).
    for cat in ("translation", "rotation"):
        assert res["per_category"][cat]["recall_ideal"] == 1.0
        assert res["per_category"][cat]["n_ref_nonideal"] == 0
        assert res["per_category"][cat]["recall_nonideal"] != res["per_category"][cat]["recall_nonideal"]  # NaN

    # Pooled keys (existing behavior) must be untouched by this addition --
    # match the sibling test's pooled numbers above.
    assert abs(res["per_category"]["stretch"]["recall"] - 0.58796) < 1e-3
    assert abs(res["per_category"]["bend"]["recall"] - 0.89286) < 1e-3


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
