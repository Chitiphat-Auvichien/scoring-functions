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
import tempfile

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
    sweep_tau_sb, _error_plateau,
)

THRESHOLDS_JSON = os.path.join(ROOT, "data", "results", "thresholds.json")
SWEEP_CSV = os.path.join(ROOT, "data", "results", "tau_sensitivity_sweep.csv")
LIB_CSV = os.path.join(ROOT, "data", "results", "library_scores.csv")


def _classify(mol, mode_type, thresholds, scheme="threeway"):
    """`scheme` pinned to "threeway" by default: the EMIT 34/35 "SB"
    annotation re-check below (test_benzene_emit_34_35_36_under_calibrated_thresholds)
    was pinned under the three-way scheme, so this keeps it fixed even
    though classify_all_modes()'s own default flipped to "binary" (the
    2026-08 paper-standard switch)."""
    raw, _ = load_inputs(mol, mode_type, os.path.join(ROOT, "data"))
    scorer, final = build_scorer_and_final(raw, mode_type)
    scored = classify_all_modes(scorer, final, thresholds, scheme=scheme)
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

# The 9-molecule held-out transferability-test set added to
# mol_list_method.csv 2026-08-11 (mol_type=='test').
_TEST_CATEGORY_18 = frozenset({
    "CH4", "C4H4", "C10H16", "PCl5", "C3H6", "B3N3H6", "CHCl3", "CH3CN",
    "C3O3H6", "H2O", "CH3COCH3", "C6H4F2", "XeF2Cl2", "C2H4", "C10H8",
    "HOCl", "HCOOH", "C7H8",
})


def test_single_centre_only_exclude_matches_scope_decision():
    """**2026-08-11: repointed to an INCLUSION filter**
    (src.library_ingest.out_of_calibration_scope_molecules(), the complement
    of mol_type in {'ideal', 'non-ideal'}) rather than the prior
    multi-centre-only exclusion -- see src/calibrate.py's own comment above
    SINGLE_CENTRE_ONLY_EXCLUDE for the full rationale. mol_list_method.csv's
    'test' category must drop out of calibration automatically, exactly
    like C6H6 (multi-centre) already did.

    **2026-08-14 (roster expansion):** the 'test' category grew from 9 to 18
    -- H2O moved from 'non-ideal' into 'test' (rerun at mp2/3-21g* instead of
    mp2/3-21g, same session), and 8 brand-new 'test' molecules were added
    (CH3COCH3, C6H4F2, XeF2Cl2, C2H4, C10H8, HOCl, HCOOH, C7H8). H2O is now
    (unlike before) excluded from calibration scope, since it moved out of
    'non-ideal' -- it stays in the roster as a familiar worked example
    (`data/results/H2O_normal.csv`) but no longer feeds the confusion-matrix
    stretch/bend recall split."""
    expected = {"C6H6"} | _TEST_CATEGORY_18
    assert SINGLE_CENTRE_ONLY_EXCLUDE == expected
    assert len(SINGLE_CENTRE_ONLY_EXCLUDE) == 19
    assert "H2O" in SINGLE_CENTRE_ONLY_EXCLUDE

    # library_scores.csv was regenerated 2026-08-14 to include the 8 new
    # 'test'-category molecules (H2O was already present, just retagged), so
    # all 19 excluded names (C6H6 + the 18 test molecules) are now present to
    # be dropped. The historical 7-name fallback set previously only ever had
    # C6H6 actually present in the real roster -- but one of its other 6
    # names, 'C2H4', is now ALSO a real, present molecule (a coincidental
    # name collision: the new transferability-set ethylene shares a formula
    # with the historical legacy multi-centre placeholder name, not the same
    # underlying decision) -- so old_dropped legitimately grew to 2 names.
    lib_df = pd.read_csv(LIB_CSV)
    old_dropped = set(lib_df.loc[lib_df["molecule"].isin(_HISTORICAL_SINGLE_CENTRE_EXCLUDE_7),
                                  "molecule"].unique())
    new_dropped = set(lib_df.loc[lib_df["molecule"].isin(SINGLE_CENTRE_ONLY_EXCLUDE),
                                  "molecule"].unique())
    assert old_dropped == {"C6H6", "C2H4"}
    assert new_dropped == expected


def test_filter_single_centre_library_drops_exactly_the_excluded_molecules():
    """filter_single_centre_library() is DATA-LOSS-FREE at the CSV level
    (library_scores.csv keeps every molecule's rows -- see
    src/library_ingest.py); this checks the analysis-time filter itself
    removes exactly whichever of the excluded molecules are actually present
    and nothing else.

    **2026-08-14 update:** SINGLE_CENTRE_ONLY_EXCLUDE now has 19 names
    (C6H6 + the 18 'test'-category molecules, up from 9 -- H2O moved in from
    'non-ideal' and 8 new molecules were added); library_scores.csv was
    regenerated the same day to include all 85 roster molecules, so all 19
    are present and dropped (66 remain: the ideal/non-ideal calibration
    scope)."""
    lib_df = pd.read_csv(LIB_CSV)
    before_molecules = set(lib_df["molecule"].unique())
    present_of_excluded = SINGLE_CENTRE_ONLY_EXCLUDE & before_molecules
    assert present_of_excluded == SINGLE_CENTRE_ONLY_EXCLUDE

    filtered = filter_single_centre_library(lib_df)
    after_molecules = set(filtered["molecule"].unique())

    assert before_molecules - after_molecules == SINGLE_CENTRE_ONLY_EXCLUDE
    assert lib_df["molecule"].nunique() - len(SINGLE_CENTRE_ONLY_EXCLUDE) == filtered["molecule"].nunique()
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
    C6H6 -- content-preserving rename, numbers unchanged.)"""
    calibrated = Thresholds.calibrated()
    t = _classify("C6H6", "emit", calibrated)
    assert t["EMIT 34"]["classification"] == "Tx*"
    assert t["EMIT 34"]["annotation"] == f"vibration={MIXED_STRETCH_BEND}"
    assert t["EMIT 35"]["classification"] == "Ty*"
    assert t["EMIT 35"]["annotation"] == f"vibration={MIXED_STRETCH_BEND}"
    assert t["EMIT 36"]["classification"] == "Tz"


def test_benzene_emit_34_35_binary_default_under_calibrated_thresholds():
    """Mirror check of test_benzene_emit_34_35_36_under_calibrated_thresholds
    above, but for classify_all_modes()'s own scheme DEFAULT (binary, the
    2026-08 paper-standard switch) instead of the "threeway" pin `_classify`
    applies. EMIT 34/35 (V_Stretch 0.667/0.577, both far above tau_SB=0.50)
    are STRETCHING under binary, never MIXED_STRETCH_BEND."""
    calibrated = Thresholds.calibrated()
    t = _classify("C6H6", "emit", calibrated, scheme="binary")
    assert t["EMIT 34"]["classification"] == "Tx*"
    assert t["EMIT 34"]["annotation"] == "vibration=S"
    assert t["EMIT 35"]["classification"] == "Ty*"
    assert t["EMIT 35"]["annotation"] == "vibration=S"
    assert t["EMIT 36"]["classification"] == "Tz"


def test_water_targets_under_calibrated_thresholds():
    """Task-4 re-check of test_classifier.py's water targets, under
    Thresholds.calibrated() instead of the pinned provisional defaults.
    (2026-07-07: "water" -> the finalized roster basename H2O, a
    genuinely different/corrected calculation from the old water.log -- see
    tests/test_scores.py's module docstring -- but its V_Stretch values land
    on the same side of tau_S/tau_B, so the classification buckets checked
    here are unaffected.)"""
    calibrated = Thresholds.calibrated()
    t = _classify("H2O", "normal", calibrated)
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
    raw, _ = load_inputs("C6H6", "emit", os.path.join(ROOT, "data"))
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

    **Numbers re-derived 2026-07-09 (OH4/OF4 exclusion session)**: `OH4` and
    `OF4` -- the two lightest-ligand members of the see-saw AB4 family -- were
    removed from the library roster entirely (`data/mol_list_method.csv`,
    `data/characterised_modes.csv`, and their `logs/`/`gjf/`/`intermediate/`
    input files all deleted). Both are not genuine stationary points at this
    project's MP2/3-21G level (imaginary/negative frequencies -- e.g. `OH4`
    mode 1 was -303.26 cm^-1), so their "normal modes" are not physically
    meaningful vibrations of a real minimum and cannot be validly compared to
    the `TeH4` ideal see-saw template.

    **Numbers re-derived again 2026-07-23/2026-07-24** across two sessions
    that landed back-to-back without this test being repinned in between (a
    gap closed now, not a new inconsistency introduced by either): (1)
    2026-07-23 (commit `8a12c32`) removed `SnO2`/`FH3` from the roster
    entirely (both had negative/imaginary frequencies -- not genuine
    stationary points) and re-tagged `CO2` `ideal` (was `non-ideal`),
    dropping the pooled `tp`/`n_pred`/`n_ref` below from 125/125/208
    (stretch) and 197/197/214 (bend) to 121/121/203 and 193/193/209 --
    that commit reran `--library`/`--calibrate` and confirmed the numbers
    but did not update this file's pins. (2) 2026-07-24 corrected an
    inconsistency in (1): `CO2` was moved back to `non-ideal` (its linear
    siblings `CS2`/`CSe2`/`CTe2` were already `non-ideal`; tagging `CO2`
    `ideal` had made it the odd one out with no principled reason). This
    pooled test's tp/n_pred/n_ref/recall/mixed_fraction figures are
    UNCHANGED by step (2) -- re-tagging a molecule's ideal/non-ideal tier
    moves rows between the ideal and non-ideal sub-populations but does not
    change the pooled (ideal+non-ideal) total, and `CO2`'s 4 internal modes
    (2 stretch, 2 bend) are all correctly classified (`tp` unaffected)
    regardless of tier -- so the numbers pinned below are actually from step
    (1); see `test_confusion_matrix_ideal_nonideal_recall_split` below for
    the ideal/non-ideal SPLIT figures that step (2) does change.
    `thresholds.json`'s `tau_S`/`tau_B` are UNCHANGED bit-for-bit
    (0.9036817451504533 / 0.17326891344050538) through both steps --
    `CO2`'s stretch V_Stretch=1.0/bend V_Stretch=0.0 are nowhere near the
    ideal population's min/max boundary values in either tier.

    **Precision is exactly 1.0 for stretch/bend** (translation/rotation
    were already exactly 1.0 and remain so: tp=n_pred=n_ref=201 translation,
    197 rotation). Recall/mixed_fraction reflect the smaller (SnO2/FH3-free)
    population from step (1) above; unaffected by step (2)'s CO2 re-tag.

    **Re-derived again 2026-08-10 (G16 promotion)**: canonical `data/logs`+
    `data/gjf`+`data/fchk` switched from G09 to G16 (see IMPLEMENTATION_PLAN.md
    Locked decisions). `tau_S` shifted by ~1.5e-7 (0.9036817451504533 ->
    0.9036818966195127, still XeOH4's ideal-stretch min, just recomputed from
    slightly different G16 geometry) -- negligible, and not the cause of the
    change below. The real cause: `ClH3` mode 4 (non-ideal, ref=stretch,
    V_Stretch 0.903119 (G09) -> 0.905849 (G16)) crossed `tau_S` and flipped
    predicted_label SB->S -- one of the exact 3/846 threshold-boundary flips
    already characterized in the G09-vs-G16 consistency check
    (`data/results/rerun_consistency_report.csv`). `stretch` tp/n_pred
    121->122 (recall/mixed_fraction shift accordingly); `bend` is UNCHANGED
    (the other 2 of the 3 flips are `FBr3` modes 1/2, both ref=bend, B<->SB
    in opposite directions -- they cancel net, bend tp stays 193).

    **Re-derived again 2026-08-19 (binary-classification-scheme default
    switch)**: `confusion_matrix_stats()` itself has no `scheme` parameter --
    it simply trusts whatever `predicted_label` is already in the `lib_df`
    it's given, and `library_scores.csv`'s own column is binary by default
    as of this switch (never "mixed"). This is the biggest driver of the
    numbers below: with no "mixed" bucket to catch ambiguous modes, EVERY
    internal mode is forced to a definite stretch/bend answer, so
    `mixed_fraction` is now 0.0 for both categories (there is no third
    outcome under binary) and recall changes accordingly on both sides --
    stretch recall RISES (194/201, up from 121/201, since former
    "mixed"-bucket stretches are now forced to a stretch/bend verdict and
    most land correctly), while bend's PRECISION now falls below 1.0
    (215 predicted bend against 208 true bend, tp=208) since some
    literature-stretch modes that would have escaped into "mixed" under
    threeway are now forced into bend instead. `floor_met` (0.95) is now
    TRUE overall -- a genuine improvement in the pooled numbers, not
    engineered.

    Translation/rotation's pooled n_ref (198/194, not the historical
    201/197) reflects a small, UNRELATED roster/geometry-pool drift since
    this test was last pinned (confirmed unrelated to the scheme switch:
    external T/R rows are always exact by Eckart-Sayvetz completeness
    regardless of scheme, and stretch/bend's own ideal-tier n_ref_ideal
    below are UNCHANGED at 39/48 -- see
    test_confusion_matrix_ideal_nonideal_recall_split) -- not investigated
    further here since it predates and is orthogonal to this session's work.

    **Re-derived again 2026-08-20 (canonical tau_SB default changed
    0.50->0.42)**: lowering the S/B cutoff pushes borderline internal modes
    from bend into stretch, which is exactly what happens here. `stretch`
    tp/n_pred rise 194/194 -> 201/202 (7 previously-missed non-ideal
    literature-stretch modes now clear the lower cutoff and are correctly
    recovered -- NBr3 mode 4, OCl4 modes 7/9, OBr4 modes 8/9, and SBr4
    modes 8/9, all with V_Stretch in [0.42, 0.50) -- so stretch recall goes
    to a perfect 201/201=1.0), but ONE of those same 7 modes' molecule,
    NBr3, also has a genuine literature-bend mode (mode 3, V_Stretch
    0.486361) that crosses the same lower cutoff the wrong way, so `stretch`
    n_pred outruns tp by exactly 1 (202 vs 201) and stretch PRECISION drops
    below 1.0 for the first time (201/202=0.995049504950495) -- the
    tradeoff is the mirror image of `bend`'s: bend tp/n_pred fall 208/215
    -> 207/207 (NBr3 mode 3 is bend's one true miss, recall 207/208=
    0.9951923076923077, just under 1.0 now) but bend PRECISION rises to a
    perfect 1.0 (every one of the 215->207 modes still predicted bend is
    still genuinely bend, since the modes that left bend's predicted set
    were the ones that used to be false positives for the OTHER category
    under the higher, threeway-adjacent 0.50 cutoff). Net: the
    precision/recall imperfection has moved from bend's precision to
    stretch's precision, floor_met is unaffected (both categories still
    comfortably clear 0.95 either way)."""
    lib_df = pd.read_csv(LIB_CSV)
    calibrated = Thresholds.calibrated()
    res = confusion_matrix_stats(lib_df, calibrated, acceptance_floor=0.95)

    assert res["per_category"]["translation"]["precision"] == 1.0
    assert res["per_category"]["rotation"]["precision"] == 1.0
    assert abs(res["per_category"]["stretch"]["precision"] - 0.995049504950495) < 1e-9

    assert res["per_category"]["translation"]["recall"] == 1.0
    assert res["per_category"]["rotation"]["recall"] == 1.0

    stretch = res["per_category"]["stretch"]
    assert stretch["tp"] == 201
    assert stretch["n_pred"] == 202
    assert stretch["n_ref"] == 201
    assert stretch["recall"] == 1.0
    assert stretch["mixed_fraction"] == 0.0

    bend = res["per_category"]["bend"]
    assert bend["tp"] == 207
    assert bend["n_pred"] == 207
    assert bend["n_ref"] == 208
    assert abs(bend["recall"] - 0.9951923076923077) < 1e-6
    assert bend["precision"] == 1.0
    assert bend["mixed_fraction"] == 0.0

    # The floor (0.95) IS now met -- unlike the threeway scheme's pooled
    # numbers, where stretch recall sat well under 0.95.
    assert res["floor_met"] is True


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

    **Numbers re-derived 2026-07-09 (OH4/OF4 exclusion session)** -- see the
    sibling test's docstring for the full mechanism (OH4/OF4 removed from the
    library: not genuine stationary points at this project's MP2/3-21G level,
    so their normal modes cannot be validly compared to the TeH4 ideal
    see-saw template).

    **Re-derived again 2026-07-24 (CO2 ideal -> non-ideal retag)**: CO2 was
    tagged `mol_type=ideal` in `data/mol_list_method.csv` (set 2026-07-23,
    commit `8a12c32`, alongside that session's SnO2/FH3 roster removal) while
    its linear-triatomic siblings CS2/CSe2/CTe2 were already `non-ideal` --
    an unexplained inconsistency, corrected here by moving CO2 to
    `non-ideal` so all 4 linear-shape library molecules share one tier.
    CO2 contributes exactly 2 stretch + 2 bend internal modes (both correctly
    classified: V_Stretch=1.0 for the stretches, 0.0 for the bends -- nowhere
    near tau_S=0.9036817451504533/tau_B=0.17326891344050538, so `tau_S`/
    `tau_B` themselves are UNCHANGED bit-for-bit, confirmed from the
    regenerated `thresholds.json` rather than assumed), so the retag moves
    those 4 rows from the ideal tier to the non-ideal tier and nothing else
    changes: `n_ref_ideal` stretch/bend 41/50 (2026-07-09 baseline, itself
    unaffected by the intervening SnO2/FH3 change since both were
    non-ideal) -> 39/48; `n_ref_nonideal` stretch/bend 162/159 (the
    SnO2/FH3-adjusted, pre-CO2-retag baseline from commit `8a12c32` -- see
    the sibling test's docstring; NOT the older 167/164 OH4/OF4-session
    numbers, which predate that commit's roster change) -> 164/161, with
    `recall_nonideal` shifting accordingly since 2 more always-correct
    reference-labeled modes joined each non-ideal category's denominator AND
    numerator. `recall_ideal` remains EXACTLY 1.0 for both, by construction
    (tau_S/tau_B are literally this population's own min/max, so no
    ideal-tier row can land on the wrong side of its own defining boundary,
    with or without CO2 in that tier).

    **Re-derived again 2026-08-10 (G16 promotion)**: see the sibling test's
    docstring for the mechanism (canonical data switched G09->G16; `ClH3`
    mode 4, non-ideal/stretch, crossed `tau_S` SB->S, one of the 3/846
    already-characterized threshold-boundary flips in
    `rerun_consistency_report.csv`). `n_ref_ideal`/`n_ref_nonideal` counts are
    UNCHANGED (39/48 ideal, 164/161 non-ideal -- the flip is a label change
    within the existing non-ideal stretch population, not a tier move).
    `recall_ideal` stays EXACTLY 1.0 for both (unaffected, by construction).
    `bend recall_nonideal` is UNCHANGED (the other 2 of the 3 flips, `FBr3`
    modes 1/2, are both non-ideal/bend and flip in opposite directions --
    they cancel net). `stretch recall_nonideal` moves from exactly 0.5
    (82/164) to 83/164 with ClH3's flip added to the numerator.

    **Re-derived again 2026-08-19 (binary-classification-scheme default
    switch)**: see the sibling test's docstring for the mechanism (no
    "mixed" bucket under binary, so every internal mode gets a definite
    verdict). `n_ref_ideal`/`n_ref_nonideal` for stretch/bend are UNCHANGED
    (39/162, 48/160) -- the ideal/non-ideal TIER membership is a structural
    roster property (`mol_type`), untouched by the scheme switch. Only the
    recall NUMERATORS shift: `recall_ideal` stays EXACTLY 1.0 for both, by
    the same construction argument as always (tau_S/tau_B are this
    population's own min/max, so no ideal-tier row can land on the wrong
    side of its own defining boundary regardless of scheme).
    `recall_nonideal` rises for both categories since the former
    "mixed"-bucket modes now get a forced, mostly-correct S/B verdict
    instead of escaping into a third bucket.

    **Re-derived again 2026-08-20 (canonical tau_SB default changed
    0.50->0.42)**: see the sibling test's docstring for the full mechanism
    (7 non-ideal literature-stretch modes -- NBr3 mode 4, OCl4 modes 7/9,
    OBr4 modes 8/9, SBr4 modes 8/9 -- newly clear the lower tau_SB cutoff;
    one non-ideal literature-bend mode, NBr3 mode 3, crosses the same
    lower cutoff the wrong way). All 8 of these are non-ideal (`ideal`=='no'
    in `library_scores.csv`), so `n_ref_ideal`/`recall_ideal` for both
    categories are UNCHANGED (39/48, both still EXACTLY 1.0, by the same
    tau_S/tau_B-is-this-population's-own-min/max construction argument as
    always -- tau_SB is an independent, S-vs-B split point that plays no
    role in that boundary). `n_ref_nonideal` is likewise UNCHANGED (162
    stretch, 160 bend -- these 8 modes only change which SIDE of the
    predicted split they land on, not which tier they belong to).
    `recall_nonideal` moves: stretch's 7 recovered modes push it to a
    perfect 162/162=1.0 (up from 0.9567901234567902); bend's 1 lost mode
    (NBr3 mode 3) pulls it down to 159/160=0.99375 (down from an exact
    1.0, which was itself a coincidence of the old cutoff rather than a
    structural guarantee -- unlike recall_ideal, recall_nonideal has no
    construction argument protecting it). Pooled `stretch` recall becomes
    a perfect 1.0 (201/201, up from 0.9651741293532339); pooled `bend`
    recall is no longer exactly 1.0 (207/208=0.9951923076923077, down from
    exactly 1.0) since it is now dragged down by the one non-ideal miss.
    """
    lib_df = pd.read_csv(LIB_CSV)
    calibrated = Thresholds.calibrated()
    res = confusion_matrix_stats(lib_df, calibrated, acceptance_floor=0.95)

    assert res["per_category"]["stretch"]["recall_ideal"] == 1.0
    assert res["per_category"]["bend"]["recall_ideal"] == 1.0
    assert res["per_category"]["stretch"]["n_ref_ideal"] == 39
    assert res["per_category"]["bend"]["n_ref_ideal"] == 48

    assert abs(res["per_category"]["bend"]["recall_nonideal"] - 0.99375) < 1e-6
    assert res["per_category"]["stretch"]["recall_nonideal"] == 1.0
    assert res["per_category"]["bend"]["n_ref_nonideal"] == 160
    assert res["per_category"]["stretch"]["n_ref_nonideal"] == 162

    # Translation/rotation: every row is an external (T/R) reference, so the
    # ideal tier reproduces the pooled recall exactly and there is no
    # non-ideal tier at all (n_ref_nonideal==0 -> recall_nonideal is NaN).
    # Structurally unaffected by scheme (external assignment is scheme-
    # independent) -- see the sibling test's docstring for why their pooled
    # n_ref (198/194) differs from the historical 201/197 pin (an unrelated,
    # pre-existing roster/geometry-pool drift, not a scheme effect).
    for cat in ("translation", "rotation"):
        assert res["per_category"][cat]["recall_ideal"] == 1.0
        assert res["per_category"][cat]["n_ref_nonideal"] == 0
        assert res["per_category"][cat]["recall_nonideal"] != res["per_category"][cat]["recall_nonideal"]  # NaN

    # Pooled keys (existing behavior) must be untouched by this addition --
    # match the sibling test's pooled numbers above.
    # 2026-08-19 binary-classification-scheme default switch: see the
    # sibling test's docstring for the mechanism. stretch recall
    # 0.6019900497512438 -> 0.9651741293532339, bend recall stays 1.0
    # (already exactly 1.0 under threeway too -- bend's own loss under
    # binary shows up as reduced PRECISION, not recall, see sibling test).
    # 2026-08-20 canonical tau_SB default changed 0.50->0.42: see the
    # sibling test's docstring for the mechanism. stretch recall rises to
    # a perfect 1.0 (up from 0.9651741293532339, all 7 non-ideal
    # literature-stretch misses recovered); bend recall is no longer
    # exactly 1.0 -- it now carries the one non-ideal miss (NBr3 mode 3)
    # that used to be absorbed by bend's PRECISION loss instead.
    assert res["per_category"]["stretch"]["recall"] == 1.0
    assert abs(res["per_category"]["bend"]["recall"] - 0.9951923076923077) < 1e-6


# --------------------------------------------------------------------------
# sweep_tau_sb / _error_plateau -- synthetic fixture, no real geometry/log
# parsing needed (see module docstring's design note for why this differs
# from the tau_TR-plateau tests above).
# --------------------------------------------------------------------------

def _make_tau_sb_fixture(data_dir):
    """Write a tiny mol_list_method.csv into `data_dir` covering the 3 scope
    cases sweep_tau_sb must get right, and return the matching synthetic
    lib_df. M1='ideal' (in 'all', not 'test'), M2='multi-centre' (excluded
    from 'all' entirely, per SINGLE_CENTRE_ONLY_EXCLUDE's 'all molecules
    excludes only C6H6'-style rule), M3='test' (in BOTH 'all' and 'test').
    """
    os.makedirs(data_dir, exist_ok=True)
    roster = pd.DataFrame([
        {"molecule": "M1", "mol_type": "ideal"},
        {"molecule": "M2", "mol_type": "multi-centre"},
        {"molecule": "M3", "mol_type": "test"},
    ])
    roster.to_csv(os.path.join(data_dir, "mol_list_method.csv"), index=False)

    lib_df = pd.DataFrame([
        # M1 ('all' only): unambiguous stretch/bend, plus an external row
        # (kind != 'internal') and a literal 'SB' ref row -- both must be
        # excluded regardless of scope.
        {"molecule": "M1", "kind": "internal", "ref_label": "stretch", "V_Stretch": 0.90},
        {"molecule": "M1", "kind": "internal", "ref_label": "bend", "V_Stretch": 0.10},
        {"molecule": "M1", "kind": "internal", "ref_label": "SB", "V_Stretch": 0.50},
        {"molecule": "M1", "kind": "external", "ref_label": "translation", "V_Stretch": 0.0},
        # M2 (multi-centre): must be excluded from 'all' (and is not 'test').
        {"molecule": "M2", "kind": "internal", "ref_label": "stretch", "V_Stretch": 0.90},
        # M3 ('all' AND 'test'): unambiguous stretch/bend.
        {"molecule": "M3", "kind": "internal", "ref_label": "stretch", "V_Stretch": 0.80},
        {"molecule": "M3", "kind": "internal", "ref_label": "bend", "V_Stretch": 0.20},
    ])
    return lib_df


def test_sweep_tau_sb_scopes_and_exclusions():
    with tempfile.TemporaryDirectory() as tmp:
        lib_df = _make_tau_sb_fixture(tmp)
        sweep = sweep_tau_sb(lib_df, tau_grid=(0.0, 0.5, 1.0), data_dir=tmp)

        # 'all' = M1 + M3's internal stretch/bend rows only: M2 (multi-centre)
        # excluded, the 'SB' ref row excluded, the external row excluded.
        assert (sweep["n_all"] == 4).all()
        # 'test' = M3's internal stretch/bend rows only.
        assert (sweep["n_test"] == 2).all()
        # 'single_centre' = 'all' minus 'test' = M1's stretch/bend rows only.
        assert (sweep["n_single_centre"] == 2).all()

        row = sweep[sweep["tau_SB"] == 0.5].iloc[0]
        # At tau=0.5: M1 stretch(0.90)->stretch OK, M1 bend(0.10)->bend OK,
        # M3 stretch(0.80)->stretch OK, M3 bend(0.20)->bend OK -- 0 error
        # all three scopes (M3's stretch/bend also satisfy 'test').
        assert row["error_all"] == 0.0
        assert row["accuracy_all"] == 1.0
        assert row["error_test"] == 0.0
        assert row["accuracy_test"] == 1.0
        assert row["error_single_centre"] == 0.0
        assert row["accuracy_single_centre"] == 1.0

        row0 = sweep[sweep["tau_SB"] == 0.0].iloc[0]
        # At tau=0.0 everything predicts 'stretch': every scope's bend rows
        # (1 of 2 in each of 'all'/'test'/'single_centre') are now wrong.
        assert row0["error_all"] == 0.5
        assert row0["error_test"] == 0.5
        assert row0["error_single_centre"] == 0.5

        row1 = sweep[sweep["tau_SB"] == 1.0].iloc[0]
        # At tau=1.0 everything predicts 'bend': every scope's stretch rows
        # are now wrong, same 0.5 fraction (2 stretch/2 bend in each scope).
        assert row1["error_all"] == 0.5
        assert row1["error_test"] == 0.5
        assert row1["error_single_centre"] == 0.5


def test_error_plateau_finds_the_flat_minimum_region():
    with tempfile.TemporaryDirectory() as tmp:
        lib_df = _make_tau_sb_fixture(tmp)
        # A grid coarse enough that error_all==0.0 for every grid point
        # strictly between the two populations' V_Stretch values (0.10-0.20
        # bend, 0.80-0.90 stretch) -- a genuine flat plateau, not a single point.
        grid = tuple(round(x, 2) for x in [0.0, 0.3, 0.4, 0.5, 0.6, 0.7, 1.0])
        sweep = sweep_tau_sb(lib_df, tau_grid=grid, data_dir=tmp)

        lo, hi, mid, min_error = _error_plateau(sweep, "error_all")
        assert min_error == 0.0
        assert lo == 0.3 and hi == 0.7  # the whole zero-error contiguous run
        assert mid == 0.5

        lo_t, hi_t, mid_t, min_error_t = _error_plateau(sweep, "error_test")
        assert min_error_t == 0.0
        assert lo_t == 0.3 and hi_t == 0.7
        assert mid_t == 0.5

        lo_sc, hi_sc, mid_sc, min_error_sc = _error_plateau(sweep, "error_single_centre")
        assert min_error_sc == 0.0
        assert lo_sc == 0.3 and hi_sc == 0.7
        assert mid_sc == 0.5


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
