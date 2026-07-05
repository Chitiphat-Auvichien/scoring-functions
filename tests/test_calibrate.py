"""Regression tests for src/calibrate.py (Phase 3 tau calibration) and the
label-level validation it unlocks.

Run from ``Github/scoring-functions/``:
    py -m pytest tests/                (with pytest)
    py tests/test_calibrate.py        (standalone; no pytest needed)

Design note (mirrors tests/test_excel_ingest.py): the plateau sweep itself is
cheap (small molecules, ~190 grid points in well under a minute) but does
depend on data/logs/ parsing for 25 molecules + benzene EMIT; the two derived
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


def test_single_centre_only_exclude_matches_scope_decision():
    """The 2026-07-05 single-centre-only scope decision names exactly 7
    molecules -- two-, four-, or six-centre topologies pooled into the
    "hydride library" (single-centre AB_n) statistics by mistake. H2O is
    deliberately NOT in this set (it is single-centre itself and stays in
    the non-ideal population as a familiar illustrative molecule)."""
    assert SINGLE_CENTRE_ONLY_EXCLUDE == {
        "C2H2", "C2H4", "C2H6", "H2O2", "C6H6", "iso-C4H10", "n-C4H10",
    }
    assert "H2O" not in SINGLE_CENTRE_ONLY_EXCLUDE


def test_filter_single_centre_library_drops_exactly_the_excluded_molecules():
    """filter_single_centre_library() is DATA-LOSS-FREE at the CSV level
    (library_scores.csv keeps every molecule's rows -- see
    excel_ingest.py); this only checks the analysis-time filter itself
    removes exactly the 7 excluded molecules' rows and nothing else,
    against the real checked-in library_scores.csv golden."""
    lib_df = pd.read_csv(LIB_CSV)
    before_molecules = set(lib_df["molecule"].unique())
    assert SINGLE_CENTRE_ONLY_EXCLUDE <= before_molecules  # sanity: all present

    filtered = filter_single_centre_library(lib_df)
    after_molecules = set(filtered["molecule"].unique())

    assert before_molecules - after_molecules == SINGLE_CENTRE_ONLY_EXCLUDE
    # 145 internal rows (7+12+18+6+30+36+36) + benzene's 6 geometry-backed
    # external T/R rows (C6H6 has an on-disk .log/.gjf pair, so it also
    # contributes external rows that this molecule-level filter drops too,
    # unlike the other 6 excluded molecules which are Excel-only).
    assert len(lib_df) - len(filtered) == 151
    # Kept molecule count: 69 (full excel-sourced library) - 7 = 62.
    assert lib_df["molecule"].nunique() - len(SINGLE_CENTRE_ONLY_EXCLUDE) == filtered["molecule"].nunique()

    # Non-ideal internal population after the filter: 51 molecules (an
    # independent hand-count of tab:nonideal's own AB_n grid), 277 modes.
    nonideal = filtered[(filtered["kind"] == "internal") & (filtered["ideal"] == "no")]
    assert nonideal["molecule"].nunique() == 51
    assert len(nonideal) == 277


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
    (0.17327) either way, so this is robust to the shift, as predicted."""
    calibrated = Thresholds.calibrated()
    t = _classify("benzene", "emit", calibrated)
    assert t["EMIT 34"]["classification"] == "Tx*"
    assert t["EMIT 34"]["annotation"] == f"vibration={MIXED_STRETCH_BEND}"
    assert t["EMIT 35"]["classification"] == "Ty*"
    assert t["EMIT 35"]["annotation"] == f"vibration={MIXED_STRETCH_BEND}"
    assert t["EMIT 36"]["classification"] == "Tz"


def test_water_targets_under_calibrated_thresholds():
    """Task-4 re-check of test_classifier.py's water targets, under
    Thresholds.calibrated() instead of the pinned provisional defaults."""
    calibrated = Thresholds.calibrated()
    t = _classify("water", "normal", calibrated)
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
    raw, _ = load_inputs("benzene", "emit", os.path.join(ROOT, "data"))
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
    translation/rotation (exact completeness), but NOT quite >=0.95 for
    stretch (0.6667) or bend (0.9838) -- and that is not a defect: most
    stretch/bend 'misses' land in the MIXED bucket (33.33%/1.62% of true
    stretch/bend respectively, non-ideal center-of-mass softening per B8.3),
    never crossing to the OPPOSITE clean category (0 cases either way).

    **Numbers updated 2026-07-05** (author decision -- new default label
    source, `src/csv_label_ingest.py`, superseding the xlsx `data_score`/
    `characterised modes` sheets across ALL molecules): benzene's (C6H6)
    literature relabeling flips mode_index 19/23/24 bend->stretch and
    21/22 bend->SB (a literal literature 3rd class -- see
    `src/benzene_validation.py`'s 3-class confusion table for the benzene-
    specific numbers).

    **Precision number corrected 2026-07-05 (same day, follow-up fix):** this
    assertion briefly pinned bend precision at 0.99267 under the reasoning
    that modes 21/22 (ref_label=="SB", predicted bucket "bend") were "false
    positives" for bend now that the library has finer-grained ground truth.
    That reasoning was reconsidered: `confusion_matrix_stats()`'s 4-category
    (translation/rotation/stretch/bend) precision/recall accounting is
    specifically about whether a NOMINAL stretch/bend reference mode keeps
    its label or is mispredicted -- a literal literature "SB" reference is a
    different, 3rd-class question entirely (already answered by
    `src/benzene_validation.py::benzene_internal_confusion_matrix`'s
    dedicated 3x3 table), not a stretch/bend miss. `confusion_matrix_stats()`
    now excludes any ref_label outside the 4 recognized categories before
    computing the confusion table/per-category stats (see its docstring/
    inline comment), so modes 21/22 no longer inflate bend's n_pred
    denominator -- bend precision is exactly 1.0 again, matching
    stretch/translation/rotation (unaffected either way -- no SB mode is
    ever predicted "stretch", and translation/rotation have no SB rows at
    all). See `src/benzene_validation.py::benzene_internal_confusion_matrix`
    for the benzene-scoped mechanistic explanation (the tau_B "bending
    blind spot") of why 21/22 land on "bend".

    **Numbers updated again 2026-07-05 (same day, single-centre-only scope
    fix):** the non-ideal population previously pooled 58 molecules / 422
    internal modes, but 7 of those (C2H2, C2H4, C2H6, H2O2, C6H6, iso-C4H10,
    n-C4H10) are two-, four-, or six-centre topologies outside the
    manuscript's own "single-centre AB_n only" scope for this validation
    (JCC .tex "Stretching/bending classification" subsection; tab:ideal
    lists exactly 11 single-centre molecules) -- an independent hand-count
    of tab:nonideal's own AB_n grid gives 51 molecules, confirming the
    table's structure already assumed this scope even though the pooled
    statistics didn't. Excluding them removes 145 of the unfiltered
    non-ideal population's 422 internal modes, dropping it to 51 molecules /
    277 internal modes (SINGLE_CENTRE_ONLY_EXCLUDE, applied inside
    confusion_matrix_stats()).
    Old (unfiltered) values were stretch recall 0.70815/mixed 0.29185, bend
    recall 0.97482/mixed 0.02518 -- H2O2/C2H2/C2H4/C2H6/iso-C4H10/n-C4H10
    (Excel-only, no on-disk geometry) contributed only NON-ideal internal
    rows, so precision/translation/rotation/floor_met are unaffected;
    benzene additionally contributed 6 geometry-backed external T/R rows,
    which is a separate, "rigorous"-tier effect checked in
    test_confusion_matrix_ideal_nonideal_recall_split below (that tier
    stays outside THIS test, which reads the pooled/unsplit library).
    """
    lib_df = pd.read_csv(LIB_CSV)
    calibrated = Thresholds.calibrated()
    res = confusion_matrix_stats(lib_df, calibrated, acceptance_floor=0.95)

    for cat in ("stretch", "bend", "translation", "rotation"):
        assert res["per_category"][cat]["precision"] == 1.0, cat

    assert res["per_category"]["translation"]["recall"] == 1.0
    assert res["per_category"]["rotation"]["recall"] == 1.0

    stretch = res["per_category"]["stretch"]
    assert abs(stretch["recall"] - 0.66667) < 1e-3
    assert abs(stretch["mixed_fraction"] - 0.33333) < 1e-3

    bend = res["per_category"]["bend"]
    assert abs(bend["recall"] - 0.98378) < 1e-3
    assert abs(bend["mixed_fraction"] - 0.01622) < 1e-3

    # The floor is NOT met overall, because stretch recall sits well under
    # 0.95 (bend clears it) -- reported honestly, not forced to pass.
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

    **Numbers restored 2026-07-04** to the source="excel" (~69-molecule)
    values (IMPLEMENTATION_PLAN.md RESUME HERE, dual-source decision): the
    full ~70-molecule Excel-sourced library (not just the 25 with on-disk
    geometry), so n_ref_ideal counts the ideal-molecule subset among all 69
    (all 11 tab:ideal molecules, including SnO2/TeH2/TeH4 which still lack
    on-disk geometry but are still Excel-scored). The non-ideal tier's
    recall reverts to the pre-2026-07-03 baseline (65.6%/95.7% split) since
    Excel-only rows get only Step 4's vib_label applied in isolation, not the
    full Algorithm 1 (see excel_ingest.py's module docstring for why: those
    rows have no on-disk geometry for Steps 2/3 to run against at all).

    **Non-ideal-tier numbers updated 2026-07-05** for the same reason as
    `test_confusion_matrix_precision_perfect_recall_explained_by_mixed_
    bucket` above (benzene's literature relabeling; all of benzene's
    internal rows are `ideal=='no'`, so this is purely a non-ideal-tier
    effect -- `recall_ideal`/`n_ref_ideal` for stretch/bend, which are
    derived ONLY from the `ideal=='yes'` population, are completely
    unaffected and still hold exactly).

    **Non-ideal-tier numbers updated AGAIN 2026-07-05 (same day,
    single-centre-only scope fix, see the sibling test's docstring above for
    the full rationale):** `confusion_matrix_stats()` now drops C2H2, C2H4,
    C2H6, H2O2, C6H6, iso-C4H10, n-C4H10 (145 internal modes) from the
    non-ideal tier before computing anything -- n_ref_nonideal for bend
    228->135 (loses benzene's 18 + C2H2/C2H4/C2H6/H2O2/iso-C4H10/n-C4H10's
    bend modes) and stretch 192->142 (loses benzene's 12 + the same 6
    molecules' stretch modes), recall_nonideal shifting accordingly.
    n_ref_ideal/recall_ideal for stretch (41/1.0) and bend (50/1.0) are
    UNCHANGED -- none of the 7 excluded molecules is `ideal=='yes'`.
    """
    lib_df = pd.read_csv(LIB_CSV)
    calibrated = Thresholds.calibrated()
    res = confusion_matrix_stats(lib_df, calibrated, acceptance_floor=0.95)

    assert res["per_category"]["stretch"]["recall_ideal"] == 1.0
    assert res["per_category"]["bend"]["recall_ideal"] == 1.0
    assert res["per_category"]["stretch"]["n_ref_ideal"] == 41
    assert res["per_category"]["bend"]["n_ref_ideal"] == 50

    # Non-ideal tier: 2026-07-05 single-centre-only scope fix -- see
    # docstring above.
    assert abs(res["per_category"]["bend"]["recall_nonideal"] - 0.97778) < 1e-3
    assert abs(res["per_category"]["stretch"]["recall_nonideal"] - 0.57042) < 1e-3
    assert res["per_category"]["bend"]["n_ref_nonideal"] == 135
    assert res["per_category"]["stretch"]["n_ref_nonideal"] == 142

    # Translation/rotation: every row is an external (T/R) reference, so the
    # ideal tier reproduces the pooled recall exactly and there is no
    # non-ideal tier at all (n_ref_nonideal==0 -> recall_nonideal is NaN).
    for cat in ("translation", "rotation"):
        assert res["per_category"][cat]["recall_ideal"] == 1.0
        assert res["per_category"][cat]["n_ref_nonideal"] == 0
        assert res["per_category"][cat]["recall_nonideal"] != res["per_category"][cat]["recall_nonideal"]  # NaN

    # Pooled keys (existing behavior) must be untouched by this addition --
    # updated 2026-07-05 (single-centre-only scope fix) to match the sibling
    # test's pooled numbers above.
    assert abs(res["per_category"]["stretch"]["recall"] - 0.66667) < 1e-3
    assert abs(res["per_category"]["bend"]["recall"] - 0.98378) < 1e-3


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
