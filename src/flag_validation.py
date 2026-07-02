"""Phase-6 (recommended-before-submission) systematic flag precision/recall
validation -- IMPLEMENTATION_PLAN.md Phase 6, item 1: "Flag precision/recall
over ALL 36 benzene EMIT modes (and library externals) against the C2
projection reference." Promotes the prior anecdotal EMIT 2/9/34-36 spot-check
into a full 36-mode confusion count, since benzene EMIT is this manuscript's
only remaining stress test of the classifier's MIXED_EXTERNAL_WITH_VIBRATION
flag (Gramicidin's scale demonstration was deferred to the companion paper,
Decision 5).

What is being compared
-----------------------
- classify_all_modes()'s PREDICTED flag: does Algorithm 1's two-gate purity
  test (src/classifier.py Step 3) label a mode MIXED_EXTERNAL_WITH_VIBRATION
  ("predicted positive") or not ("predicted negative" -- this collapses
  CLEAN_TRANSLATION/CLEAN_ROTATION/STRETCHING/BENDING/MIXED_STRETCH_BEND into
  one bucket, since a mode never even assigned an external slot in Step 2 has
  no opportunity to be flagged at all -- see the FN-mechanism note below).
- The projection's (src/projection.py, eq:emitproj) GROUND TRUTH: does the
  mode genuinely have fractional external+vibration character?

Ground-truth criterion (no single canonical one exists a priori for this --
unlike the calibration's exact-by-construction normal-mode T/R ground truth;
this is the same honesty the task itself calls for). Define, per EMIT mode,

    M_ext = max(C2_Tx, C2_Ty, C2_Tz, C2_Rx, C2_Ry, C2_Rz)   (projection fractions,
                                                              benzene_EMIT_contributions.csv)

    ground truth = MIXED  iff  GT_EXT_LO < M_ext < GT_EXT_HI   (genuine fractional
                                                                 external/vibration split)
    ground truth = CLEAN  iff  M_ext <= GT_EXT_LO  or  M_ext >= GT_EXT_HI
                                                        (negligible external character,
                                                         or ~complete external dominance)

Thresholds: GT_EXT_HI = 0.95 mirrors tau_TR itself (the classifier's own
Step-3 purity bar) -- reusing the same number keeps "dominant" consistently
defined across ground truth and classifier. GT_EXT_LO = 0.05 is a small,
symmetric "negligible" floor, roughly 10-50x the ~1e-4-1e-3 numerical
orthonormality residual documented in projection.py's module docstring, so
real (if small) coupling is not misclassified as noise. Verified empirically
this session: no benzene EMIT mode's M_ext exceeds ~0.7674 (the EMIT 34/35/36
triad) -- the GT_EXT_HI=0.95 "clean-external, M_ext~1" branch never actually
fires for these 36 modes (EMIT eigenvectors are not constructed to be pure
external modes the way real geometry-backed normal-mode T/R references are;
see library_external_flag_confusion for the case where M_ext DOES hit 1.0
exactly, by construction). All 19 "ground truth CLEAN" modes here are
CLEAN via the M_ext<=0.05 (purely internal) branch, not the M_ext>=0.95
branch -- documented, not hidden.

Result (see IMPLEMENTATION_PLAN.md Changelog / final report for the full
write-up): precision = 1.0 (the flag never fires on a genuinely clean mode),
recall = 5/17 = 0.294 (it misses most genuinely mixed modes). The
false-negative mechanism generalizes the previously-documented EMIT-36 blind
spot (Decision X, "amplitude-invariant, score-indistinguishable-from-pure-
translation") to a SECOND, independent, and more widespread cause: Step 2's
plain one-to-one linear_sum_assignment only ever assigns exactly n_T+n_R=6
of the 36 modes an external slot at all (by construction -- the algorithm
spec's global assignment, not a bug); the other 30 modes fall straight to
Step 4 and can NEVER be flagged MIXED_EXTERNAL_WITH_VIBRATION regardless of
how much genuine external character their projection shows (e.g. EMIT 1, 2,
5, 7, 8, 10-14, 18 all have 7-39% external character by projection but are
Step-2 assignment "losers" for their slot, not just amplitude-degenerate
like EMIT 36). EMIT 2 vs EMIT 9 is the sharpest illustration: EMIT 2 has
MORE genuine external (Ry) character by projection (38.7%) than EMIT 9
(14.1%), yet EMIT 9 wins the Ry slot (its s[Ry] SCORE is larger, 0.215 vs
0.143 -- the documented score/projection ranking inversion), so EMIT 9 is
correctly flagged (TP) while EMIT 2, despite being MORE mixed, is missed
(FN). This is evidence the low recall is systematic to Step 2's one-to-one
assignment mechanism, not an isolated edge case.
"""
import os

import numpy as np
import pandas as pd

from src.classifier import Thresholds, MIXED_EXTERNAL_WITH_VIBRATION

# Ground-truth thresholds on the projection-derived max external fraction.
GT_EXT_LO = 0.05
GT_EXT_HI = 0.95

_EXTERNAL_C2_COLS = ("C2_Tx", "C2_Ty", "C2_Tz", "C2_Rx", "C2_Ry", "C2_Rz")


def _confusion_from_bools(pred_pos, gt_pos):
    """Standard 2x2 confusion-matrix stats from two boolean arrays."""
    pred_pos = np.asarray(pred_pos, dtype=bool)
    gt_pos = np.asarray(gt_pos, dtype=bool)
    tp = int(np.sum(pred_pos & gt_pos))
    fp = int(np.sum(pred_pos & ~gt_pos))
    fn = int(np.sum(~pred_pos & gt_pos))
    tn = int(np.sum(~pred_pos & ~gt_pos))
    precision = tp / (tp + fp) if (tp + fp) else float("nan")
    recall = tp / (tp + fn) if (tp + fn) else float("nan")
    return {"TP": tp, "FP": fp, "FN": fn, "TN": tn,
            "precision": precision, "recall": recall, "n": tp + fp + fn + tn}


def ground_truth_label(m_ext, ext_lo=GT_EXT_LO, ext_hi=GT_EXT_HI):
    """MIXED iff ext_lo < m_ext < ext_hi, else CLEAN. See module docstring."""
    return "MIXED" if (ext_lo < m_ext < ext_hi) else "CLEAN"


def benzene_emit_flag_confusion(data_dir="data", thresholds=None,
                                 ext_lo=GT_EXT_LO, ext_hi=GT_EXT_HI):
    """Systematic per-mode flag confusion over all 36 benzene EMIT modes.

    Recomputes both the classifier labels (classify_all_modes) and the
    projection fractions (project_emit) fresh from data/logs + data/EMIT,
    rather than reading the committed CSVs, so this stays correct if either
    upstream module changes without a matching CSV regeneration.

    Returns (detail_df, stats):
      detail_df -- one row per EMIT mode: Mode, classifier_label, M_ext,
        ground_truth ('MIXED'/'CLEAN'), predicted_positive,
        ground_truth_positive, cell ('TP'/'FP'/'FN'/'TN').
      stats -- _confusion_from_bools() dict + 'ext_lo'/'ext_hi'.
    """
    from main import load_inputs, build_scorer_and_final, run_projection_pipeline
    from src.classifier import classify_all_modes

    thresholds = thresholds or Thresholds.calibrated()

    raw, _ = load_inputs("benzene", "emit", data_dir)
    scorer, final = build_scorer_and_final(raw, "emit")
    scored = classify_all_modes(scorer, final, thresholds)
    by_name = {m["name"]: m for m in scored}

    df_contrib, _, _ = run_projection_pipeline("benzene", data_dir, thresholds, write=False)
    df_contrib = df_contrib.copy()
    df_contrib["M_ext"] = df_contrib[list(_EXTERNAL_C2_COLS)].max(axis=1)

    rows = []
    for _, r in df_contrib.iterrows():
        name = r["Mode"]
        m = by_name[name]
        m_ext = float(r["M_ext"])
        gt_pos = ext_lo < m_ext < ext_hi
        pred_pos = m["classification"] == MIXED_EXTERNAL_WITH_VIBRATION
        if pred_pos and gt_pos:
            cell = "TP"
        elif pred_pos and not gt_pos:
            cell = "FP"
        elif gt_pos and not pred_pos:
            cell = "FN"
        else:
            cell = "TN"
        rows.append({
            "Mode": name,
            "classifier_label": m["classification"],
            "M_ext": m_ext,
            "ground_truth": "MIXED" if gt_pos else "CLEAN",
            "predicted_positive": pred_pos,
            "ground_truth_positive": gt_pos,
            "cell": cell,
        })
    detail_df = pd.DataFrame(rows)
    stats = _confusion_from_bools(detail_df["predicted_positive"], detail_df["ground_truth_positive"])
    stats["ext_lo"], stats["ext_hi"] = ext_lo, ext_hi
    return detail_df, stats


def library_external_flag_confusion(data_dir="data", lib_df=None):
    """Systematic clean-vs-mixed-external check on the 25 geometry-backed
    library molecules' REAL normal-mode T/R references (the task's
    parenthetical "and library externals").

    Ground truth here is degenerate by construction, not merely assumed: a
    real normal mode's ideal T/R reference IS one of the projection
    reference basis Q's own columns (src/projection.py), so it is
    orthonormal to the other 3N-1 columns to machine/print precision --
    Eckart-Sayvetz completeness -- meaning M_ext=1.0 exactly for every one of
    these rows without needing to actually run project_emit on them. Ground
    truth is therefore CLEAN for 100% of these rows; this function verifies
    (rather than assumes) that the classifier's predicted_label is never
    MIXED_EXTERNAL_WITH_VIBRATION for any of them (FP=0) -- the same fact
    src/calibrate.py's tau_TR sensitivity sweep already established
    indirectly (100% accuracy at every grid point), re-expressed here as a
    directly comparable flag-confusion count alongside the benzene EMIT
    table above.

    Returns (ext_df, stats) where ext_df is the filtered
    (kind=='external' & has_geometry) subset of library_scores.csv and stats
    is the _confusion_from_bools() dict (TP=FN=0 by construction; only
    FP/TN are informative here).
    """
    if lib_df is None:
        lib_df = pd.read_csv(os.path.join(data_dir, "results", "library_scores.csv"))
    ext = lib_df[(lib_df["kind"] == "external") & (lib_df["has_geometry"])].copy()
    pred_pos = (ext["predicted_label"] == MIXED_EXTERNAL_WITH_VIBRATION).to_numpy()
    gt_pos = np.zeros(len(ext), dtype=bool)  # exact completeness -> ground truth always CLEAN
    stats = _confusion_from_bools(pred_pos, gt_pos)
    stats["n_molecules"] = int(ext["molecule"].nunique())
    return ext, stats


def run_flag_validation_pipeline(data_dir="data", thresholds=None, write=True):
    """Headless entry point (Phase-6 recommended item). Writes
    data/results/benzene_EMIT_flag_confusion.csv (per-mode detail; the
    library-external check is a pure re-read of the already-committed
    library_scores.csv and is not separately persisted).

    Returns (benzene_detail_df, benzene_stats, library_detail_df, library_stats).
    """
    benzene_detail, benzene_stats = benzene_emit_flag_confusion(data_dir, thresholds)
    library_detail, library_stats = library_external_flag_confusion(data_dir)
    if write:
        out = os.path.join(data_dir, "results", "benzene_EMIT_flag_confusion.csv")
        benzene_detail.to_csv(out, index=False, float_format="%.6f")
    return benzene_detail, benzene_stats, library_detail, library_stats
