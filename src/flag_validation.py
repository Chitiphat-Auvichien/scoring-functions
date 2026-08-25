"""Systematic flag precision/recall validation over all 36 benzene EMIT modes
(and library externals) against the projection reference (src/projection.py).

Compares classify_all_modes()'s mixed-external flag (Step 3's two-gate purity
test; an axis label with trailing "*", e.g. "Tx*") against a projection-based
ground truth: per EMIT mode, M_ext = max(C2_Tx..C2_Rz). Ground truth is MIXED
iff GT_EXT_LO < M_ext < GT_EXT_HI, else CLEAN. GT_EXT_HI=0.95 reuses tau_TR
(the classifier's own purity bar); GT_EXT_LO=0.05 is a small symmetric floor
above the ~1e-4 numerical orthonormality noise floor (projection.py).

Finding (pre-2026-08-25, PRESERVED HISTORICAL RESULT, see restructuring note
below): precision = 1.0 (flag never fires on a genuinely clean mode), recall
= 5/17 = 0.294 (misses most genuinely mixed modes). Root cause: Step 2's
one-to-one linear_sum_assignment only ever assigns 6 of 36 modes an external
slot at all (by construction, not a bug) -- the other 30 modes reach Step 4
and can never be flagged regardless of true external character. E.g. EMIT 2
has more genuine Ry character than EMIT 9 by projection (38.7% vs 14.1%) but
EMIT 9 wins the Ry slot on raw score (0.215 vs 0.143), so EMIT 9 is flagged
(TP) while EMIT 2 is missed (FN) -- the recall gap is systematic to Step 2's
assignment mechanism, not an isolated edge case.

RESTRUCTURING NOTE (2026-08-25): src/classifier.py's Step 3 (T/R
identification) no longer has any purity gate at all -- there is no more
"mixed-external flag" for this module to validate. `is_mixed_external` on a
freshly-produced `tr_label` is now definitionally always False, so
`benzene_emit_flag_confusion`'s `predicted_positive` column is 0 everywhere
by construction (TP=FP=0, recall collapses to 0 rather than the historical
0.294) -- this module's original question ("does the flag fire correctly?")
no longer has a live subject. Kept functional (not deleted) purely so old
call sites/tests don't crash, and because `ground_truth_label`'s projection-
based MIXED/CLEAN classification is still a valid, independent diagnostic on
its own -- but its precision/recall numbers against `predicted_positive` are
no longer a meaningful pipeline-behavior claim. This was already excluded
from the manuscript (2026-07-09, see IMPLEMENTATION_PLAN.md) as an internal
diagnostic only; not revisited further here (out of scope for the Step-3
restructuring pass -- flag this module for removal or a real rewrite in a
future session if it is still wanted).
"""
import os

import numpy as np
import pandas as pd

from src.classifier import Thresholds, is_mixed_external, predicted_category_column
from src.library_ingest import load_library_scores

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
    from src.library_ingest import resolve_log_basename

    thresholds = thresholds or Thresholds.calibrated()

    # Resolve benzene's roster name (C6H6) to its on-disk basename.
    base = resolve_log_basename("C6H6", data_dir)
    raw, _ = load_inputs(base, "emit", data_dir)
    scorer, final = build_scorer_and_final(raw, "emit")
    scored = classify_all_modes(scorer, final, thresholds)
    by_name = {m["name"]: m for m in scored}

    df_contrib, _, _ = run_projection_pipeline(base, data_dir, thresholds, write=False)
    df_contrib = df_contrib.copy()
    df_contrib["M_ext"] = df_contrib[list(_EXTERNAL_C2_COLS)].max(axis=1)

    rows = []
    for _, r in df_contrib.iterrows():
        name = r["Mode"]
        m = by_name[name]
        # 2026-08-25: no combined "classification" field anymore -- and
        # since Step 3 no longer gates, is_mixed_external() on a
        # freshly-produced tr_label is always False (see module docstring's
        # restructuring note). Reconstructed only for the CSV's own
        # informational "classifier_label" column.
        classifier_label = m["tr_label"] or m["vib_label"]
        m_ext = float(r["M_ext"])
        gt_pos = ext_lo < m_ext < ext_hi
        pred_pos = is_mixed_external(classifier_label)
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
            "classifier_label": classifier_label,
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
    """Clean-vs-mixed-external check on the library molecules' real normal-
    mode T/R references. Ground truth is degenerate by construction: a real
    T/R reference is one of Q's own columns, orthonormal to the rest by
    Eckart-Sayvetz completeness, so M_ext=1.0 and ground truth is CLEAN for
    every row -- this verifies (not assumes) the classifier's flag never
    fires here (FP=0).

    Returns (ext_df, stats): ext_df is the (kind=='external' & has_geometry)
    subset of library_scores.csv; stats is _confusion_from_bools() output
    (TP=FN=0 by construction; only FP/TN are informative).
    """
    lib_df = load_library_scores(data_dir, lib_df)
    ext = lib_df[(lib_df["kind"] == "external") & (lib_df["has_geometry"])].copy()
    predicted_label = predicted_category_column(
        ext["predicted_tr_label"], ext["predicted_vib_label"])
    pred_pos = predicted_label.apply(is_mixed_external).to_numpy()
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
