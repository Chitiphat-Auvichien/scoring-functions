"""Phase-3 threshold calibration: derive tau_S/tau_B from the library's
literature stretch/bend labels, sweep tau_TR for a stability plateau, and
freeze the final Thresholds to data/results/thresholds.json.

Design
------
tau_S / tau_B (Step 4, internal stretch/bend split)
    Derived from the IDEAL-molecule subset of the ingested library
    (``library_scores.csv``'s ``ideal == 'yes'`` rows) ONLY, not the full
    ideal+non-ideal set. Rationale: group theory guarantees that in an ideal
    (high-symmetry) molecule every mode's irreducible representation is
    either uniquely stretching or uniquely bending, so the ideal-molecule
    s[V_S] values realize the framework's true, noise-free step function.
    Verified empirically this session: the ideal-labeled stretch population's
    minimum is 0.90368 and the ideal-labeled bend population's maximum is
    0.17327 -- a clean, NON-OVERLAPPING gap of width ~0.73 -- so a threshold
    pair drawn from these two boundary values separates the calibration
    population with zero error, by construction (verified by an assertion,
    not assumed). Calibrating on the ideal-only population and then
    evaluating against the FULL library (ideal+non-ideal) is the honest
    procedure: non-ideal modes are exactly the population the manuscript
    predicts will migrate into the mixed bucket (center-of-mass-driven
    softening, B8.3), so folding them into the calibration set would
    contaminate the threshold with the very effect it exists to detect.
    tau_S/tau_B are taken as the EXACT empirical boundary values (no ad hoc
    rounding to e.g. the nearest 0.05) -- this is the tightest boundary that
    still achieves exact separation on the calibration population, and is
    fully reproducible from the data with no subjective rounding choice.

tau_TR (Step 3, external purity gate 1)
    Swept over a grid against TWO evaluation sets, evaluated once per grid
    point (Step 2's Hungarian assignment does not depend on tau_TR, only
    Step 3's clean/mixed decision does, but classify_all_modes is re-run in
    full each time for simplicity -- these are small molecules, so the cost
    is negligible):
      1. every geometry-backed library molecule's REAL normal modes (water,
         CO2, benzene, + the hydride-library molecules that ship a .log/.gjf
         pair -- 25 molecules total). Ground truth: because normal-mode
         translation/rotation references are constructed directly from
         geometry (Eckart-Sayvetz), completeness is EXACT for them -- every
         one of the n_T+n_R ideal references in this set MUST classify clean
         at any sane tau_TR. This gives a genuine accuracy metric.
      2. benzene's 36 EMIT modes -- the paper's sole strongly-mixed-mode
         stress test with real, already-characterized behavior (EMIT 2/9's
         Ry-inversion, EMIT 34/35/36's flag/blind-spot triad). No independent
         ground truth exists for all 36 EMIT modes (that full-set validation
         is Phase 6, out of scope here), so EMIT contributes to the
         label-change-fraction signal only, not to the accuracy metric.
    Plateau criterion (stated explicitly, per the task): scanning tau_TR on a
    0.005 grid from 0.05 to 0.999, the plateau is the LONGEST CONTIGUOUS run
    of grid points for which (a) the label-change fraction relative to the
    immediately preceding grid point is EXACTLY ZERO across the combined
    evaluation set (25 molecules' external slots + benzene's 36 EMIT modes),
    AND (b) accuracy against the normal-mode ground truth is at its run
    maximum. A step size of 0.005 is fine enough that any genuine tau_TR
    sensitivity shows up as a nonzero change fraction at that resolution, so
    a strictly-zero run is a meaningful "nothing changes here" band, not an
    artifact of coarse sampling. tau_TR is then frozen at 0.95 (the value
    used as a worked example in the algorithm spec, PDF Section B6.2) if 0.95
    falls inside the plateau, else at the plateau's midpoint (documented
    either way in the output JSON's ``plateau_tau_TR_range``).

Outputs
-------
data/results/thresholds.json       -- frozen {tau_TR, tau_S, tau_B} + the
                                       derivation stats + plateau range.
data/results/tau_sensitivity_sweep.csv -- the full per-grid-point
                                       (tau_TR, accuracy, label_change_fraction)
                                       curve, for fig:sensitivity.
"""
import json
import os

import numpy as np
import pandas as pd

from src.classifier import Thresholds, is_clean_external
from src.excel_ingest import (
    build_library_scores, resolve_log_basename, _EXTERNAL_SLOTS,
)

DEFAULT_TAU_GRID = tuple(round(x, 4) for x in np.arange(0.05, 0.9991, 0.005))


def derive_stretch_bend_thresholds(lib_df):
    """tau_S / tau_B from the ideal-molecule subset's stretch/bend V_Stretch
    distributions (see module docstring). Returns (tau_S, tau_B, stats).
    """
    internal = lib_df[lib_df["kind"] == "internal"]
    ideal = internal[internal["ideal"] == "yes"]
    stretch = ideal.loc[ideal["ref_label"] == "stretch", "V_Stretch"].dropna()
    bend = ideal.loc[ideal["ref_label"] == "bend", "V_Stretch"].dropna()
    if len(stretch) == 0 or len(bend) == 0:
        raise ValueError("Ideal-molecule stretch or bend population is empty; "
                          "cannot calibrate tau_S/tau_B.")

    tau_S = float(stretch.min())
    tau_B = float(bend.max())
    if not (tau_B < tau_S):
        raise ValueError(
            f"Ideal stretch/bend V_Stretch populations overlap (bend max "
            f"{tau_B} >= stretch min {tau_S}) -- the zero-overlap assumption "
            "calibration relies on does not hold for this library; revisit "
            "the threshold-derivation method rather than silently rounding "
            "past the overlap.")

    stats = {
        "ideal_stretch_n": int(len(stretch)), "ideal_stretch_min": tau_S,
        "ideal_stretch_mean": float(stretch.mean()),
        "ideal_bend_n": int(len(bend)), "ideal_bend_max": tau_B,
        "ideal_bend_mean": float(bend.mean()),
        "gap_width": tau_S - tau_B,
    }
    return tau_S, tau_B, stats


def _load_geometry_pool(lib_df, data_dir="data"):
    """(scorer, final) pairs for every geometry-backed library molecule, built
    once and reused across the whole tau_TR grid."""
    from main import load_inputs, build_scorer_and_final

    pool = {}
    molecules = sorted(lib_df.loc[lib_df["has_geometry"], "molecule"].unique())
    for mol in molecules:
        base = resolve_log_basename(mol, data_dir)
        if base is None:
            continue
        raw, _ = load_inputs(base, "normal", data_dir)
        pool[mol] = build_scorer_and_final(raw, "normal")
    return pool


def _load_benzene_emit(data_dir="data"):
    from main import load_inputs, build_scorer_and_final
    raw, _ = load_inputs("benzene", "emit", data_dir)
    return build_scorer_and_final(raw, "emit")


def sweep_tau_tr(lib_df, tau_S, tau_B, data_dir="data", tau_grid=DEFAULT_TAU_GRID):
    """Sweep tau_TR over `tau_grid`; return a DataFrame with one row per grid
    point: tau_TR, accuracy (library normal-mode T/R ground truth), and
    label_change_fraction (combined library-external + benzene-EMIT set,
    relative to the previous grid point; NaN for the first point).
    """
    from src.classifier import classify_all_modes

    pool = _load_geometry_pool(lib_df, data_dir)
    scorer_e, final_e = _load_benzene_emit(data_dir)

    records = []
    prev_labels = None
    for tau_TR in tau_grid:
        th = Thresholds(tau_TR=tau_TR, tau_S=tau_S, tau_B=tau_B)
        labels = {}
        correct = 0
        total = 0
        for mol, (scorer, final) in pool.items():
            scored = classify_all_modes(scorer, final, th)
            for m in scored:
                if m["name"] in _EXTERNAL_SLOTS:
                    total += 1
                    is_clean = is_clean_external(m["classification"])
                    correct += int(is_clean)
                    labels[(mol, m["name"])] = m["classification"]
        scored_e = classify_all_modes(scorer_e, final_e, th)
        for m in scored_e:
            labels[("benzene_EMIT", m["name"])] = m["classification"]

        accuracy = correct / total if total else float("nan")
        if prev_labels is None:
            change_frac = 0.0
        else:
            n_changed = sum(1 for k, v in labels.items() if prev_labels.get(k) != v)
            change_frac = n_changed / len(labels)
        records.append({
            "tau_TR": tau_TR, "accuracy": accuracy,
            "label_change_fraction": change_frac, "n_external_ground_truth": total,
        })
        prev_labels = labels
    return pd.DataFrame(records)


def find_plateau(sweep_df, change_tol=0.0):
    """Longest contiguous run of grid points with label_change_fraction <=
    change_tol AND accuracy at the run's maximum. Returns (tau_lo, tau_hi).
    """
    max_acc = sweep_df["accuracy"].max()
    ok = (sweep_df["label_change_fraction"] <= change_tol) & \
         (sweep_df["accuracy"] >= max_acc - 1e-9)

    best_start = best_len = cur_start = cur_len = 0
    for i, v in enumerate(ok.tolist()):
        if v:
            if cur_len == 0:
                cur_start = i
            cur_len += 1
            if cur_len > best_len:
                best_len, best_start = cur_len, cur_start
        else:
            cur_len = 0
    if best_len == 0:
        raise ValueError("No plateau found: label_change_fraction is never 0 "
                          "at the grid's accuracy maximum -- the tau_TR grid "
                          "or evaluation set needs revisiting.")
    lo = float(sweep_df["tau_TR"].iloc[best_start])
    hi = float(sweep_df["tau_TR"].iloc[best_start + best_len - 1])
    return lo, hi


def freeze_tau_tr(sweep_df, preferred=0.95):
    lo, hi = find_plateau(sweep_df)
    if lo <= preferred <= hi:
        return preferred, (lo, hi)
    return round((lo + hi) / 2, 4), (lo, hi)


def calibrate(lib_df, data_dir="data", tau_grid=DEFAULT_TAU_GRID, preferred_tau_tr=0.95):
    """Full Phase-3 calibration. Returns (Thresholds, result_dict, sweep_df)."""
    tau_S, tau_B, sb_stats = derive_stretch_bend_thresholds(lib_df)
    sweep_df = sweep_tau_tr(lib_df, tau_S, tau_B, data_dir, tau_grid)
    tau_TR, plateau = freeze_tau_tr(sweep_df, preferred_tau_tr)

    thresholds = Thresholds(tau_TR=tau_TR, tau_S=tau_S, tau_B=tau_B)
    result = {
        "tau_TR": tau_TR, "tau_S": tau_S, "tau_B": tau_B,
        "plateau_tau_TR_range": list(plateau),
        "plateau_criterion": (
            "Longest contiguous run of the tau_TR grid (step 0.005, 0.05-0.999) "
            "with zero label changes vs. the previous grid point AND accuracy "
            "at its run maximum, evaluated over {25 geometry-backed library "
            "molecules' real normal-mode T/R references [ground truth: always "
            "clean, exact completeness] + benzene's 36 EMIT modes [no "
            "independent ground truth; label-change signal only]}."
        ),
        "stretch_bend_derivation": sb_stats,
    }
    return thresholds, result, sweep_df


def confusion_matrix_stats(lib_df, thresholds, acceptance_floor=0.95):
    """Clean-category confusion matrix + per-category precision/recall
    (fig:confusion's underlying numbers), evaluated over the WHOLE ingested
    library (both Excel-only and geometry-backed rows).

    Reference labels: 'stretch'/'bend' for internal rows, 'translation'/
    'rotation' for external rows (the latter are ground truth by
    construction -- Eckart-Sayvetz normal-mode T/R references are exact).

    Predicted labels:
      - internal rows WITHOUT geometry: vib_label(V_Stretch, thresholds)
        applied directly (Step 4 only). This is not a simplification for
        lack of data -- data_score's internal rows are already Gaussian's
        own T/R-projected-out vibrational modes, so there is no external
        character for Step 2/3 to catch even if geometry were available;
        see excel_ingest.py's module docstring for the full argument.
      - internal/external rows WITH geometry: the already-computed
        'predicted_label' column from attach_geometry_classification()
        (full Algorithm 1, Steps 2-4) -- not recomputed here.

    Returns {'acceptance_floor', 'per_category': {...}, 'floor_met',
    'confusion_table': DataFrame} -- the last is the raw reference-label x
    predicted-bucket contingency table for fig:confusion.

    ADDITIVE (2026-07-02, formula-auditor + lead-author recommendation):
    each `per_category[cat]` entry also carries `recall_ideal`/
    `recall_nonideal` (+ `n_ref_ideal`/`n_ref_nonideal`), a formal, tested
    version of the ideal-vs-non-ideal ground-truth-strength split that
    `src/figures.py::plot_confusion_matrix` already applies ad hoc (as of
    commit 424a666) directly to `library_scores.csv` for `fig:confusion`'s
    two-tier layout. The pooled `precision`/`recall` keys above are UNCHANGED
    (existing callers/tests are unaffected) -- this only adds new keys, and
    intentionally does not touch `src/figures.py` (that figure is already
    correct and wired into the manuscript; this just gives the same split a
    single formal, computed home instead of only living inside plotting code).
    Tier masks are copied verbatim from `plot_confusion_matrix` for identical
    semantics: `ideal` tier = every external (T/R) row REGARDLESS of its
    `ideal` tag (Eckart-Sayvetz completeness makes those exact regardless)
    OR any internal row with `ideal=='yes'`; `nonideal` tier = internal rows
    with `ideal=='no'` only. For `translation`/`rotation`, every row is
    `kind=='external'`, so `recall_ideal` reproduces the pooled `recall`
    exactly and `recall_nonideal` is NaN (n_ref_nonideal=0, no external row
    ever has `ideal=='no'`) -- expected, not a bug. For `stretch`/`bend`,
    `recall_ideal` is guaranteed to be exactly 1.0: tau_S/tau_B
    (`derive_stretch_bend_thresholds`) are LITERALLY the min/max of this same
    `ideal=='yes'` population, so by construction no ideal-tier stretch/bend
    row can land on the wrong side of its own defining boundary.
    """
    from src.classifier import vib_label, classification_bucket

    # Bucket lookup is logic-based (classification_bucket(), src/classifier.py),
    # not a flat dict keyed by exact label: clean/mixed-external labels are now
    # 6 distinct axis-specific strings each (e.g. "Tx".."Rz", "Tx*".."Rz*")
    # rather than the 2 fixed CLEAN_TRANSLATION/CLEAN_ROTATION/
    # MIXED_EXTERNAL_WITH_VIBRATION constants a flat dict used to key on.

    df = lib_df.copy()
    df = df[df["ref_label"].notna()]  # every row here has a stretch/bend/
                                       # translation/rotation reference label

    # EXCLUDE any reference label outside the 4 recognized categories
    # (2026-07-05 fix, flagged after commit 69d549e gave benzene modes 21/22
    # a genuine literal literature "SB" (mixed) ground-truth label for the
    # first time). This function's confusion table and per-category
    # precision/recall/retention accounting is specifically a 4-category
    # (translation/rotation/stretch/bend) contingency check. Left in `df`, a
    # foreign ref_label whose predicted bucket happens to land on one of
    # those 4 (both mode 21 and 22 predict bucket "bend") would silently
    # inflate that category's n_pred -- precision's denominator -- without
    # ever being able to contribute a true positive, corrupting bend
    # precision from 1.000 to 0.99267 for a reason that has nothing to do
    # with classifier error. A literal literature "SB" row answers a
    # different question ("does the literature call this mode genuinely
    # mixed?" -- see src/benzene_validation.py::benzene_internal_confusion_matrix
    # for that benzene-scoped 3-class table) than "did a nominal
    # stretch/bend keep its label or migrate to mixed under non-ideal mass
    # effects?", which is what this function's retention/precision
    # accounting measures. Any other currently-unrecognized ref_label would
    # hit the same silent miscount, so this is a general guard (not a
    # benzene-only special case) and does not redesign the function or add a
    # 3rd category to its own bucket vocabulary.
    _KNOWN_REF_LABELS = ("stretch", "bend", "translation", "rotation")
    df = df[df["ref_label"].isin(_KNOWN_REF_LABELS)]

    def _predict(row):
        if row["kind"] == "internal" and not row["has_geometry"]:
            return classification_bucket(vib_label(row["V_Stretch"], thresholds))
        return classification_bucket(row["predicted_label"])

    df = df.copy()
    df["_pred_bucket"] = df.apply(_predict, axis=1)

    confusion_table = pd.crosstab(df["ref_label"], df["_pred_bucket"])

    # Ideal/non-ideal ground-truth-strength tiers -- identical masks to
    # src/figures.py::plot_confusion_matrix's ad hoc split (see docstring
    # above); computed once here and reused for every category below.
    ideal_tier_mask = (df["kind"] == "external") | (df["ideal"] == "yes")
    nonideal_tier_mask = (df["kind"] == "internal") & (df["ideal"] == "no")

    per_category = {}
    for ref in ("stretch", "bend", "translation", "rotation"):
        ref_mask = df["ref_label"] == ref
        pred_mask = df["_pred_bucket"] == ref
        n_ref = int(ref_mask.sum())
        n_pred = int(pred_mask.sum())
        tp = int((ref_mask & pred_mask).sum())
        recall = tp / n_ref if n_ref else float("nan")
        precision = tp / n_pred if n_pred else float("nan")
        entry = {"tp": tp, "n_ref": n_ref, "n_pred": n_pred,
                 "precision": precision, "recall": recall}
        if ref in ("stretch", "bend"):
            entry["mixed_fraction"] = float((df.loc[ref_mask, "_pred_bucket"] == "mixed").mean())

        for tier_name, tier_mask in (("ideal", ideal_tier_mask), ("nonideal", nonideal_tier_mask)):
            tier_ref_mask = ref_mask & tier_mask
            n_ref_tier = int(tier_ref_mask.sum())
            tp_tier = int((tier_ref_mask & pred_mask).sum())
            entry[f"n_ref_{tier_name}"] = n_ref_tier
            entry[f"recall_{tier_name}"] = tp_tier / n_ref_tier if n_ref_tier else float("nan")

        per_category[ref] = entry

    floor_met = all(
        per_category[c]["precision"] >= acceptance_floor and
        per_category[c]["recall"] >= acceptance_floor
        for c in ("stretch", "bend", "translation", "rotation")
    )
    return {
        "acceptance_floor": acceptance_floor,
        "per_category": per_category,
        "floor_met": floor_met,
        "confusion_table": confusion_table,
    }


def run_calibration_pipeline(data_dir="data",
                              xlsx_path=None,
                              tau_grid=DEFAULT_TAU_GRID,
                              preferred_tau_tr=0.95,
                              write=True,
                              source="excel"):
    """Headless entry point: ingest the library fresh, calibrate, and
    (optionally) write data/results/thresholds.json +
    data/results/tau_sensitivity_sweep.csv. Returns
    (thresholds, result_dict, sweep_df, (path_json, path_sweep)).

    `source` ("excel" default, or "gaussian") is passed straight through to
    src.excel_ingest.build_library_scores() -- see its docstring for the
    2026-07-04 dual-source contract. Default matches run_ingest_pipeline()'s
    default so --library and --calibrate calibrate against the same
    population unless told otherwise.
    """
    xlsx_path = xlsx_path or os.path.join(data_dir, "vibrational-scoring-functions.xlsx")
    lib_df = build_library_scores(xlsx_path, data_dir, source=source)
    thresholds, result, sweep_df = calibrate(lib_df, data_dir, tau_grid, preferred_tau_tr)

    path_json = os.path.join(data_dir, "results", "thresholds.json")
    path_sweep = os.path.join(data_dir, "results", "tau_sensitivity_sweep.csv")
    if write:
        with open(path_json, "w") as f:
            json.dump(result, f, indent=2)
        sweep_df.to_csv(path_sweep, index=False)
    return thresholds, result, sweep_df, (path_json, path_sweep)
