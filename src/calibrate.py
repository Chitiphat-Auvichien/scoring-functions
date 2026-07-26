"""Threshold calibration: derive tau_S/tau_B from the library's literature
stretch/bend labels, sweep tau_TR for a stability plateau, and freeze the
final Thresholds to data/results/thresholds.json.

tau_S / tau_B (Step 4 stretch/bend split): derived from the IDEAL-molecule
subset only (``ideal == 'yes'``), since group theory guarantees an ideal
molecule's modes are unambiguously stretch or bend. Verified empirical gap:
ideal-stretch min = 0.90368, ideal-bend max = 0.17327 (width ~0.73, zero
overlap) -- tau_S/tau_B are these exact boundary values, no rounding.
Non-ideal modes are excluded from calibration because they are exactly the
population expected to migrate into the mixed bucket; including them would
contaminate the threshold with the effect it's meant to detect.

tau_TR (Step 3 purity gate): swept over a 0.005 grid (0.05-0.999) against
(1) every geometry-backed library molecule's real normal-mode T/R references
(exact ground truth by Eckart-Sayvetz completeness -- always clean) and
(2) benzene's 36 EMIT modes (label-change signal only, no independent ground
truth). The plateau is the longest contiguous run where the label-change
fraction vs. the previous grid point is exactly zero AND accuracy is at its
run maximum. tau_TR is frozen at 0.95 if it falls in the plateau, else at
the plateau's midpoint.

Outputs: data/results/thresholds.json (frozen thresholds + derivation stats
+ plateau range), data/results/tau_sensitivity_sweep.csv (full sweep curve,
fig:sensitivity).
"""
import json
import os

import numpy as np
import pandas as pd

from src.classifier import Thresholds, is_clean_external
from src.library_ingest import (
    build_library_scores, resolve_log_basename, _EXTERNAL_SLOTS,
    multi_centre_molecules,
)

DEFAULT_TAU_GRID = tuple(round(x, 4) for x in np.arange(0.05, 0.9991, 0.005))

# Single-centre-only scope filter: the hydride-library validation (tau_S/
# tau_B derivation + ideal/non-ideal confusion stats) is scoped to
# single-centre AB_n topologies only. Excludes multi-/two-centre molecules
# (currently just C6H6, benzene -- it has its own separate confusion matrix
# in src/benzene_validation.py, unaffected by this filter). H2O (single-
# centre AB2) is NOT excluded. Sourced from
# src.library_ingest.multi_centre_molecules() (mol_list_method.csv's
# 'mol_type' column) so it stays roster-driven, not hardcoded.
#
# IMPORTANT: if data/mol_list_method.csv can't be read at import time (e.g.
# cwd isn't the repo root), this silently falls back to a hardcoded 7-name
# frozenset from the pre-roster scope decision -- keeps the module importable
# rather than crashing, but a future editor should know this fallback exists.
try:
    SINGLE_CENTRE_ONLY_EXCLUDE = multi_centre_molecules()
except Exception:
    SINGLE_CENTRE_ONLY_EXCLUDE = frozenset({
        "C2H2", "C2H4", "C2H6", "H2O2", "C6H6", "iso-C4H10", "n-C4H10",
    })


def filter_single_centre_library(lib_df):
    """Drop rows belonging to molecules outside the single-centre AB_n scope
    (SINGLE_CENTRE_ONLY_EXCLUDE). Idempotent."""
    return lib_df[~lib_df["molecule"].isin(SINGLE_CENTRE_ONLY_EXCLUDE)].copy()


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


def _geometry_pool_molecules(lib_df, data_dir="data"):
    """Molecule names with has_geometry True and a roster-resolvable
    basename -- the population _load_geometry_pool() will load, without
    actually parsing anything."""
    molecules = sorted(lib_df.loc[lib_df["has_geometry"], "molecule"].unique())
    return [mol for mol in molecules if resolve_log_basename(mol, data_dir) is not None]


def _load_geometry_pool(lib_df, data_dir="data"):
    """(scorer, final) pairs for every geometry-backed library molecule, built
    once and reused across the whole tau_TR grid."""
    from main import load_inputs, build_scorer_and_final

    pool = {}
    for mol in _geometry_pool_molecules(lib_df, data_dir):
        base = resolve_log_basename(mol, data_dir)
        raw, _ = load_inputs(base, "normal", data_dir)
        pool[mol] = build_scorer_and_final(raw, "normal")
    return pool


def _load_benzene_emit(data_dir="data"):
    """Load benzene's EMIT modes via its roster-resolved basename (C6H6);
    data/EMIT/*_EMIT.txt is keyed by the same basename as .log/.gjf."""
    from main import load_inputs, build_scorer_and_final
    base = resolve_log_basename("C6H6", data_dir)
    raw, _ = load_inputs(base, "emit", data_dir)
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

    n_pool = len(_geometry_pool_molecules(lib_df, data_dir))
    thresholds = Thresholds(tau_TR=tau_TR, tau_S=tau_S, tau_B=tau_B)
    result = {
        "tau_TR": tau_TR, "tau_S": tau_S, "tau_B": tau_B,
        "plateau_tau_TR_range": list(plateau),
        "plateau_criterion": (
            "Longest contiguous run of the tau_TR grid (step 0.005, 0.05-0.999) "
            "with zero label changes vs. the previous grid point AND accuracy "
            f"at its run maximum, evaluated over {{{n_pool} geometry-backed library "
            "molecules' real normal-mode T/R references [ground truth: always "
            "clean, exact completeness] + benzene's 36 EMIT modes [no "
            "independent ground truth; label-change signal only]}."
        ),
        "stretch_bend_derivation": sb_stats,
    }
    return thresholds, result, sweep_df


def confusion_matrix_stats(lib_df, thresholds, acceptance_floor=0.95):
    """Clean-category confusion matrix + per-category precision/recall
    (fig:confusion's numbers), restricted to the single-centre AB_n scope
    (SINGLE_CENTRE_ONLY_EXCLUDE, applied first). Reference labels:
    'stretch'/'bend' for internal rows, 'translation'/'rotation' for
    external rows (exact ground truth, Eckart-Sayvetz). Predicted labels
    come from the already-computed 'predicted_label' column
    (score_geometry_molecule(), full Algorithm 1).

    Returns {'acceptance_floor', 'per_category': {...}, 'floor_met',
    'confusion_table': DataFrame} -- the contingency table for fig:confusion.

    Each `per_category[cat]` entry also carries `recall_ideal`/
    `recall_nonideal` (+ `n_ref_ideal`/`n_ref_nonideal`), an ideal-vs-non-
    ideal ground-truth-strength split: `ideal` tier = every external row
    (always exact) OR internal rows with `ideal=='yes'`; `nonideal` tier =
    internal rows with `ideal=='no'`. NON-OBVIOUS INVARIANT: for
    `stretch`/`bend`, `recall_ideal` is guaranteed exactly 1.0 by
    construction -- tau_S/tau_B are literally the min/max of this same
    `ideal=='yes'` population, so no ideal-tier row can land on the wrong
    side of its own defining boundary.
    """
    from src.classifier import vib_label, classification_bucket

    # Applied first, regardless of whether the caller already pre-filtered.
    lib_df = filter_single_centre_library(lib_df)

    # classification_bucket() does the label->bucket mapping logically
    # (src/classifier.py), not via a flat dict, since clean/mixed-external
    # labels are 6 distinct axis-specific strings each (e.g. "Tx".."Rz*").

    df = lib_df.copy()
    df = df[df["ref_label"].notna()]

    # Exclude any ref_label outside the 4 recognized categories: a literal
    # literature "SB" (genuine mixed) label, e.g. benzene modes 21/22, would
    # otherwise inflate a category's n_pred (precision's denominator)
    # without ever contributing a true positive -- this previously corrupted
    # bend precision from 1.000 to 0.99267 (commit 69d549e) before this
    # filter was added. "SB" answers a different question (does the
    # literature call this mode genuinely mixed?) than this function's
    # retention/precision accounting (did a stretch/bend keep its label
    # under non-ideal mass effects?) -- see src/benzene_validation.py for
    # the former's own 3-class table.
    _KNOWN_REF_LABELS = ("stretch", "bend", "translation", "rotation")
    df = df[df["ref_label"].isin(_KNOWN_REF_LABELS)]

    def _predict(row):
        if row["kind"] == "internal" and not row["has_geometry"]:
            return classification_bucket(vib_label(row["V_Stretch"], thresholds))
        return classification_bucket(row["predicted_label"])

    df = df.copy()
    df["_pred_bucket"] = df.apply(_predict, axis=1)

    confusion_table = pd.crosstab(df["ref_label"], df["_pred_bucket"])

    # Ideal/non-ideal ground-truth-strength tiers, see docstring above.
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
                              tau_grid=DEFAULT_TAU_GRID,
                              preferred_tau_tr=0.95,
                              write=True):
    """Headless entry point: ingest the library fresh, calibrate, and
    (optionally) write data/results/thresholds.json +
    data/results/tau_sensitivity_sweep.csv. Returns
    (thresholds, result_dict, sweep_df, (path_json, path_sweep)).

    Ingests via src.library_ingest.build_library_scores() -- the single,
    roster-driven (data/mol_list_method.csv) pipeline; matches
    run_ingest_pipeline()'s own population so --library and --calibrate
    calibrate against the same data.
    """
    lib_df = build_library_scores(data_dir)
    thresholds, result, sweep_df = calibrate(lib_df, data_dir, tau_grid, preferred_tau_tr)

    path_json = os.path.join(data_dir, "results", "thresholds.json")
    path_sweep = os.path.join(data_dir, "results", "tau_sensitivity_sweep.csv")
    if write:
        with open(path_json, "w") as f:
            json.dump(result, f, indent=2)
        sweep_df.to_csv(path_sweep, index=False)
    return thresholds, result, sweep_df, (path_json, path_sweep)
