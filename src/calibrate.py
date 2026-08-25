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

tau_SB (binary-scheme S/B split, see sweep_tau_sb/run_tau_sb_error_analysis
below) is a SEPARATE, ADVISORY-ONLY analysis: it sweeps a single cutoff
against classification error (not tau_TR-style label-change stability) over
two scopes -- "all" (single-centre roster union the test tier, i.e.
everything except C6H6) and "test" (the 18-molecule held-out tier only) --
and reports where the error is minimized, WITHOUT overwriting the frozen
tau_SB default. Outputs: data/results/tau_sb_sensitivity_sweep.csv and an
advisory "tau_SB_error_sweep" block merged into thresholds.json.
"""
import json
import os

import numpy as np
import pandas as pd

from src.classifier import Thresholds, vib_label, predicted_category_column
from src.scoring import get_v_weighting
from src.library_ingest import (
    build_library_scores, resolve_log_basename, _EXTERNAL_SLOTS,
    multi_centre_molecules, out_of_calibration_scope_molecules,
    test_tier_molecules,
)

DEFAULT_TAU_GRID = tuple(round(x, 4) for x in np.arange(0.05, 0.9991, 0.005))

# Grid for the tau_SB (binary-scheme) error sweep -- coarser (0.01 step,
# 0.00-1.00) than DEFAULT_TAU_GRID since this sweep is a cheap post-hoc
# pass over the already-computed V_Stretch column (no re-scoring), not a
# full library re-classification per grid point like sweep_tau_tr.
DEFAULT_TAU_SB_GRID = tuple(round(x, 3) for x in np.arange(0.0, 1.001, 0.01))

# Calibration-scope filter: the hydride-library validation (tau_S/tau_B
# derivation + ideal/non-ideal confusion stats) is scoped to single-centre
# AB_n topologies only. As of 2026-08-11 this is expressed as an INCLUSION
# filter -- mol_type in {'ideal', 'non-ideal'} -- rather than an exclusion
# blacklist, so any mol_type outside that pair (multi-centre C6H6, and the
# 9-molecule 'test' transferability set added this session) drops out of
# calibration automatically, with no hardcoded special-case needed for a
# future category either. out_of_calibration_scope_molecules() returns that
# inclusion filter's complement (currently C6H6 + the 9 test molecules = 10
# names). C6H6 has its own separate confusion matrix in
# src/benzene_validation.py, unaffected by this filter. H2O (single-centre
# AB2, mol_type=='non-ideal') is NOT excluded.
#
# IMPORTANT: if data/mol_list_method.csv can't be read at import time (e.g.
# cwd isn't the repo root), this silently falls back to a hardcoded frozenset
# covering both the pre-roster multi-centre scope decision AND the 9 known
# test molecules -- keeps the module importable rather than crashing, but a
# future editor should know this fallback exists.
try:
    SINGLE_CENTRE_ONLY_EXCLUDE = out_of_calibration_scope_molecules()
except Exception:
    SINGLE_CENTRE_ONLY_EXCLUDE = frozenset({
        "C2H2", "C2H4", "C2H6", "H2O2", "C6H6", "iso-C4H10", "n-C4H10",
        "CH4", "C4H4", "C10H16", "PCl5", "C3H6", "B3N3H6", "CHCl3",
        "CH3CN", "C3O3H6",
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

    2026-08-25: Step 3 (T/R identification) no longer depends on tau_TR (or
    scheme) at all -- see src/classifier.py's restructuring -- so the
    Hungarian assignment itself is computed ONCE per molecule (not once per
    grid point, a genuine simplification over the pre-2026-08-25 version of
    this function), and tau_TR is applied as a purely POST-HOC filter over
    each mode's already-fixed `tr_label`/`tr_score`: a library T/R reference
    row counts "clean" at a given tau_TR iff it won its OWN slot (tr_label
    == its own name -- true by construction for a real geometry-backed
    reference, Eckart-Sayvetz completeness) AND |tr_score| >= tau_TR; a
    benzene EMIT mode counts "flagged" iff it won ANY slot (tr_label is not
    None) AND |tr_score| >= tau_TR. This reproduces the same "how stable is
    the purity call as tau_TR moves" question the pre-2026-08-25 two-gate
    version answered, now expressed as a diagnostic lens rather than a
    pipeline gate (tau_TR itself no longer changes any label).
    """
    from src.classifier import classify_all_modes

    pool = _load_geometry_pool(lib_df, data_dir)
    scorer_e, final_e = _load_benzene_emit(data_dir)

    # Match whatever weighting the pool's scorers were built with, so a
    # `--v-weighting none --calibrate` run doesn't trip its own guard.
    # scheme="threeway" pinned explicitly: this sweep determines the
    # THREE-WAY tau_S/tau_B boundaries specifically, regardless of
    # classify_all_modes()'s own default (binary as of the 2026-08 scheme
    # switch) -- flipping that default must not silently change what this
    # sweep measures. tau_TR itself is irrelevant to Step 3's assignment, so
    # any placeholder value works here; the real sweep happens below.
    th = Thresholds(tau_S=tau_S, tau_B=tau_B, v_weighting=get_v_weighting())

    tr_by_mol_slot = {}  # (mol, slot) -> (tr_label, tr_score) for T/R ground-truth rows
    for mol, (scorer, final) in pool.items():
        scored = classify_all_modes(scorer, final, th, scheme="threeway")
        for m in scored:
            if m["name"] in _EXTERNAL_SLOTS:
                tr_by_mol_slot[(mol, m["name"])] = (m["tr_label"], m["tr_score"])
    total = len(tr_by_mol_slot)

    scored_e = classify_all_modes(scorer_e, final_e, th, scheme="threeway")
    tr_emit = {m["name"]: (m["tr_label"], m["tr_score"]) for m in scored_e}

    records = []
    prev_labels = None
    for tau_TR in tau_grid:
        labels = {}
        correct = 0
        for (mol, slot), (tr_label, tr_score) in tr_by_mol_slot.items():
            is_clean = (tr_label == slot and tr_score is not None
                        and abs(tr_score) >= tau_TR)
            correct += int(is_clean)
            labels[(mol, slot)] = is_clean
        for name, (tr_label, tr_score) in tr_emit.items():
            is_flagged = tr_label is not None and abs(tr_score) >= tau_TR
            labels[("benzene_EMIT", name)] = is_flagged

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
    weighting = get_v_weighting()
    thresholds = Thresholds(tau_TR=tau_TR, tau_S=tau_S, tau_B=tau_B,
                            v_weighting=weighting)
    result = {
        "tau_TR": tau_TR, "tau_S": tau_S, "tau_B": tau_B,
        # tau_SB is NOT computed by this sweep (see Thresholds' docstring) --
        # it is carried through from the Thresholds dataclass default (the
        # sole source of truth for it) purely so thresholds.json stays a
        # complete, self-documenting record of what a run was classified
        # under. 2026-08-20: previously omitted here entirely, which meant
        # any hand-edit of thresholds.json's "tau_SB" key was silently wiped
        # by the next real --calibrate run (this dict fully overwrites the
        # file) -- see Thresholds' docstring for the full incident writeup.
        "tau_SB": thresholds.tau_SB,
        # Which eq:vscore definition these cut points were read off. Scoring
        # runs under a different weighting are refused, not silently relabelled.
        "v_weighting": weighting,
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


def confusion_matrix_stats(lib_df, thresholds, acceptance_floor=0.95, scheme="binary"):
    """Clean-category confusion matrix + per-category precision/recall
    (fig:confusion's numbers), restricted to the single-centre AB_n scope
    (SINGLE_CENTRE_ONLY_EXCLUDE, applied first). Reference labels:
    'stretch'/'bend' for internal rows, 'translation'/'rotation' for
    external rows (exact ground truth, Eckart-Sayvetz). Predicted category
    (2026-08-25 restructuring, see src/classifier.py's module docstring):
    the canonical rule `predicted_category(predicted_tr_label,
    vib_label(V_Stretch, thresholds, scheme))` -- the Step-3 winning slot if
    this mode was assigned one (`predicted_tr_label`), else a FRESH Step-2
    vib_label recomputed from `V_Stretch` under `thresholds`/`scheme` (not
    the already-stored `predicted_vib_label` column, so a caller passing
    alternate thresholds -- e.g. a sweep -- gets genuinely re-derived
    labels, not the stale ones baked into library_scores.csv).

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
    from src.classifier import classification_bucket

    # Applied first, regardless of whether the caller already pre-filtered.
    lib_df = filter_single_centre_library(lib_df)

    # classification_bucket() does the label->bucket mapping logically
    # (src/classifier.py), not via a flat dict, since external slot labels
    # are 6 distinct axis-specific strings each (e.g. "Tx".."Rz").

    df = lib_df.copy()
    df = df[df["ref_label"].notna()]

    # Exclude any ref_label outside the 4 recognized categories: a literal
    # literature "SB" (genuine mixed) label, e.g. benzene modes 21/22, would
    # otherwise inflate a category's n_pred (precision's denominator)
    # without ever contributing a true positive. "SB" answers a different
    # question (does the literature call this mode genuinely mixed?) than
    # this function's retention/precision accounting (did a stretch/bend
    # keep its label under non-ideal mass effects?) -- see
    # src/benzene_validation.py for the former's own 3-class table.
    _KNOWN_REF_LABELS = ("stretch", "bend", "translation", "rotation")
    df = df[df["ref_label"].isin(_KNOWN_REF_LABELS)]

    df = df.copy()
    fresh_vib = df["V_Stretch"].map(lambda v: vib_label(v, thresholds, scheme))
    pred_cat = predicted_category_column(df["predicted_tr_label"], fresh_vib)
    df["_pred_bucket"] = pred_cat.map(classification_bucket)

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


# --------------------------------------------------------------------------
# tau_SB (binary-scheme) error-vs-threshold sweep.
#
# Distinct from sweep_tau_tr/calibrate() above: this sweeps a single cutoff
# (tau_SB, the binary scheme's S/B split -- src/classifier.py's
# vib_label_binary) against classification ERROR relative to the literature
# reference labels, not against tau_TR label-change stability. It is purely
# advisory: it reports where the error is minimized over two scopes so a
# human can decide, by hand, whether to adopt that value as the new default
# via --tau-sb or by editing thresholds.json -- it never overwrites the
# frozen tau_SB field itself.
# --------------------------------------------------------------------------

def sweep_tau_sb(lib_df, tau_grid=DEFAULT_TAU_SB_GRID, data_dir="data"):
    """Sweep a single S/B cutoff (tau_SB) over `tau_grid`, scoring
    classification error against the literature ref_label for two molecule
    scopes: "all" (single-centre roster union the 18-molecule test tier --
    everything except C6H6, the only multi-centre molecule) and "test" (the
    18 held-out mol_type=='test' molecules only, the genuinely out-of-sample
    check).

    Restricted to internal rows with a known binary ground truth
    (ref_label in {"stretch", "bend"}) -- this drops literal "SB" reference
    rows (a different question: "is this mode genuinely mixed?", not
    answered by a binary S/B classifier) and any unlabeled row, mirroring
    confusion_matrix_stats's _KNOWN_REF_LABELS exclusion rationale.

    No re-scoring needed: V_Stretch is already threshold-independent (see
    derive_stretch_bend_thresholds's precedent), so this is a cheap post-hoc
    pass over the already-computed library_scores.csv column.

    Returns a DataFrame with columns tau_SB, error_all, accuracy_all, n_all,
    error_test, accuracy_test, n_test, error_single_centre,
    accuracy_single_centre, n_single_centre. The "single_centre" scope is
    "all" minus the test tier -- i.e. exactly filter_single_centre_library's
    scope (mol_type in {'ideal', 'non-ideal'}), the calibration-only
    population with no held-out transferability molecules mixed in.
    """
    df = lib_df[(lib_df["kind"] == "internal") & (lib_df["ref_label"].isin(("stretch", "bend")))]

    multi_centre = multi_centre_molecules(data_dir)
    test_mols = test_tier_molecules(data_dir)
    all_df = df[~df["molecule"].isin(multi_centre)]
    test_df = df[df["molecule"].isin(test_mols)]
    single_centre_df = all_df[~all_df["molecule"].isin(test_mols)]

    n_all = len(all_df)
    n_test = len(test_df)
    n_single_centre = len(single_centre_df)
    ref_all = all_df["ref_label"].to_numpy()
    v_all = all_df["V_Stretch"].to_numpy(dtype=float)
    ref_test = test_df["ref_label"].to_numpy()
    v_test = test_df["V_Stretch"].to_numpy(dtype=float)
    ref_sc = single_centre_df["ref_label"].to_numpy()
    v_sc = single_centre_df["V_Stretch"].to_numpy(dtype=float)

    records = []
    for tau in tau_grid:
        pred_all = np.where(v_all >= tau, "stretch", "bend")
        error_all = float((pred_all != ref_all).mean()) if n_all else float("nan")
        pred_test = np.where(v_test >= tau, "stretch", "bend")
        error_test = float((pred_test != ref_test).mean()) if n_test else float("nan")
        pred_sc = np.where(v_sc >= tau, "stretch", "bend")
        error_sc = float((pred_sc != ref_sc).mean()) if n_single_centre else float("nan")
        records.append({
            "tau_SB": tau,
            "error_all": error_all, "accuracy_all": 1 - error_all if n_all else float("nan"),
            "n_all": n_all,
            "error_test": error_test, "accuracy_test": 1 - error_test if n_test else float("nan"),
            "n_test": n_test,
            "error_single_centre": error_sc,
            "accuracy_single_centre": 1 - error_sc if n_single_centre else float("nan"),
            "n_single_centre": n_single_centre,
        })
    return pd.DataFrame(records)


def _error_plateau(sweep_df, error_col):
    """Longest contiguous run of grid points within 1e-9 of `error_col`'s
    minimum (mirrors find_plateau's contiguous-run logic, but keyed off
    minimum error rather than the tau_TR sweep's change-fraction+accuracy
    criterion, since this sweep has no "label-change" signal of its own).
    Returns (tau_lo, tau_hi, midpoint, min_error).
    """
    min_error = sweep_df[error_col].min()
    ok = sweep_df[error_col] <= min_error + 1e-9

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
        raise ValueError(f"No plateau found for {error_col}: this should be unreachable "
                          "since the grid's own minimum always satisfies its own tolerance.")
    lo = float(sweep_df["tau_SB"].iloc[best_start])
    hi = float(sweep_df["tau_SB"].iloc[best_start + best_len - 1])
    midpoint = round((lo + hi) / 2, 4)
    return lo, hi, midpoint, float(min_error)


def run_tau_sb_error_analysis(data_dir="data", tau_grid=DEFAULT_TAU_SB_GRID, write=True):
    """Headless entry point: ingest the library fresh, sweep tau_SB, and
    (optionally) write data/results/tau_sb_sensitivity_sweep.csv plus an
    ADVISORY "tau_SB_error_sweep" block merged into
    data/results/thresholds.json -- WITHOUT touching the frozen tau_SB field
    itself (that stays whatever Thresholds's default/JSON value already is;
    this analysis only reports where the error is minimized, it does not
    silently adopt that value).

    Returns (sweep_df, plateau_all, plateau_test, path_sweep,
    plateau_single_centre), where plateau_all/plateau_test/
    plateau_single_centre are (tau_lo, tau_hi, midpoint, min_error) tuples
    from _error_plateau.
    """
    lib_df = build_library_scores(data_dir)
    sweep_df = sweep_tau_sb(lib_df, tau_grid, data_dir)
    plateau_all = _error_plateau(sweep_df, "error_all")
    plateau_test = _error_plateau(sweep_df, "error_test")
    plateau_single_centre = _error_plateau(sweep_df, "error_single_centre")

    path_sweep = os.path.join(data_dir, "results", "tau_sb_sensitivity_sweep.csv")
    path_json = os.path.join(data_dir, "results", "thresholds.json")
    if write:
        sweep_df.to_csv(path_sweep, index=False)

        result = {}
        if os.path.exists(path_json):
            with open(path_json) as f:
                result = json.load(f)
        lo_all, hi_all, mid_all, err_all = plateau_all
        lo_test, hi_test, mid_test, err_test = plateau_test
        lo_sc, hi_sc, mid_sc, err_sc = plateau_single_centre
        result["tau_SB_error_sweep"] = {
            "grid": "0.00-1.00 step 0.01",
            "optimal_tau_all": mid_all, "min_error_all": err_all,
            "n_all": int(sweep_df["n_all"].iloc[0]) if len(sweep_df) else 0,
            "optimal_tau_test": mid_test, "min_error_test": err_test,
            "n_test": int(sweep_df["n_test"].iloc[0]) if len(sweep_df) else 0,
            "optimal_tau_single_centre": mid_sc, "min_error_single_centre": err_sc,
            "n_single_centre": int(sweep_df["n_single_centre"].iloc[0]) if len(sweep_df) else 0,
        }
        with open(path_json, "w") as f:
            json.dump(result, f, indent=2)
    return sweep_df, plateau_all, plateau_test, path_sweep, plateau_single_centre
