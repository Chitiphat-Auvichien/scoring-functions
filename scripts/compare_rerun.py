"""Rerun consistency check: does this project's classification pipeline
reproduce the same T/R/V scores and predicted labels when the raw
Gaussian calculation is redone on different software/hardware?

Commit 33d9efe added data/logs_rerun/ + data/gjf_rerun/ + data/fchk_rerun/ --
the same 68-molecule roster as data/mol_list_method.csv, same per-molecule
method/basis (mol_list_method.csv's `current_method` column), but the raw
Gaussian calculations were redone on Gaussian 16 on a different machine,
instead of the original Gaussian 09 run in data/logs/ + data/gjf/.
data/fchk_rerun/ is not read here -- nothing in the scoring/classification
pipeline consumes .fchk files.

Usage (from Github/scoring-functions/):
    py scripts/compare_rerun.py

What this does
--------------
1. (Re)builds a throwaway mirror data directory, data/_rerun_mirror/
   (gitignored, idempotent, safe to delete/re-run):
     - logs/, gjf/  -- NTFS directory junctions (mklink /J, no admin rights
       needed) into data/logs_rerun / data/gjf_rerun.
     - EMIT/        -- junctioned into the CANONICAL data/EMIT. EMIT modes
       were NOT part of this rerun batch; junctioning the canonical,
       unaffected EMIT directory just lets code paths that need raw EMIT
       input (src.calibrate.sweep_tau_tr's benzene-EMIT sensitivity check)
       run against the mirror without a separate carve-out.
     - mol_list_method.csv / characterised_modes.csv / ref-label_citation.csv
       -- plain copies (literature/roster metadata, unaffected by G09 vs G16).
     - results/, figures/, intermediate/ -- empty, generated fresh.

2. Re-scores the full 68-molecule roster through the REAL engine
   (src.library_ingest.run_ingest_pipeline) against the mirror, reusing the
   FROZEN thresholds in data/results/thresholds.json (tau_TR=0.95,
   tau_S=0.90368, tau_B=0.17327) -- NOT recalibrated from the rerun data, so
   this is an apples-to-apples same-classifier comparison. Every
   warning/exception raised during this call (parse failures, roster/disk
   mismatches, label-join skip-reports) is captured and printed -- the
   primary signal for "are there errors in the rerun data".

3. Runs main.run_scoring_pipeline("C6H6", "normal", ...) against the mirror
   to produce C6H6_normal.csv from rerun data (the only per-molecule CSV any
   rerun-affected figure needs -- plot_benzene_normal_modes; EMIT-based
   figures reuse the canonical C6H6_EMIT.csv unmodified, see step 4).

4. Regenerates every figure whose inputs are rerun-derived, via direct
   plot_*(...) calls with path overrides (no figures.py source edit needed
   beyond plot_rigorous_tier_check's new csv_path= parameter, added
   alongside this script). EMIT-based and CPU-timing/PED figures are
   deliberately NOT regenerated (their inputs are unaffected by this rerun
   batch) -- see NOT_REGENERATED_REASONS below for the full list and why,
   including one genuine figures.py code-limitation finding
   (plot_benzene_internal_confusion cannot be redirected at all -- see its
   entry).

5. Compares canonical data/results/library_scores.csv against the mirror's
   freshly-scored data/_rerun_mirror/results/library_scores.csv, merged on
   (molecule, mode_index, kind). freq / V_Stretch / d_CA are compared
   directly (already sign-invariant by construction). Tx..Rz are compared
   BY ABSOLUTE VALUE per axis (abs(canonical) - abs(rerun)) -- a mode's
   displacement vector is only defined up to an overall sign, and some G09
   vs G16 vector components are benignly sign-flipped (confirmed by hand
   before this script was written); comparing raw signed values would
   generate false-positive deltas from that known, harmless convention
   difference. predicted_label / predicted_annotation are compared as exact
   strings -- Step 3's purity gate is itself built on |score| (see
   IMPLEMENTATION_PLAN.md's authoritative spec), so labels are expected to
   be sign-invariant and match exactly; any mismatch here is a genuine
   finding, not sign-convention noise.

6. Writes data/results/rerun_consistency_report.csv (one row per mode,
   flagged/mismatched rows sorted first) and prints a summary: total label
   mismatches, the largest deltas per numeric column (with which
   molecule/mode), and anything captured in step 2.

7. Copies the successfully-regenerated rerun figures into the TRACKED
   data/figures/rerun/ directory (mirrors canonical figure basenames) for
   visual side-by-side comparison.

Never writes into data/results/*.csv or data/figures/*.{pdf,png} other than
the two new outputs named above (rerun_consistency_report.csv,
data/figures/rerun/*) -- the pre-existing canonical files are read-only
inputs here.
"""
import os
import shutil
import subprocess
import sys
import warnings

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd

from src.classifier import Thresholds

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MIRROR = os.path.join(REPO_ROOT, "data", "_rerun_mirror")
MIRROR_RESULTS = os.path.join(MIRROR, "results")
MIRROR_FIGURES = os.path.join(MIRROR, "figures")
TRACKED_RERUN_FIGURES_DIR = os.path.join(REPO_ROOT, "data", "figures", "rerun")

# link name (inside the mirror) -> real target directory it should junction to.
_JUNCTIONS = {
    "logs": os.path.join(REPO_ROOT, "data", "logs_rerun"),
    "gjf": os.path.join(REPO_ROOT, "data", "gjf_rerun"),
    "EMIT": os.path.join(REPO_ROOT, "data", "EMIT"),  # canonical -- not rerun, see module docstring
}

_STATIC_FILES = ["mol_list_method.csv", "characterised_modes.csv", "ref-label_citation.csv"]

_TR_AXES = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz"]

# Comparison tolerances (step 5/6). freq: same 0.05 cm^-1 magnitude
# src/library_ingest.py's attach_labels() already uses for its own
# frequency-agreement gate. max_TR_abs_delta: scores live in [-1, 1], so
# 1e-3 is tight but not over-sensitive to double-precision/print-format noise.
_FREQ_TOL = 0.05
_TR_TOL = 1e-3


# --------------------------------------------------------------------------
# Step 1: build the mirror
# --------------------------------------------------------------------------

def _ensure_junction(link_path, target_path):
    """Idempotently (re)create an NTFS directory junction at `link_path`
    pointing at `target_path` (mklink /J -- no admin rights required, unlike
    a symlink). Safe to re-run: if a junction already sits at `link_path` it
    is removed (os.rmdir on a junction only detaches the link -- verified by
    hand before this script was written; it never touches the TARGET
    directory's contents) and recreated; a REAL directory sitting at
    `link_path` is left alone and raises, so this can never silently delete
    real data.
    """
    target_path = os.path.abspath(target_path)
    if not os.path.isdir(target_path):
        raise FileNotFoundError(f"junction target does not exist: {target_path}")
    if os.path.exists(link_path):
        if os.path.isjunction(link_path):
            os.rmdir(link_path)
        else:
            raise RuntimeError(
                f"{link_path} exists and is NOT a junction -- refusing to touch it "
                "(expected only scripts/compare_rerun.py-managed junctions here).")
    result = subprocess.run(
        ["cmd", "/c", "mklink", "/J", link_path, target_path],
        capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"mklink /J failed for {link_path} -> {target_path}: "
                            f"{result.stdout}{result.stderr}")


def build_mirror():
    """(Re)build data/_rerun_mirror/ -- idempotent, safe to re-run."""
    os.makedirs(MIRROR, exist_ok=True)
    for name, target in _JUNCTIONS.items():
        _ensure_junction(os.path.join(MIRROR, name), target)
    for fname in _STATIC_FILES:
        shutil.copyfile(os.path.join(REPO_ROOT, "data", fname), os.path.join(MIRROR, fname))
    for sub in ("results", "figures", "intermediate"):
        os.makedirs(os.path.join(MIRROR, sub), exist_ok=True)
    print(f"Mirror ready at {MIRROR}")
    print(f"  logs -> {_JUNCTIONS['logs']}")
    print(f"  gjf  -> {_JUNCTIONS['gjf']}")
    print(f"  EMIT -> {_JUNCTIONS['EMIT']} (canonical, reused as-is -- EMIT was not rerun)")


# --------------------------------------------------------------------------
# Step 2: real-engine rerun ingest (frozen thresholds, warnings captured)
# --------------------------------------------------------------------------

def run_rerun_ingest():
    from src.library_ingest import run_ingest_pipeline, check_roster_disk_consistency, load_mol_roster

    roster = load_mol_roster(MIRROR)
    missing, orphaned = check_roster_disk_consistency(roster, MIRROR)
    print(f"\nRoster/disk consistency (mirror): {len(missing)} missing, {len(orphaned)} orphaned "
          f"(of {len(roster)} roster rows).")
    if missing:
        print("  MISSING roster row(s) with no matching .log+.gjf pair on disk (real finding):")
        for mol, base in missing:
            print(f"    - {mol} ({base})")
    if orphaned:
        print(f"  Orphaned on-disk basename(s) (expected, out-of-roster files): {orphaned}")

    thresholds = Thresholds.calibrated()  # frozen data/results/thresholds.json -- NOT recalibrated
    print(f"Using FROZEN thresholds (not recalibrated): tau_TR={thresholds.tau_TR}, "
          f"tau_S={thresholds.tau_S}, tau_B={thresholds.tau_B}")

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        df_lib, path, skip_report = run_ingest_pipeline(data_dir=MIRROR, thresholds=thresholds)

    print(f"Rerun library_scores.csv: {len(df_lib)} rows -> {path}")
    if skip_report:
        print(f"  {len(skip_report)} molecule(s) had ref_label/ref_key join skipped "
              "(frequency mismatch vs. characterised_modes.csv) -- see warnings below.")

    print(f"\nWarnings/exceptions captured during rerun ingest ({len(caught)}):")
    if not caught:
        print("  (none)")
    for rec in caught:
        print(f"  [{rec.category.__name__}] {rec.message}")

    return df_lib, thresholds, missing, orphaned, caught


# --------------------------------------------------------------------------
# Step 3: per-molecule C6H6 normal-mode CSV from rerun data
# --------------------------------------------------------------------------

def run_rerun_c6h6_normal(thresholds):
    from main import run_scoring_pipeline
    df, path = run_scoring_pipeline("C6H6", "normal", data_dir=MIRROR, thresholds=thresholds)
    print(f"\nRerun C6H6_normal.csv -> {path} ({len(df)} rows)")
    return path


# --------------------------------------------------------------------------
# Step 4: figure regeneration (rerun-derived inputs only)
# --------------------------------------------------------------------------

# Figures deliberately NOT regenerated for the rerun, and why -- see the
# module docstring's step 4 for the summary; full reasoning per entry below.
NOT_REGENERATED_REASONS = {
    "fig_benzene (plot_benzene_stress_test)": (
        "EMIT-based (C6H6_EMIT.csv's C2_* projected-contribution columns); EMIT "
        "inputs were not part of this rerun batch -- reuses the canonical "
        "data/results/C6H6_EMIT.csv unmodified, so regenerating would just "
        "reproduce the identical canonical figure."),
    "fig_benzeneemitcounts (plot_benzene_emit_counts)": (
        "EMIT-based (C6H6_EMIT.csv's label value_counts) -- same reasoning as fig_benzene."),
    "fig_benzeneconfusion (plot_benzene_internal_confusion)": (
        "CODE LIMITATION, not a scope choice: this function takes no data_dir/"
        "library_csv/lib_df parameter -- it calls "
        "src.benzene_validation.benzene_normal_reference_detail() with NO "
        "arguments, which hardcodes data_dir='data' internally. It cannot be "
        "redirected at the mirror without a second figures.py source edit, and "
        "this task's scope authorizes exactly one (plot_rigorous_tier_check). Its "
        "SI companion, fig_benzene_precision_recall, WAS regenerated below "
        "(src.benzene_validation.run_benzene_internal_confusion() IS "
        "data_dir-parameterized) -- so the underlying precision/recall numbers "
        "ARE available for the rerun even though the heatmap figure itself is not."),
    "fig_cputime (plot_cpu_time_benchmark, linear+log)": (
        "CPU-timing benchmark, unrelated to G09-vs-G16 QM-calc consistency -- "
        "also depends on data/results/cpu_time_benchmark.csv, which was never "
        "regenerated for the rerun logs (out of this task's scope)."),
    "fig_gaussian_nbasis (plot_gaussian_nbasis_scaling, linear+log)": (
        "CPU/basis-size scaling diagnostic, unrelated to this consistency check -- "
        "same reasoning as fig_cputime."),
    "fig_ped_vs_vscore / fig_ped_vs_bondscore": (
        "VEDA4 PED-based (data/results/combined_ped_vs_scores.csv, built from "
        "data/ved/*.ved/.vdf) -- PED inputs were not part of this rerun batch."),
}


def regenerate_rerun_figures(c6h6_normal_path):
    from src.figures import (
        plot_benzene_normal_modes, plot_confusion_matrix, plot_confusion_retention_migration,
        plot_rigorous_tier_check, plot_bond_scores, plot_boxplots, plot_mode_mixing,
        plot_irrep_coupling, plot_benzene_confusion_precision_recall, plot_sensitivity,
    )
    from src.benzene_validation import run_benzene_internal_confusion
    from src.calibrate import sweep_tau_tr

    lib_csv = os.path.join(MIRROR_RESULTS, "library_scores.csv")
    cm_csv = os.path.join(MIRROR, "characterised_modes.csv")
    mol_list_csv = os.path.join(MIRROR, "mol_list_method.csv")

    built = {}
    failed = {}

    jobs = [
        ("fig_benzene_normal", lambda: plot_benzene_normal_modes(
            normal_csv=c6h6_normal_path, out_dir=MIRROR_FIGURES)),
        ("fig_confusion", lambda: plot_confusion_matrix(
            library_csv=lib_csv, out_dir=MIRROR_FIGURES)),
        ("fig_confusion_retention_migration", lambda: plot_confusion_retention_migration(
            library_csv=lib_csv, out_dir=MIRROR_FIGURES)),
        ("fig_rigorous_tier_check", lambda: plot_rigorous_tier_check(
            library_csv=lib_csv, out_dir=MIRROR_FIGURES,
            csv_path=os.path.join(MIRROR_RESULTS, "rigorous_tier_consistency_table.csv"))),
        ("fig_bondscores", lambda: plot_bond_scores(
            library_csv=lib_csv, out_dir=MIRROR_FIGURES)),
        ("fig_boxplots", lambda: plot_boxplots(
            library_csv=lib_csv, out_dir=MIRROR_FIGURES)),
        ("fig_modemixing", lambda: plot_mode_mixing(
            library_csv=lib_csv, out_dir=MIRROR_FIGURES)),
        ("fig_irrep_coupling", lambda: plot_irrep_coupling(
            characterised_modes_csv=cm_csv, library_scores_csv=lib_csv,
            mol_list_csv=mol_list_csv, out_dir=MIRROR_FIGURES)),
    ]

    print("\nRegenerating rerun-derived figures:")
    for label, fn in jobs:
        try:
            result = fn()
            built[label] = result
            print(f"  built {label} -> {result['pdf']}")
        except Exception as e:
            failed[label] = str(e)
            print(f"  [FAILED] {label}: {e}")

    # benzene internal-confusion precision/recall companion: its inputs
    # (benzene_internal_confusion_matrix/summary.csv) come from
    # src.benzene_validation.run_benzene_internal_confusion(data_dir=...),
    # which IS data_dir-parameterized -- run it against the mirror directly,
    # independent of the (non-redirectable) heatmap figure itself.
    try:
        run_benzene_internal_confusion(data_dir=MIRROR)
        result = plot_benzene_confusion_precision_recall(
            matrix_csv=os.path.join(MIRROR_RESULTS, "benzene_internal_confusion_matrix.csv"),
            summary_csv=os.path.join(MIRROR_RESULTS, "benzene_internal_confusion_summary.csv"),
            out_dir=MIRROR_FIGURES)
        built["fig_benzene_precision_recall"] = result
        print(f"  built fig_benzene_precision_recall -> {result['pdf']}")
    except Exception as e:
        failed["fig_benzene_precision_recall"] = str(e)
        print(f"  [FAILED] fig_benzene_precision_recall: {e}")

    # tau_TR sensitivity sweep: diagnostic re-run of the sweep CURVE against
    # the rerun's geometry pool + the reused canonical benzene EMIT modes,
    # using the FROZEN tau_S/tau_B (not re-derived). This is NOT a
    # recalibration of the thresholds actually used to score above -- those
    # stay frozen at data/results/thresholds.json throughout this script.
    try:
        th = Thresholds.calibrated()
        lib_df = pd.read_csv(lib_csv)
        sweep_df = sweep_tau_tr(lib_df, th.tau_S, th.tau_B, data_dir=MIRROR)
        sweep_path = os.path.join(MIRROR_RESULTS, "tau_sensitivity_sweep.csv")
        sweep_df.to_csv(sweep_path, index=False)
        thresholds_json_mirror = os.path.join(MIRROR_RESULTS, "thresholds.json")
        shutil.copyfile(os.path.join(REPO_ROOT, "data", "results", "thresholds.json"),
                         thresholds_json_mirror)
        result = plot_sensitivity(sweep_csv=sweep_path, thresholds_json=thresholds_json_mirror,
                                   out_dir=MIRROR_FIGURES)
        built["fig_sensitivity"] = result
        print(f"  built fig_sensitivity -> {result['pdf']}")
    except Exception as e:
        failed["fig_sensitivity"] = str(e)
        print(f"  [FAILED] fig_sensitivity: {e}")

    print("\nFigures NOT regenerated for the rerun (see reasons):")
    for label, reason in NOT_REGENERATED_REASONS.items():
        print(f"  [SKIPPED] {label}: {reason}")

    return built, failed


def copy_built_figures(built):
    os.makedirs(TRACKED_RERUN_FIGURES_DIR, exist_ok=True)
    copied = []
    for result in built.values():
        for key in ("pdf", "png"):
            src = result.get(key)
            if src and os.path.exists(src):
                dst = os.path.join(TRACKED_RERUN_FIGURES_DIR, os.path.basename(src))
                shutil.copyfile(src, dst)
                copied.append(dst)
    print(f"\nCopied {len(copied)} rerun figure file(s) -> {TRACKED_RERUN_FIGURES_DIR}")
    return copied


# --------------------------------------------------------------------------
# Step 5/6: comparison + report
# --------------------------------------------------------------------------

def _load_lib(path):
    return pd.read_csv(path, dtype={"mode_index": str})


def compare(canonical_path, rerun_path):
    """Merge canonical vs. rerun library_scores.csv on
    (molecule, mode_index, kind) and compute the delta/match columns
    described in the module docstring. Returns one row per union of both
    sides (an unpaired row -- present on only one side -- is itself flagged
    via `present_both=False`, not silently dropped)."""
    canon = _load_lib(canonical_path)
    rerun = _load_lib(rerun_path)
    key = ["molecule", "mode_index", "kind"]
    merged = canon.merge(rerun, on=key, how="outer", suffixes=("_canon", "_rerun"), indicator=True)

    merged["freq_delta"] = merged["freq_rerun"] - merged["freq_canon"]
    merged["V_Stretch_delta"] = merged["V_Stretch_rerun"] - merged["V_Stretch_canon"]
    merged["d_CA_delta"] = merged["d_CA_rerun"] - merged["d_CA_canon"]

    max_axis_abs_delta = pd.Series(0.0, index=merged.index)
    for ax in _TR_AXES:
        c, r = f"{ax}_canon", f"{ax}_rerun"
        d = (merged[c].abs() - merged[r].abs()).abs()
        merged[f"{ax}_abs_delta"] = d
        max_axis_abs_delta = np.maximum(max_axis_abs_delta, d.fillna(0.0))
    merged["max_TR_abs_delta"] = max_axis_abs_delta

    merged["label_match"] = merged["predicted_label_canon"] == merged["predicted_label_rerun"]
    ann_c, ann_r = merged["predicted_annotation_canon"], merged["predicted_annotation_rerun"]
    merged["annotation_match"] = (ann_c == ann_r) | (ann_c.isna() & ann_r.isna())

    merged["present_both"] = merged["_merge"] == "both"
    merged["any_mismatch"] = (
        (~merged["present_both"])
        | (~merged["label_match"])
        | (~merged["annotation_match"])
        | (merged["freq_delta"].abs() > _FREQ_TOL)
        | (merged["max_TR_abs_delta"] > _TR_TOL)
    )

    report_cols = key + [
        "present_both",
        "freq_canon", "freq_rerun", "freq_delta",
        "V_Stretch_canon", "V_Stretch_rerun", "V_Stretch_delta",
        "d_CA_canon", "d_CA_rerun", "d_CA_delta",
    ] + [f"{ax}_abs_delta" for ax in _TR_AXES] + [
        "max_TR_abs_delta",
        "predicted_label_canon", "predicted_label_rerun", "label_match",
        "predicted_annotation_canon", "predicted_annotation_rerun", "annotation_match",
        "any_mismatch",
    ]
    report = merged[report_cols].copy()
    report = report.sort_values(["any_mismatch", "max_TR_abs_delta"], ascending=[False, False])
    return report


def print_summary(report):
    n = len(report)
    n_mismatch = int(report["any_mismatch"].sum())
    n_label_mismatch = int((~report["label_match"]).sum())
    n_annotation_mismatch = int((~report["annotation_match"]).sum())
    n_unpaired = int((~report["present_both"]).sum())

    print("\n" + "=" * 70)
    print("Rerun consistency verdict")
    print("=" * 70)
    print(f"Total (molecule, mode_index, kind) rows compared: {n}")
    print(f"Rows flagged (any_mismatch): {n_mismatch}")
    print(f"  predicted_label mismatches: {n_label_mismatch}  (target: 0)")
    print(f"  predicted_annotation mismatches: {n_annotation_mismatch}")
    print(f"  rows unpaired (present in only canonical or only rerun): {n_unpaired}")

    if n_label_mismatch:
        print("\n  predicted_label mismatches (genuine finding, not sign noise):")
        for _, row in report.loc[~report["label_match"] & report["present_both"]].iterrows():
            print(f"    {row['molecule']} mode {row['mode_index']} ({row['kind']}): "
                  f"canonical={row['predicted_label_canon']!r} rerun={row['predicted_label_rerun']!r}")

    if n_unpaired:
        print("\n  Unpaired rows:")
        for _, row in report.loc[~report["present_both"]].iterrows():
            side = "canonical only" if pd.isna(row["freq_rerun"]) else "rerun only"
            print(f"    {row['molecule']} mode {row['mode_index']} ({row['kind']}): {side}")

    print("\n  Largest deltas:")
    for col, name in [("freq_delta", "freq (cm^-1)"), ("V_Stretch_delta", "V_Stretch"),
                       ("d_CA_delta", "d_CA"), ("max_TR_abs_delta", "max |T/R| axis (abs-compared)")]:
        sub = report.dropna(subset=[col])
        if sub.empty:
            print(f"    {name}: (no comparable rows)")
            continue
        idx = sub[col].abs().idxmax()
        row = sub.loc[idx]
        print(f"    {name}: {row[col]:.6g}  ({row['molecule']} mode {row['mode_index']}, {row['kind']})")


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def main():
    print("=" * 70)
    print("Rerun consistency check: G09 (canonical, data/logs+gjf) vs.")
    print("G16 (rerun, commit 33d9efe, data/logs_rerun+gjf_rerun)")
    print("=" * 70)

    build_mirror()
    df_lib_rerun, thresholds, missing, orphaned, ingest_warnings = run_rerun_ingest()
    c6h6_normal_path = run_rerun_c6h6_normal(thresholds)
    built, failed_figures = regenerate_rerun_figures(c6h6_normal_path)
    copy_built_figures(built)

    canonical_path = os.path.join(REPO_ROOT, "data", "results", "library_scores.csv")
    rerun_path = os.path.join(MIRROR_RESULTS, "library_scores.csv")
    report = compare(canonical_path, rerun_path)

    report_path = os.path.join(REPO_ROOT, "data", "results", "rerun_consistency_report.csv")
    report.to_csv(report_path, index=False, float_format="%.6f")
    print(f"\nWrote {len(report)}-row consistency report -> {report_path}")

    print_summary(report)

    if missing:
        print(f"\n*** REAL FINDING: {len(missing)} roster molecule(s) missing from the rerun "
              "disk (logs_rerun/gjf_rerun) -- see 'Roster/disk consistency' above. ***")
    if failed_figures:
        print(f"\n*** {len(failed_figures)} figure(s) FAILED to regenerate against rerun data: "
              f"{list(failed_figures.keys())} -- see [FAILED] lines above. ***")

    return report


if __name__ == "__main__":
    main()
