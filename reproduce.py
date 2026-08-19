"""Regenerate every data/results/*.csv and data/figures/* from the inputs.

One ordered pass over the whole pipeline, so a change to the scoring definition
can be propagated without anyone having to remember the order by hand. Pure
sequencing: every stage calls the existing headless entry point, so there is no
logic here that is not already covered by that module's own tests.

    py reproduce.py                     # everything, under the default weighting
    py reproduce.py --v-weighting none  # the original unweighted definition
    py reproduce.py --dry-run           # list the stages without running them
    py reproduce.py --only figures      # one stage
    py reproduce.py --skip benchmark    # everything except one

The order is not arbitrary; these are the couplings that force it:

  * --calibrate must sit BETWEEN two library passes. library_scores.csv's
    predicted_label column depends on the thresholds, but the tau_S/tau_B
    derivation reads only V_Stretch, which does not. So: score the library to
    get the distribution, calibrate off it, then re-score so the labels agree
    with the cut points that were just frozen.
  * merge_ped_scores writes back INTO <mol>_normal.csv, so it has to follow the
    per-molecule scoring pass, not precede it.
  * projection merges into <mol>_EMIT.csv, which must already exist.
  * benzene_validation has no main.py flag of its own, and several figures read
    the CSVs only it produces -- omitting it leaves those figures silently
    stale. That is exactly why it is a named stage here.
  * the CPU-time benchmark has to precede the figures too: four of them read
    cpu_time_benchmark.csv, and since it times OUR classifier rather than
    Gaussian, a change to the per-bond arithmetic moves it as well.

Writes data/results/_run_manifest.json recording which weighting, which commit
and which stages produced the current outputs -- the provenance record for the
per-molecule CSVs, which carry no column of their own.
"""

import argparse
import json
import os
import subprocess
import sys
import time

from src.scoring import V_WEIGHTINGS, DEFAULT_V_WEIGHTING, set_v_weighting, get_v_weighting

REPO_ROOT = os.path.dirname(os.path.abspath(__file__))
MANIFEST_NAME = "_run_manifest.json"

# Molecules that get their own <mol>_normal.csv. Derived from what is already
# on disk rather than hardcoded, so the regenerated set matches the set being
# replaced; --molecules overrides it.
EMIT_MOLECULES = ["C6H6"]


def results_dir(data_dir):
    return os.path.join(data_dir, "results")


def discover_molecules(data_dir):
    """Molecules that currently have a <mol>_normal.csv, i.e. the set this run
    is expected to reproduce. Falls back to the archived copy so the list
    survives the archive step that empties data/results/."""
    for d in (results_dir(data_dir),
              os.path.join(results_dir(data_dir), "archive_unweighted")):
        if not os.path.isdir(d):
            continue
        found = sorted(f[:-len("_normal.csv")] for f in os.listdir(d)
                       if f.endswith("_normal.csv"))
        if found:
            return found
    return []


# --- stages -----------------------------------------------------------------

def stage_library(ctx):
    """Score all roster molecules -> library_scores.csv (provisional labels)."""
    from src.library_ingest import run_ingest_pipeline
    df, path, skip_report = run_ingest_pipeline(data_dir=ctx["data_dir"])
    ctx["outputs"].append(path)
    return f"{len(df)} rows -> {path}" + (
        f" ({len(skip_report)} label-join skips)" if skip_report else "")


def stage_calibrate(ctx):
    """Freeze tau_TR/tau_S/tau_B off that distribution -> thresholds.json,
    plus the advisory (non-overwriting) tau_SB error-vs-threshold sweep."""
    from src.calibrate import run_calibration_pipeline, run_tau_sb_error_analysis
    th, _result, sweep, (path_json, path_sweep) = run_calibration_pipeline(
        data_dir=ctx["data_dir"])
    ctx["outputs"] += [path_json, path_sweep]

    sb_sweep, plateau_all, plateau_test, path_sb_sweep = run_tau_sb_error_analysis(
        data_dir=ctx["data_dir"])
    ctx["outputs"].append(path_sb_sweep)
    _lo_all, _hi_all, mid_all, err_all = plateau_all
    _lo_test, _hi_test, mid_test, err_test = plateau_test

    return (f"tau_TR={th.tau_TR} tau_S={th.tau_S} tau_B={th.tau_B} "
            f"({th.v_weighting}) -> {path_json}; {len(sweep)}-row sweep; "
            f"tau_SB sweep ({len(sb_sweep)} rows) -> {path_sb_sweep}; "
            f"suggested tau_SB (ADVISORY, active default stays {th.tau_SB}): "
            f"all={mid_all} (err={err_all:.4f}), test={mid_test} (err={err_test:.4f})")


def stage_library2(ctx):
    """Re-score the library so predicted_label matches the frozen thresholds."""
    return stage_library(ctx)


def stage_molecules(ctx):
    """Per-molecule normal-mode scoring + classification CSVs."""
    from main import run_scoring_pipeline
    done = []
    for mol in ctx["molecules"]:
        _df, path = run_scoring_pipeline(mol, "normal", data_dir=ctx["data_dir"])
        ctx["outputs"].append(path)
        done.append(mol)
    return f"{len(done)} molecules: {', '.join(done)}"


def stage_emit(ctx):
    """EMIT-mode scoring CSVs."""
    from main import run_scoring_pipeline
    for mol in EMIT_MOLECULES:
        _df, path = run_scoring_pipeline(mol, "emit", data_dir=ctx["data_dir"])
        ctx["outputs"].append(path)
    return ", ".join(EMIT_MOLECULES)


def stage_projection(ctx):
    """EMIT -> normal-mode projection, merged into <mol>_EMIT.csv."""
    from main import run_projection_pipeline
    for mol in EMIT_MOLECULES:
        _df, _full, paths = run_projection_pipeline(mol, data_dir=ctx["data_dir"])
        ctx["outputs"] += list(paths)
    return ", ".join(EMIT_MOLECULES)


def stage_ped(ctx):
    """VEDA4 PED merge: back into <mol>_normal.csv, the combined table, and the
    per-molecule full PED tables."""
    from ped.merge_ped_scores import (build_combined_table, full_ped_table,
                                      _read_roster, _resolve_ved_paths, _REPO_ROOT)
    molecules = _read_roster(_REPO_ROOT)
    combined, per_molecule = build_combined_table(
        molecules, _REPO_ROOT, skip_missing=True)
    out = os.path.join(ctx["data_dir"], "results", "combined_ped_vs_scores.csv")
    combined.to_csv(out, index=False, float_format="%.4f")
    ctx["outputs"].append(out)

    # The full per-coordinate PED tables are a separate --full-table pass in
    # merge_ped_scores' own CLI, one molecule at a time; without them the
    # <mol>_full_ped_table.csv files silently go missing from a rebuild.
    n_full = 0
    for mol in molecules:
        ved, vdf, missing = _resolve_ved_paths(mol, ctx["data_dir"])
        if missing:
            continue
        table = full_ped_table(ved, vdf)
        path = os.path.join(ctx["data_dir"], "results", f"{mol}_full_ped_table.csv")
        table.to_csv(path, index=False, float_format="%.4f")
        ctx["outputs"].append(path)
        n_full += 1

    return (f"{len(per_molecule)}/{len(molecules)} merged, {len(combined)} combined "
            f"rows, {n_full} full PED tables")


def stage_benzene_validation(ctx):
    """Benzene diagnostic CSVs -- inputs to two of the manuscript figures."""
    from src import benzene_validation as bv
    paths = []
    _d, _s, p = bv.run_benzene_normal_validation(data_dir=ctx["data_dir"])
    paths.append(p)
    _ct, _pc, p = bv.run_benzene_internal_confusion(data_dir=ctx["data_dir"])
    paths.append(p)
    _bd, _pr, p = bv.run_benzene_bond_diagnostic(data_dir=ctx["data_dir"])
    paths.append(p)
    _sb, p = bv.run_benzene_sb_vs_stretch_bond_diagnostic(data_dir=ctx["data_dir"])
    paths.append(p)
    _we, p = bv.run_benzene_worked_examples(data_dir=ctx["data_dir"])
    paths.append(p)
    ctx["outputs"] += [str(p) for p in paths]
    return f"{len(paths)} diagnostic CSVs"


# The manuscript's images/ is a manual copy of data/figures/, not a live link.
# Three of the figures it actually \includegraphics are hand-annotated variants
# built in a vector editor from the generated PDFs; nothing here can rebuild
# them, so they are reported as stale rather than silently left alone.
ANNOTATED_FIGURES = {
    "fig_bondscores": "fig_bondscores_annot.pdf",
    "fig_benzene_normal": "fig_benzene_normal_annot.pdf",
    "fig_benzene_emit_counts": "fig_benzene_emit_counts_annot.pdf",
}
DEFAULT_MANUSCRIPT_IMAGES = os.path.normpath(os.path.join(
    REPO_ROOT, "..", "..", "JCC", "JCC_man_scoring", "images"))


def sync_manuscript(figures_dir, images_dir, dry_run=False):
    """Copy regenerated figures into the manuscript's images/ and report what
    moved. Returns (copied, skipped, stale_annotated)."""
    import filecmp
    import shutil

    copied, skipped, stale = [], [], []
    for name in sorted(os.listdir(figures_dir)):
        if not name.endswith((".pdf", ".png")):
            continue
        src = os.path.join(figures_dir, name)
        dst = os.path.join(images_dir, name)
        if not os.path.exists(dst):
            skipped.append(name)          # not used by the manuscript
            continue
        if filecmp.cmp(src, dst, shallow=False):
            continue
        if not dry_run:
            shutil.copy2(src, dst)
        copied.append(name)
        stem = os.path.splitext(name)[0]
        if stem in ANNOTATED_FIGURES:
            stale.append(ANNOTATED_FIGURES[stem])
    return copied, skipped, sorted(set(stale))


def stage_sync_manuscript(ctx):
    """Copy regenerated figures into the manuscript images/ (manual copy)."""
    images = ctx.get("images_dir") or DEFAULT_MANUSCRIPT_IMAGES
    if not os.path.isdir(images):
        return f"SKIPPED: no manuscript images dir at {images}"
    copied, skipped, stale = sync_manuscript(
        os.path.join(ctx["data_dir"], "figures"), images)
    ctx["outputs"] += [os.path.join(images, n) for n in copied]
    msg = f"{len(copied)} updated, {len(skipped)} not used by the manuscript"
    if stale:
        msg += ("\n    STALE, needs manual re-annotation in a vector editor "
                "(their .svg sources are in images/): " + ", ".join(stale))
    return msg


def stage_benchmark(ctx):
    """CPU-time benchmark -> cpu_time_benchmark.csv (four figures read it)."""
    from scripts.benchmark_cpu_time import main as benchmark_main
    benchmark_main()
    path = os.path.join(results_dir(ctx["data_dir"]), "cpu_time_benchmark.csv")
    ctx["outputs"].append(path)
    return path


def stage_figures(ctx):
    """Every manuscript figure (PDF + PNG), plus the tier-consistency table."""
    from src.figures import regenerate_all
    results = regenerate_all(verbose=False)
    for r in results.values():
        ctx["outputs"] += [r["pdf"], r["png"]]
    return f"{len(results)} figures"


# Order matters -- see the module docstring. 'benchmark' is opt-in via --only
# because it is slow and measures timing rather than scores.
STAGES = [
    ("library", stage_library),
    ("calibrate", stage_calibrate),
    ("library2", stage_library2),
    ("molecules", stage_molecules),
    ("emit", stage_emit),
    ("projection", stage_projection),
    ("ped", stage_ped),
    ("benzene_validation", stage_benzene_validation),
    # Must precede 'figures': four of them read cpu_time_benchmark.csv, and it
    # measures OUR classifier, so the weighting change moves it too.
    ("benchmark", stage_benchmark),
    ("figures", stage_figures),
]
# Opt-in: writes outside data/, so never part of a default run.
OPTIONAL_STAGES = [
    ("sync_manuscript", stage_sync_manuscript),
]
OPTIONAL_NAMES = {name for name, _ in OPTIONAL_STAGES}
ALL_STAGES = STAGES + OPTIONAL_STAGES
STAGE_NAMES = [name for name, _ in ALL_STAGES]


def select_stages(only=None, skip=()):
    """Stages to run, always in the canonical order. Optional stages are
    excluded from a default run and must be asked for by name."""
    chosen = []
    for name, fn in ALL_STAGES:
        wanted = name in only if only else name not in OPTIONAL_NAMES
        if wanted and name not in skip:
            chosen.append((name, fn))
    return chosen


def git_commit():
    try:
        out = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO_ROOT,
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return None


def write_manifest(ctx):
    path = os.path.join(results_dir(ctx["data_dir"]), MANIFEST_NAME)
    manifest = {
        "v_weighting": get_v_weighting(),
        "git_commit": git_commit(),
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "stages_run": ctx["stages_run"],
        "molecules": ctx["molecules"],
        "outputs": sorted(set(os.path.relpath(p, REPO_ROOT)
                              if os.path.isabs(p) else p
                              for p in ctx["outputs"])),
    }
    with open(path, "w") as f:
        json.dump(manifest, f, indent=2)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Regenerate all results and figures from the inputs.")
    ap.add_argument("--v-weighting", choices=list(V_WEIGHTINGS),
                    default=DEFAULT_V_WEIGHTING, dest="v_weighting",
                    help="V-score bond weighting to reproduce under (default: "
                         f"{DEFAULT_V_WEIGHTING}).")
    ap.add_argument("--data-dir", default="data")
    ap.add_argument("--only", nargs="+", choices=STAGE_NAMES,
                    help="Run only these stages (in the canonical order).")
    ap.add_argument("--skip", nargs="+", choices=STAGE_NAMES, default=[],
                    help="Run everything except these stages.")
    ap.add_argument("--molecules", nargs="+",
                    help="Override the per-molecule list (default: whatever "
                         "already has a <mol>_normal.csv).")
    ap.add_argument("--images-dir", default=None,
                    help="Manuscript images/ for the sync_manuscript stage "
                         f"(default: {DEFAULT_MANUSCRIPT_IMAGES}).")
    ap.add_argument("--dry-run", action="store_true",
                    help="List the stages that would run, then stop.")
    args = ap.parse_args(argv)

    set_v_weighting(args.v_weighting)

    selected = select_stages(args.only, args.skip)

    molecules = args.molecules or discover_molecules(args.data_dir)
    ctx = {"data_dir": args.data_dir, "molecules": molecules,
           "images_dir": args.images_dir, "outputs": [], "stages_run": []}

    print(f"reproduce.py: v_weighting={get_v_weighting()}, "
          f"{len(molecules)} molecule(s), {len(selected)} stage(s)")
    if args.dry_run:
        for i, (name, fn) in enumerate(selected, 1):
            print(f"  {i}. {name:20s} {(fn.__doc__ or '').strip().splitlines()[0]}")
        return 0

    t0 = time.time()
    for i, (name, fn) in enumerate(selected, 1):
        start = time.time()
        print(f"[{i}/{len(selected)}] {name} ...", flush=True)
        summary = fn(ctx)
        ctx["stages_run"].append(name)
        print(f"    {summary}  ({time.time() - start:.1f}s)", flush=True)

    path = write_manifest(ctx)
    print(f"Done in {time.time() - t0:.1f}s. Manifest -> {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
