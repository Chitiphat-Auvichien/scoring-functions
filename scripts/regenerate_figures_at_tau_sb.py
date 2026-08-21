"""Regenerate the manuscript's ``tau_SB``-sensitive figures AND results
CSVs at an ALTERNATE ``tau_SB`` cutoff, for comparison against the
canonical ``tau_SB=0.50`` default -- without touching ``data/figures/``,
``data/results/``, or ``data/results/thresholds.json`` in any way.

Cheap by design: ``V_Stretch`` (s[V_S]) itself does not depend on
``tau_SB`` at all -- only which SIDE of an already-computed score a mode
falls on changes. Every affected figure function in ``src/figures.py``
exposes a ``tau_SB=None`` override kwarg that cheaply re-derives labels
from the already-computed ``V_Stretch`` column
(``src.classifier.vib_label_binary`` / ``rescheme_internal_label`` /
``vib_label``), and this script applies the identical pattern directly to
the underlying results CSVs (``library_scores.csv`` and everything
downstream of it). This script never re-runs ``reproduce.py``,
``--library``, ``-m <mol>``, or any Gaussian-log/geometry parsing.

Usage
-----
    python scripts/regenerate_figures_at_tau_sb.py 0.42
    python scripts/regenerate_figures_at_tau_sb.py 0.42 --out-suffix tau_sb_0.42

Produces:
  - ``data/figures_<out-suffix>/`` -- a COMPLETE, self-contained figure set
    (tau_SB-sensitive figures regenerated, everything else copied verbatim
    from ``data/figures/``).
  - ``data/results_<out-suffix>/`` -- a SUBSET of ``data/results/``: only
    the tau_SB-sensitive CSVs (``library_scores.csv``, ``C6H6_normal.csv``,
    ``combined_ped_vs_scores.csv``, 2 of the 5 ``benzene_*`` diagnostic CSV
    pairs, ``plot_transferability_confusion_binary``'s 3 companion CSVs),
    clearly named to match their canonical counterparts. Everything
    tau_SB-independent in ``data/results/`` is deliberately NOT copied here
    -- see the generated README for the full audit/reasoning.
"""
from __future__ import annotations

import argparse
import dataclasses
import os
import shutil
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src import figures as F
from src.classifier import (
    Thresholds, rescheme_internal_label, vib_label, vib_label_binary,
)

# ==========================================================================
# FIGURES
# ==========================================================================

# --------------------------------------------------------------------------
# tau_SB-SENSITIVE figures: regenerated via each function's own tau_SB=
# override kwarg. 3 already had it (the binary-scheme CANONICAL confusion
# figures); 3 were extended with it as part of this script (see
# src/figures.py's plot_benzene_normal_modes/plot_mode_mixing/
# plot_ped_vs_vscore docstrings).
# --------------------------------------------------------------------------
SENSITIVE = [
    ("fig:confusion (binary scheme, CANONICAL)", F.plot_confusion_matrix_binary),
    ("fig:transferabilityconfusion (binary scheme, CANONICAL)",
     F.plot_transferability_confusion_binary),
    ("fig:benzeneconfusion (binary scheme, CANONICAL)",
     F.plot_benzene_internal_confusion_binary),
    ("benzene-normal-modes gallery (no fig: label yet)", F.plot_benzene_normal_modes),
    ("fig:modemixing", F.plot_mode_mixing),
    ("fig:ped_vs_vscore (proposed, not yet in .tex)", F.plot_ped_vs_vscore),
]

# --------------------------------------------------------------------------
# tau_SB-INDEPENDENT figures: copied verbatim (PDF + PNG) from data/figures/
# so the alternate-tau_SB folder is a complete, self-contained set, not just
# the ones that changed. Audited by reading every function body in
# src/figures.py (not assumed) -- see this script's own module docstring
# and the accompanying README template below for the full reasoning.
#
# fig_sensitivity_binary (plot_tau_sb_sensitivity) is a deliberate special
# case: it IS built from tau_SB (it's the tau_SB error-sweep plot itself),
# but its data curve (error_all/error_test vs. every tau_SB on the grid)
# does not depend on which tau_SB happens to be "active" -- only the
# dashed "active tau_SB" vertical line's position would move. Re-annotating
# it at tau_SB=0.42 would just slide that line on top of the plot's own
# already-drawn "optimal (all)=0.42" dotted reference line (opt_all is
# read from thresholds.json's frozen tau_SB_error_sweep block regardless of
# this script's argument), which is redundant, not more informative -- so
# it is copied unchanged, not regenerated. See README below.
# --------------------------------------------------------------------------
INDEPENDENT_STEMS = [
    "fig_benzene",
    "fig_confusion_threeway",
    "fig_confusion_retention_migration",
    "fig_rigorous_tier_check",
    "fig_transferability_confusion_threeway",
    "fig_benzene_confusion_threeway",
    "fig_benzene_precision_recall",
    "fig_benzene_emit_counts",
    "fig_bondscores",
    "fig_boxplots",
    "fig_irrep_coupling",
    "fig_sensitivity",
    "fig_sensitivity_binary",
    "fig_cputime",
    "fig_cputime_log",
    "fig_cputime_bw",
    "fig_gaussian_nbasis",
    "fig_gaussian_nbasis_linear",
    "fig_cputime_scaling_comparison",
    "fig_ped_vs_bondscore",
]


def regenerate_figures(tau_sb, fig_dir, results_dir):
    """Regenerate the 6 tau_SB-sensitive figures into `fig_dir` (redirecting
    plot_transferability_confusion_binary's 3 companion CSVs into
    `results_dir`), then copy the tau_SB-independent figures verbatim from
    data/figures/. Returns (regenerated, copied, skipped) -- lists of
    (tex_label, pdf_path) / stem / stem respectively.
    """
    regenerated = []
    for tex_label, fn in SENSITIVE:
        kwargs = dict(out_dir=fig_dir, tau_SB=tau_sb)
        if fn is F.plot_transferability_confusion_binary:
            os.makedirs(results_dir, exist_ok=True)
            kwargs["matrix_csv_path"] = os.path.join(
                results_dir, "transferability_confusion_matrix.csv")
            kwargs["summary_csv_path"] = os.path.join(
                results_dir, "transferability_confusion_summary.csv")
            kwargs["misclassified_csv_path"] = os.path.join(
                results_dir, "transferability_confusion_misclassified.csv")
        result = fn(**kwargs)
        regenerated.append((tex_label, result["pdf"]))
        print(f"  [regenerated] {tex_label}")
        print(f"                -> {result['pdf']}")

    print()
    copied, skipped = [], []
    canon_dir = "data/figures"
    for stem in INDEPENDENT_STEMS:
        ok = True
        for ext in ("pdf", "png"):
            src_path = os.path.join(canon_dir, f"{stem}.{ext}")
            dst_path = os.path.join(fig_dir, f"{stem}.{ext}")
            if not os.path.exists(src_path):
                print(f"  [WARNING] missing canonical source, skipped: {src_path}")
                ok = False
                continue
            shutil.copy2(src_path, dst_path)
        (copied if ok else skipped).append(stem)
    print(f"  [copied unchanged] {len(copied)} tau_SB-independent figures "
          f"({len(skipped)} skipped/missing)")
    return regenerated, copied, skipped


# ==========================================================================
# RESULTS (data/results/*.csv)
# ==========================================================================

# --------------------------------------------------------------------------
# Audit of data/results/ (beyond the figure-adjacent CSVs already handled
# above), by reading src/library_ingest.py, src/benzene_validation.py, and
# every to_csv() call site in src/figures.py:
#
# tau_SB-SENSITIVE (regenerated below):
#   - library_scores.csv -- THE master table. predicted_label (internal
#     rows) and predicted_annotation's embedded "vibration=<label>" (any
#     mixed-external row) are both derived from vib_label(..., scheme=
#     "binary") under the hood -- both are re-derived here from the
#     already-computed V_Stretch column.
#   - C6H6_normal.csv -- see NORMAL_CSV_SCOPE_NOTE below for why only this
#     one molecule's dump is regenerated, not all 26.
#   - combined_ped_vs_scores.csv -- same vib_label_binary(V_Stretch, tau_SB)
#     re-derivation plot_ped_vs_vscore already used transiently; written
#     out as a first-class CSV here too, full canonical row scope (not the
#     figure's mol_type=='test'-only filter).
#   - benzene_normal_reference_detail.csv / _summary.csv
#     (src.benzene_validation.benzene_normal_reference_detail/_summary):
#     default scheme=None "tracks the cached CSV's own scheme" (binary) --
#     genuinely tau_SB-sensitive.
#   - benzene_worked_examples.csv (benzene_validation.benzene_worked_examples):
#     no scheme param of its own -- reads predicted_label straight off
#     whatever lib_df it's given, so it inherits tau_SB-sensitivity from the
#     already-reschemed library_scores.csv dataframe passed in below.
#   - transferability_confusion_{matrix,summary,misclassified}.csv --
#     already redirected into results_dir by
#     plot_transferability_confusion_binary itself (see regenerate_figures
#     above); listed here for completeness of the results-folder audit.
#
# tau_SB-INDEPENDENT (deliberately NOT copied -- see item 7 of the task:
# the results folder is a SUBSET, not a full mirror, unlike the figures
# folder):
#   - benzene_internal_confusion_matrix.csv / _summary.csv
#     (benzene_internal_confusion_matrix): default scheme="threeway" (tau_S/
#     tau_B, unrelated to tau_SB).
#   - benzene_mixed_bond_diagnostic.csv / benzene_mixed_degenerate_pairs.csv
#     (benzene_mixed_bond_diagnostic): default scheme="threeway" -- in fact
#     INHERENTLY three-way-only (raises ValueError under binary, since
#     MIXED_STRETCH_BEND never occurs there).
#   - benzene_sb_vs_stretch_bond_diagnostic.csv
#     (benzene_sb_vs_stretch_bond_diagnostic): default scheme="threeway",
#     same reasoning.
#   - rigorous_tier_consistency_table.csv (plot_rigorous_tier_check): built
#     off tau_S/tau_B, unrelated to tau_SB.
#   - transferability_confusion_threeway_{matrix,summary,misclassified}.csv
#     (plot_transferability_confusion, the THREEWAY sibling): tau_S/tau_B
#     scheme, unrelated to tau_SB.
#   - every other file in data/results/ (per-molecule *_normal.csv for the
#     other 25 roster molecules, *_full_ped_table.csv, cpu_time_benchmark.csv,
#     tau_sb_sensitivity_sweep.csv, tau_sensitivity_sweep.csv, the EMIT
#     CSVs, thresholds.json itself, etc.) carries no tau_SB-derived column
#     at all.
# --------------------------------------------------------------------------

NORMAL_CSV_SCOPE_NOTE = (
    "Only `C6H6_normal.csv` is regenerated among the 26 roster molecules' "
    "`<mol>_normal.csv` files -- it is the only one read by a tau_SB-"
    "sensitive figure (`plot_benzene_normal_modes`). The other 25 "
    "`_normal.csv` files are per-molecule diagnostic dumps not consumed by "
    "anything tau_SB-sensitive; regenerating all 26 would add 25 files "
    "nothing here actually reads (see `data/results/` for the canonical "
    "originals, all built at tau_SB=0.50)."
)

RESULTS_SKIPPED_NOTE_ITEMS = [
    "benzene_internal_confusion_matrix.csv / benzene_internal_confusion_summary.csv "
    "(built under scheme=\"threeway\" -- tau_S/tau_B, unrelated to tau_SB)",
    "benzene_mixed_bond_diagnostic.csv / benzene_mixed_degenerate_pairs.csv "
    "(scheme=\"threeway\"; inherently three-way-only -- MIXED_STRETCH_BEND "
    "never occurs under binary)",
    "benzene_sb_vs_stretch_bond_diagnostic.csv (scheme=\"threeway\", same reasoning)",
    "rigorous_tier_consistency_table.csv (tau_S/tau_B, unrelated to tau_SB)",
    "transferability_confusion_threeway_{matrix,summary,misclassified}.csv "
    "(the THREEWAY sibling figure's own companion CSVs -- tau_S/tau_B scheme)",
    "every other data/results/ file (per-molecule *_normal.csv for the other "
    "25 roster molecules, *_full_ped_table.csv, cpu_time_benchmark.csv, "
    "tau_sb_sensitivity_sweep.csv, tau_sensitivity_sweep.csv, the EMIT CSVs, "
    "thresholds.json itself, etc.) -- no tau_SB-derived column at all",
]


def reschemed_library_scores(tau_sb, canonical_csv="data/results/library_scores.csv"):
    """Cheaply re-derive library_scores.csv's tau_SB-dependent columns at an
    alternate tau_SB, from the already-computed (threshold-independent)
    V_Stretch column -- no re-scoring, no re-parsing.

    predicted_label: re-derived for internal (stretch/bend/mixed) rows via
    rescheme_internal_label (external T/R rows pass through unchanged, since
    tau_SB never touches Step 2/3).

    predicted_annotation: any row whose annotation already reads
    "vibration=<label>" (a Step-3 mixed-external row) has that embedded
    sub-label re-derived too, via the same vib_label(..., scheme="binary")
    call classify_all_modes itself uses to build it. This is currently a
    no-op for the canonical library (0 mixed-external rows at the current
    tau_TR=0.95), but kept correct rather than silently stale.
    """
    df = pd.read_csv(canonical_csv)
    th = dataclasses.replace(Thresholds.calibrated(), tau_SB=tau_sb)

    df["predicted_label"] = df.apply(
        lambda row: rescheme_internal_label(
            row["predicted_label"], row["V_Stretch"], th, scheme="binary"),
        axis=1)

    def _reannotate(row):
        ann = row["predicted_annotation"]
        if isinstance(ann, str) and ann.startswith("vibration="):
            return f"vibration={vib_label(row['V_Stretch'], th, scheme='binary')}"
        return ann

    df["predicted_annotation"] = df.apply(_reannotate, axis=1)
    return df


def regenerate_results(tau_sb, results_dir):
    """Regenerate the tau_SB-sensitive results CSVs (everything audited
    above) into `results_dir`. Returns a list of (description, path)
    written, in the order written. Does NOT touch data/results/.
    """
    os.makedirs(results_dir, exist_ok=True)
    written = []
    th = dataclasses.replace(Thresholds.calibrated(), tau_SB=tau_sb)

    # 1. library_scores.csv -- the master table.
    lib_df = reschemed_library_scores(tau_sb)
    lib_out = os.path.join(results_dir, "library_scores.csv")
    lib_df.to_csv(lib_out, index=False)
    written.append(("master table: predicted_label/predicted_annotation re-derived", lib_out))
    print(f"  [regenerated] library_scores.csv -> {lib_out}")

    # 2. C6H6_normal.csv only -- see NORMAL_CSV_SCOPE_NOTE.
    normal_df = pd.read_csv("data/results/C6H6_normal.csv")
    normal_df["label"] = normal_df.apply(
        lambda row: rescheme_internal_label(
            row["label"], row["V_Stretch"], th, scheme="binary"),
        axis=1)
    normal_out = os.path.join(results_dir, "C6H6_normal.csv")
    normal_df.to_csv(normal_out, index=False)
    written.append(("only roster molecule regenerated -- see scope note below", normal_out))
    print(f"  [regenerated] C6H6_normal.csv -> {normal_out}")

    # 3. combined_ped_vs_scores.csv -- full canonical row scope (not just
    # the mol_type=='test' subset plot_ped_vs_vscore filters to).
    ped_df = pd.read_csv("data/results/combined_ped_vs_scores.csv")
    ped_df["label"] = ped_df["V_Stretch"].map(lambda v: vib_label_binary(v, tau_sb))
    ped_out = os.path.join(results_dir, "combined_ped_vs_scores.csv")
    ped_df.to_csv(ped_out, index=False)
    written.append(("PED-vs-V_Stretch comparison, full canonical row scope", ped_out))
    print(f"  [regenerated] combined_ped_vs_scores.csv -> {ped_out}")

    # 4. Benzene diagnostics -- only the 2 tau_SB-sensitive pairs (see the
    # audit comment block above for why the other 3 pairs are skipped).
    from src.benzene_validation import (
        benzene_normal_reference_detail, benzene_normal_reference_summary,
        benzene_worked_examples,
    )

    detail = benzene_normal_reference_detail(lib_df=lib_df, scheme=None)
    summary = benzene_normal_reference_summary(detail)
    detail_out = os.path.join(results_dir, "benzene_normal_reference_detail.csv")
    summary_out = os.path.join(results_dir, "benzene_normal_reference_summary.csv")
    detail.to_csv(detail_out, index=False)
    summary.to_csv(summary_out, index=False)
    written.append(("per-mode reference-vs-predicted detail", detail_out))
    written.append(("per-category recall summary", summary_out))
    print(f"  [regenerated] benzene_normal_reference_detail.csv -> {detail_out}")
    print(f"  [regenerated] benzene_normal_reference_summary.csv -> {summary_out}")

    worked = benzene_worked_examples(lib_df=lib_df)
    worked_out = os.path.join(results_dir, "benzene_worked_examples.csv")
    worked.to_csv(worked_out, index=False)
    written.append(("worked-example gallery mode picks (inherits reschemed lib_df)", worked_out))
    print(f"  [regenerated] benzene_worked_examples.csv -> {worked_out}")

    # 5. thresholds_active.json -- small provenance note (not a mirror of
    # the canonical thresholds.json, which stays untouched).
    import json
    active_path = os.path.join(results_dir, "thresholds_active.json")
    with open(active_path, "w") as f:
        json.dump({
            "tau_SB": tau_sb,
            "tau_TR": th.tau_TR, "tau_S": th.tau_S, "tau_B": th.tau_B,
            "v_weighting": th.v_weighting,
            "note": ("Active thresholds for this exploratory alternate-"
                     "tau_SB run. tau_SB is the only value overridden from "
                     "Thresholds.calibrated(); tau_TR/tau_S/tau_B/"
                     "v_weighting are unchanged from canonical. The "
                     "canonical default (tau_SB=0.50) is untouched -- see "
                     "data/results/thresholds.json."),
        }, f, indent=2)
    written.append(("provenance note (tau_SB override only; not a mirror of thresholds.json)", active_path))
    print(f"  [wrote] thresholds_active.json -> {active_path}")

    return written


# ==========================================================================
# README (shared template, written into both output folders)
# ==========================================================================

README_TEMPLATE = """\
# Alternate-tau_SB output set (tau_SB={tau_sb:g})

This is an EXPLORATORY alternate-threshold output set (figures + results),
generated at tau_SB={tau_sb:g} instead of the manuscript's canonical
default, for comparison only.

**The canonical default remains tau_SB=0.50, unchanged** -- see
`data/figures/`, `data/results/`, and `data/results/thresholds.json`.
Nothing under those paths was touched by generating this output.

## Why this value

{source}

## How this was built

`V_Stretch` (s[V_S]) does not depend on tau_SB at all -- tau_SB only
decides which side of an already-computed score a mode falls on. Every
tau_SB-sensitive output below (figures and results CSVs alike) was
therefore regenerated via the cheap V_Stretch-based re-labeling path
(`src.classifier.vib_label_binary` / `rescheme_internal_label` /
`vib_label`) -- **not** by re-running `reproduce.py`, `--library`,
`-m <mol>`, or any Gaussian-log/geometry parsing.

## Figures (`data/figures_{suffix}/`)

A COMPLETE, self-contained figure set: every tau_SB-sensitive figure is
regenerated (via a `tau_SB=` override kwarg on the relevant `src.figures`
plotting function), and every tau_SB-independent figure is copied verbatim
from `data/figures/`, so this folder is not just the ones that changed.

One figure, `fig_sensitivity_binary` (the tau_SB error-sweep plot itself),
is a deliberate exception within the "copied unchanged" list: its curve
does not depend on which tau_SB is "active", and re-annotating its active-
tau_SB marker at {tau_sb:g} would only slide that line onto its own
already-drawn "optimal (all)" reference line -- so it is copied as-is
rather than regenerated with a redundant re-annotation.

### Regenerated (tau_SB-sensitive, n={n_fig_regenerated})

{fig_regenerated_list}

### Copied unchanged (tau_SB-independent, n={n_fig_copied})

{fig_copied_list}

## Results (`data/results_{suffix}/`)

Unlike the figures folder above, the results folder is a SUBSET of
`data/results/`, not a full mirror: only outputs that carry a tau_SB-
derived column are regenerated here, clearly named to match their
canonical counterparts. A full mirror of `data/results/` (which also holds
per-molecule geometry/PED dumps, CPU benchmarks, and other tau_SB-
independent files) would just duplicate untouched files and bloat the repo
for no benefit.

{normal_csv_scope_note}

### Regenerated (tau_SB-sensitive, n={n_results_regenerated})

{results_regenerated_list}

### Not included here (tau_SB-independent -- see `data/results/` for the canonical files)

{results_skipped_list}
"""


def write_readme(tau_sb, suffix, source, fig_regenerated, fig_copied,
                  results_written, fig_dir, results_dir):
    fig_regenerated_list = "\n".join(
        f"- `{lbl}` -> `{os.path.basename(pdf)}`" for lbl, pdf in fig_regenerated)
    fig_copied_list = "\n".join(f"- `{stem}.pdf` / `{stem}.png`" for stem in fig_copied)
    results_regenerated_list = "\n".join(
        f"- `{os.path.basename(path)}` -- {desc}" for desc, path in results_written)
    results_skipped_list = "\n".join(f"- {item}" for item in RESULTS_SKIPPED_NOTE_ITEMS)

    text = README_TEMPLATE.format(
        tau_sb=tau_sb, suffix=suffix, source=source,
        n_fig_regenerated=len(fig_regenerated), fig_regenerated_list=fig_regenerated_list,
        n_fig_copied=len(fig_copied), fig_copied_list=fig_copied_list,
        normal_csv_scope_note=NORMAL_CSV_SCOPE_NOTE,
        n_results_regenerated=len(results_written), results_regenerated_list=results_regenerated_list,
        results_skipped_list=results_skipped_list,
    )
    # Same content written into both output folders (single source of
    # truth, not a divergent second document) so each is discoverable on
    # its own.
    paths = []
    for out_dir in (fig_dir, results_dir):
        readme_path = os.path.join(out_dir, "README.md")
        with open(readme_path, "w") as f:
            f.write(text)
        paths.append(readme_path)
    return paths


def main():
    p = argparse.ArgumentParser(
        description="Regenerate the tau_SB-sensitive manuscript figures and "
                     "results CSVs at an alternate tau_SB cutoff, without "
                     "touching the canonical data/figures/ or data/results/.")
    p.add_argument("tau_sb", type=float, help="alternate tau_SB cutoff, e.g. 0.42")
    p.add_argument("--out-suffix", default=None,
                    help="output folder suffix (default: tau_sb_<value>, "
                         "e.g. tau_sb_0.42)")
    p.add_argument("--source", default=None,
                    help="human-readable provenance note for the generated "
                         "README (default: points at thresholds.json's "
                         "tau_SB_error_sweep.optimal_tau_all)")
    args = p.parse_args()

    tau_sb = args.tau_sb
    suffix = args.out_suffix or f"tau_sb_{tau_sb:g}"
    source = args.source or (
        "the all-molecules classification-error-minimizing value from "
        "data/results/thresholds.json's tau_SB_error_sweep block "
        "(optimal_tau_all, min_error_all ~= 2.28%, n=788 -- over all "
        "single-centre + test molecules, excluding only benzene/C6H6)."
    )

    fig_dir = os.path.join("data", f"figures_{suffix}")
    results_dir = os.path.join("data", f"results_{suffix}")
    os.makedirs(fig_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)

    print(f"Regenerating alternate-tau_SB output set: tau_SB={tau_sb}, "
          f"out-suffix={suffix!r}")
    print(f"  figures  -> {fig_dir}/")
    print(f"  results  -> {results_dir}/")
    print()

    print("== Figures ==")
    fig_regenerated, fig_copied, fig_skipped = regenerate_figures(tau_sb, fig_dir, results_dir)

    print()
    print("== Results ==")
    results_written = regenerate_results(tau_sb, results_dir)

    print()
    readme_paths = write_readme(tau_sb, suffix, source, fig_regenerated, fig_copied,
                                 results_written, fig_dir, results_dir)
    for rp in readme_paths:
        print(f"  wrote {rp}")

    print()
    print(f"Done. Figures: regenerated {len(fig_regenerated)}, copied {len(fig_copied)} "
          f"unchanged. Results: regenerated {len(results_written)} tau_SB-sensitive files.")
    print(f"Output: {fig_dir}/ , {results_dir}/")


if __name__ == "__main__":
    main()
