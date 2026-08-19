"""Regenerate the manuscript's ``tau_SB``-sensitive figures at an ALTERNATE
``tau_SB`` cutoff, for visual comparison against the canonical
``tau_SB=0.50`` default -- without touching ``data/figures/``,
``data/results/``, or ``data/results/thresholds.json`` in any way.

Cheap by design: ``V_Stretch`` (s[V_S]) itself does not depend on
``tau_SB`` at all -- only which SIDE of an already-computed score a mode
falls on changes. Every affected figure function in ``src/figures.py``
already exposes (or, for three of them, was just given) a ``tau_SB=None``
override kwarg that cheaply re-derives labels from the already-computed
``V_Stretch`` column (``src.classifier.vib_label_binary`` /
``rescheme_internal_label``). This script never re-runs ``reproduce.py``,
``--library``, ``-m <mol>``, or any Gaussian-log/geometry parsing.

Usage
-----
    python scripts/regenerate_figures_at_tau_sb.py 0.42
    python scripts/regenerate_figures_at_tau_sb.py 0.42 --out-suffix tau_sb_0.42

Produces ``data/figures_<out-suffix>/`` (a complete, self-contained figure
set -- tau_SB-sensitive figures regenerated, everything else copied
verbatim from ``data/figures/``) and, for
``plot_transferability_confusion_binary``'s companion CSVs only,
``data/results_<out-suffix>/``.
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src import figures as F

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
# already-drawn "suggested (all)=0.42" dotted reference line (opt_all is
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

README_TEMPLATE = """\
# Alternate-tau_SB figure set (tau_SB={tau_sb:g})

This is an EXPLORATORY alternate-threshold figure set, generated at
tau_SB={tau_sb:g} instead of the manuscript's canonical default, for visual
comparison only.

**The canonical default remains tau_SB=0.50, unchanged** -- see
`data/figures/` and `data/results/thresholds.json`. Nothing under
`data/figures/`, `data/results/`, or `thresholds.json` was touched by
generating this folder.

## Why this value

{source}

## How this was built

`V_Stretch` (s[V_S]) does not depend on tau_SB at all -- tau_SB only
decides which side of an already-computed score a mode falls on. Every
tau_SB-sensitive figure below was therefore regenerated via the cheap
V_Stretch-based re-labeling path (`src.classifier.vib_label_binary` /
`rescheme_internal_label`, exposed as a `tau_SB=` override kwarg on the
relevant `src.figures` plotting functions) -- **not** by re-running
`reproduce.py`, `--library`, `-m <mol>`, or any Gaussian-log/geometry
parsing. Every other (tau_SB-independent) figure is copied verbatim from
`data/figures/`, unregenerated, so this is a complete, self-contained
figure set.

One figure, `fig_sensitivity_binary` (the tau_SB error-sweep plot itself),
is a deliberate exception within the "copied unchanged" list: its curve
does not depend on which tau_SB is "active", and re-annotating its active-
tau_SB marker at {tau_sb:g} would only slide that line onto its own
already-drawn "suggested (all)" reference line -- so it is copied as-is
rather than regenerated with a redundant re-annotation.

## Regenerated (tau_SB-sensitive, n={n_regenerated})

{regenerated_list}

## Copied unchanged (tau_SB-independent, n={n_copied})

{copied_list}
"""


def main():
    p = argparse.ArgumentParser(
        description="Regenerate the tau_SB-sensitive manuscript figures at "
                     "an alternate tau_SB cutoff, without touching the "
                     "canonical data/figures/ or data/results/.")
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

    print(f"Regenerating alternate-tau_SB figure set: tau_SB={tau_sb}, "
          f"out-suffix={suffix!r}")
    print(f"  figures  -> {fig_dir}/")
    print(f"  results  -> {results_dir}/ (companion CSVs only, if any)")
    print()

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
    copied = []
    skipped = []
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

    readme_path = os.path.join(fig_dir, "README.md")
    regenerated_list = "\n".join(f"- `{lbl}` -> `{os.path.basename(pdf)}`"
                                  for lbl, pdf in regenerated)
    copied_list = "\n".join(f"- `{stem}.pdf` / `{stem}.png`" for stem in copied)
    with open(readme_path, "w") as f:
        f.write(README_TEMPLATE.format(
            tau_sb=tau_sb, source=source,
            n_regenerated=len(regenerated), regenerated_list=regenerated_list,
            n_copied=len(copied), copied_list=copied_list,
        ))
    print(f"  wrote {readme_path}")

    print()
    print(f"Done. Regenerated {len(regenerated)} tau_SB-sensitive figures, "
          f"copied {len(copied)} tau_SB-independent figures unchanged.")
    print(f"Output: {fig_dir}/ (+ {results_dir}/ companion CSVs)")


if __name__ == "__main__":
    main()
