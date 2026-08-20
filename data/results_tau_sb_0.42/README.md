# tau_SB=0.42 output set (NOW CANONICAL as of 2026-08-20)

This folder was originally built as an EXPLORATORY alternate-threshold
output set (figures + results), generated at tau_SB=0.42 for comparison
against the manuscript's then-canonical tau_SB=0.50 default.

**As of 2026-08-20, tau_SB=0.42 IS the canonical default** -- see the
top-level `"tau_SB": 0.42` key in `data/results/thresholds.json` and the
2026-08-20 entry in `IMPLEMENTATION_PLAN.md`'s "Recent history". The value
that was "alternate" here is now the default; tau_SB=0.50 is the one that
is no longer active anywhere.

## Current relationship to `data/results/` and `data/figures/`

- **This folder (`data/results_tau_sb_0.42/`) is now REDUNDANT with
  `data/results/`** for the 9 tau_SB-sensitive CSVs listed below --
  `data/results/`'s canonical copies were overwritten with these exact same
  tau_SB=0.42 values on 2026-08-20. Kept on disk for its self-contained
  history and the standalone `thresholds_active.json` provenance note, not
  because it is the only place to find these numbers anymore.
- **`data/figures_tau_sb_0.42/` has NOT yet been mirrored back into
  `data/figures/`** -- the 6 tau_SB-sensitive figures listed below
  (`fig_confusion`, `fig_transferability_confusion`, `fig_benzene_confusion`,
  `fig_benzene_normal`, `fig_modemixing`, `fig_ped_vs_vscore`) in
  `data/figures/` still reflect the OLD tau_SB=0.50 default and have not
  been regenerated. `data/figures_tau_sb_0.42/` is currently the only place
  holding up-to-date (canonical-threshold) versions of those 6 figures.
  Regenerating `data/figures/` to match is an open follow-up.

## Why this value

the all-molecules classification-error-minimizing value from data/results/thresholds.json's tau_SB_error_sweep block (optimal_tau_all, min_error_all ~= 2.28%, n=788 -- over all single-centre + test molecules, excluding only benzene/C6H6).

## How this was built

`V_Stretch` (s[V_S]) does not depend on tau_SB at all -- tau_SB only
decides which side of an already-computed score a mode falls on. Every
tau_SB-sensitive output below (figures and results CSVs alike) was
therefore regenerated via the cheap V_Stretch-based re-labeling path
(`src.classifier.vib_label_binary` / `rescheme_internal_label` /
`vib_label`) -- **not** by re-running `reproduce.py`, `--library`,
`-m <mol>`, or any Gaussian-log/geometry parsing.

## Figures (`data/figures_tau_sb_0.42/`)

A COMPLETE, self-contained figure set: every tau_SB-sensitive figure is
regenerated (via a `tau_SB=` override kwarg on the relevant `src.figures`
plotting function), and every tau_SB-independent figure is copied verbatim
from `data/figures/`, so this folder is not just the ones that changed.

One figure, `fig_sensitivity_binary` (the tau_SB error-sweep plot itself),
is a deliberate exception within the "copied unchanged" list: its curve
does not depend on which tau_SB is "active", and re-annotating its active-
tau_SB marker at 0.42 would only slide that line onto its own
already-drawn "suggested (all)" reference line -- so it is copied as-is
rather than regenerated with a redundant re-annotation.

### Regenerated (tau_SB-sensitive, n=6)

- `fig:confusion (binary scheme, CANONICAL)` -> `fig_confusion.pdf`
- `fig:transferabilityconfusion (binary scheme, CANONICAL)` -> `fig_transferability_confusion.pdf`
- `fig:benzeneconfusion (binary scheme, CANONICAL)` -> `fig_benzene_confusion.pdf`
- `benzene-normal-modes gallery (no fig: label yet)` -> `fig_benzene_normal.pdf`
- `fig:modemixing` -> `fig_modemixing.pdf`
- `fig:ped_vs_vscore (proposed, not yet in .tex)` -> `fig_ped_vs_vscore.pdf`

### Copied unchanged (tau_SB-independent, n=20)

- `fig_benzene.pdf` / `fig_benzene.png`
- `fig_confusion_threeway.pdf` / `fig_confusion_threeway.png`
- `fig_confusion_retention_migration.pdf` / `fig_confusion_retention_migration.png`
- `fig_rigorous_tier_check.pdf` / `fig_rigorous_tier_check.png`
- `fig_transferability_confusion_threeway.pdf` / `fig_transferability_confusion_threeway.png`
- `fig_benzene_confusion_threeway.pdf` / `fig_benzene_confusion_threeway.png`
- `fig_benzene_precision_recall.pdf` / `fig_benzene_precision_recall.png`
- `fig_benzene_emit_counts.pdf` / `fig_benzene_emit_counts.png`
- `fig_bondscores.pdf` / `fig_bondscores.png`
- `fig_boxplots.pdf` / `fig_boxplots.png`
- `fig_irrep_coupling.pdf` / `fig_irrep_coupling.png`
- `fig_sensitivity.pdf` / `fig_sensitivity.png`
- `fig_sensitivity_binary.pdf` / `fig_sensitivity_binary.png`
- `fig_cputime.pdf` / `fig_cputime.png`
- `fig_cputime_log.pdf` / `fig_cputime_log.png`
- `fig_cputime_bw.pdf` / `fig_cputime_bw.png`
- `fig_gaussian_nbasis.pdf` / `fig_gaussian_nbasis.png`
- `fig_gaussian_nbasis_linear.pdf` / `fig_gaussian_nbasis_linear.png`
- `fig_cputime_scaling_comparison.pdf` / `fig_cputime_scaling_comparison.png`
- `fig_ped_vs_bondscore.pdf` / `fig_ped_vs_bondscore.png`

## Results (`data/results_tau_sb_0.42/`)

Unlike the figures folder above, the results folder is a SUBSET of
`data/results/`, not a full mirror: only outputs that carry a tau_SB-
derived column are regenerated here, clearly named to match their
canonical counterparts. A full mirror of `data/results/` (which also holds
per-molecule geometry/PED dumps, CPU benchmarks, and other tau_SB-
independent files) would just duplicate untouched files and bloat the repo
for no benefit.

Only `C6H6_normal.csv` is regenerated among the 26 roster molecules' `<mol>_normal.csv` files -- it is the only one read by a tau_SB-sensitive figure (`plot_benzene_normal_modes`). The other 25 `_normal.csv` files are per-molecule diagnostic dumps not consumed by anything tau_SB-sensitive; regenerating all 26 would add 25 files nothing here actually reads (see `data/results/` for the canonical originals, all built at tau_SB=0.50).

### Regenerated (tau_SB-sensitive, n=7)

- `library_scores.csv` -- master table: predicted_label/predicted_annotation re-derived
- `C6H6_normal.csv` -- only roster molecule regenerated -- see scope note below
- `combined_ped_vs_scores.csv` -- PED-vs-V_Stretch comparison, full canonical row scope
- `benzene_normal_reference_detail.csv` -- per-mode reference-vs-predicted detail
- `benzene_normal_reference_summary.csv` -- per-category recall summary
- `benzene_worked_examples.csv` -- worked-example gallery mode picks (inherits reschemed lib_df)
- `thresholds_active.json` -- provenance note (tau_SB override only; not a mirror of thresholds.json)

### Not included here (tau_SB-independent -- see `data/results/` for the canonical files)

- benzene_internal_confusion_matrix.csv / benzene_internal_confusion_summary.csv (built under scheme="threeway" -- tau_S/tau_B, unrelated to tau_SB)
- benzene_mixed_bond_diagnostic.csv / benzene_mixed_degenerate_pairs.csv (scheme="threeway"; inherently three-way-only -- MIXED_STRETCH_BEND never occurs under binary)
- benzene_sb_vs_stretch_bond_diagnostic.csv (scheme="threeway", same reasoning)
- rigorous_tier_consistency_table.csv (tau_S/tau_B, unrelated to tau_SB)
- transferability_confusion_threeway_{matrix,summary,misclassified}.csv (the THREEWAY sibling figure's own companion CSVs -- tau_S/tau_B scheme)
- every other data/results/ file (per-molecule *_normal.csv for the other 25 roster molecules, *_full_ped_table.csv, cpu_time_benchmark.csv, tau_sb_sensitivity_sweep.csv, tau_sensitivity_sweep.csv, the EMIT CSVs, thresholds.json itself, etc.) -- no tau_SB-derived column at all
