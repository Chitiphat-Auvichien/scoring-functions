# Alternate-tau_SB figure set (tau_SB=0.42)

This is an EXPLORATORY alternate-threshold figure set, generated at
tau_SB=0.42 instead of the manuscript's canonical default, for visual
comparison only.

**The canonical default remains tau_SB=0.50, unchanged** -- see
`data/figures/` and `data/results/thresholds.json`. Nothing under
`data/figures/`, `data/results/`, or `thresholds.json` was touched by
generating this folder.

## Why this value

the all-molecules classification-error-minimizing value from data/results/thresholds.json's tau_SB_error_sweep block (optimal_tau_all, min_error_all ~= 2.28%, n=788 -- over all single-centre + test molecules, excluding only benzene/C6H6).

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
tau_SB marker at 0.42 would only slide that line onto its own
already-drawn "suggested (all)" reference line -- so it is copied as-is
rather than regenerated with a redundant re-annotation.

## Regenerated (tau_SB-sensitive, n=6)

- `fig:confusion (binary scheme, CANONICAL)` -> `fig_confusion.pdf`
- `fig:transferabilityconfusion (binary scheme, CANONICAL)` -> `fig_transferability_confusion.pdf`
- `fig:benzeneconfusion (binary scheme, CANONICAL)` -> `fig_benzene_confusion.pdf`
- `benzene-normal-modes gallery (no fig: label yet)` -> `fig_benzene_normal.pdf`
- `fig:modemixing` -> `fig_modemixing.pdf`
- `fig:ped_vs_vscore (proposed, not yet in .tex)` -> `fig_ped_vs_vscore.pdf`

## Copied unchanged (tau_SB-independent, n=20)

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
