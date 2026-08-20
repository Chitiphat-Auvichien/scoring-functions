# tau_SB=0.42 output set (CANONICAL as of 2026-08-20)

This folder was originally built as an EXPLORATORY alternate-threshold
output set (figures + results), generated at tau_SB=0.42 for comparison
against the manuscript's then-canonical tau_SB=0.50 default.

**As of 2026-08-20, tau_SB=0.42 IS the canonical default** -- see the
top-level `"tau_SB": 0.42` key in `data/results/thresholds.json` and the
2026-08-20 entries in `IMPLEMENTATION_PLAN.md`'s "Recent history". The value
that was "alternate" here is now the default; tau_SB=0.50 is the one that
is no longer active anywhere (its own retired snapshot lives in
`data/results_tau_sb_0.50/`).

## Figures (`data/figures_tau_sb_0.42/`) -- now a FULL mirror

**Updated 2026-08-20 (was a 26-file "6 regenerated + 20 copied-unchanged"
set until this session; see `data/results_tau_sb_0.42/README.md` for the
parallel story on the results side, fixed earlier the same day).** Canonical
`data/figures/` itself had not been regenerated at all since the tau_SB
switch -- it still reflected the OLD tau_SB=0.50 default for the 6
tau_SB-sensitive figures, AND (more importantly) it predated the broader
`data/results/` regeneration that happened earlier this session (the
`Thresholds.tau_SB` source-of-truth fix, the H2O EMIT staleness fix, and
several other pipeline fixes accumulated over recent commits -- see
`data/results_tau_sb_0.42/README.md` for the details). That meant even this
folder's own "6 regenerated" figures were built from the cheap
`tau_SB=`-override relabeling path on top of an already-stale
`data/results/`, not from the fully-corrected pipeline.

This session closed that gap for real: `py reproduce.py --only figures`
was run against the now fully self-consistent canonical `data/results/`,
regenerating all 26 figure stems (52 files: PDF + PNG each) in canonical
`data/figures/` from scratch. Every file in this folder is now a **verbatim,
byte-for-byte copy of canonical `data/figures/`'s top-level files** (52
files -- confirmed via `filecmp`), not a curated subset.

Comparing the new canonical figures against the OLD (pre-this-session)
canonical set (saved off before regenerating) confirms the expected
pattern: the 6 tau_SB-sensitive figures below changed the most (up to
~23% of pixels differing for `fig_benzene_confusion`, where reclassified
cells move between confusion-matrix categories), while the other 20 changed
by a much smaller amount (~1-4% of pixels, consistent with matplotlib's
inherent font-subsetting/tight-bbox non-determinism between separate
`savefig` runs plus a few benzene-EMIT and CPU-benchmark figures also
picking up other, unrelated pipeline fixes that had accumulated in
`data/results/` since `data/figures/` was last regenerated -- see the
per-figure list below). None of that reflects an error in this session's
regeneration; it reflects how stale the previous canonical `data/figures/`
had become.

- **tau_SB-sensitive (n=6):** `fig_confusion`, `fig_transferability_confusion`,
  `fig_benzene_confusion`, `fig_benzene_normal`, `fig_modemixing`,
  `fig_ped_vs_vscore` -- these move because reclassifying modes at the new
  cutoff changes which category cell/marker/count they land in.
- **tau_SB-independent (n=20):** everything else. These do not depend on
  tau_SB at all, but most also changed by a small amount relative to the OLD
  canonical figures because `data/figures/` had not been regenerated in a
  while and picked up other, unrelated fixes already present in
  `data/results/` (e.g. the CPU-time figures reflect a freshly re-run,
  inherently timing-noisy `cpu_time_benchmark.csv`; `fig_benzene` reflects
  the same C6H6 EMIT-projection pipeline that the H2O EMIT staleness fix
  touched). None of these 20 depend on the tau_SB value itself.

One figure, `fig_sensitivity_binary` (the tau_SB error-sweep plot itself),
has a curve that never depended on which tau_SB is "active" in the first
place -- its content changing here (like the rest) simply reflects the
underlying calibration re-run, not a tau_SB re-annotation.

### All figures (n=26 stems, 52 files, full mirror)

- `fig_benzene` -- `fig_benzene_confusion` (SENSITIVE) -- `fig_benzene_confusion_threeway`
- `fig_benzene_emit_counts` -- `fig_benzene_normal` (SENSITIVE) -- `fig_benzene_precision_recall`
- `fig_bondscores` -- `fig_boxplots`
- `fig_confusion` (SENSITIVE) -- `fig_confusion_retention_migration` -- `fig_confusion_threeway`
- `fig_cputime` -- `fig_cputime_bw` -- `fig_cputime_log` -- `fig_cputime_scaling_comparison`
- `fig_gaussian_nbasis` -- `fig_gaussian_nbasis_linear`
- `fig_irrep_coupling`
- `fig_modemixing` (SENSITIVE)
- `fig_ped_vs_bondscore` -- `fig_ped_vs_vscore` (SENSITIVE)
- `fig_rigorous_tier_check`
- `fig_sensitivity` -- `fig_sensitivity_binary`
- `fig_transferability_confusion` (SENSITIVE) -- `fig_transferability_confusion_threeway`

**Not included:** `data/figures/archive_threeway/` and
`data/figures/archive_unweighted/` -- these are frozen historical snapshots
along a *different* axis (scheme/weighting changes, not tau_SB) that predate
and are unrelated to this switch, exactly like `data/results/archive*/` was
excluded from `data/results_tau_sb_0.42/` for the same reason.

## Why this value

the all-molecules classification-error-minimizing value from
`data/results/thresholds.json`'s `tau_SB_error_sweep` block
(`optimal_tau_all`, `min_error_all` ~= 2.28%, n=788 -- over all
single-centre + test molecules, excluding only benzene/C6H6).

## How this was built (2026-08-20 remirror)

1. Confirmed canonical `data/results/` was already the fully-corrected,
   self-consistent tau_SB=0.42 set (fixed earlier this session -- see
   `data/results_tau_sb_0.42/README.md`).
2. Saved a copy of the then-current (stale, pre-this-step) canonical
   `data/figures/` for comparison.
3. Ran `py reproduce.py --only figures`, which calls
   `src.figures.regenerate_all()` against canonical `data/results/`,
   producing all 26 figure stems (PDF + PNG) fresh.
4. Diffed old-vs-new canonical figures (content-stream and pixel-level, not
   raw byte comparison, since matplotlib's PDF font subsetting and
   `bbox_inches="tight"` sizing are not perfectly deterministic run-to-run
   even with identical input data) to confirm the tau_SB-sensitive figures
   moved by far more than the tau_SB-independent ones.
5. Copied every resulting top-level file (52) in `data/figures/` into this
   folder verbatim (confirmed byte-identical via `filecmp`).

No manual relabeling and no `tau_SB=` override kwarg this time -- every
figure here came from the real pipeline against the real canonical results,
not a cheap re-annotation shortcut.

## Results (`data/results_tau_sb_0.42/`)

See `data/results_tau_sb_0.42/README.md` for the (already-completed, earlier
this session) parallel story on the results side -- that folder is likewise
now a full, self-contained mirror of canonical `data/results/`.
