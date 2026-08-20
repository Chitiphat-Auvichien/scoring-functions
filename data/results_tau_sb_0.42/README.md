# tau_SB=0.42 output set (CANONICAL as of 2026-08-20)

This folder was originally built as an EXPLORATORY alternate-threshold
output set (figures + results), generated at tau_SB=0.42 for comparison
against the manuscript's then-canonical tau_SB=0.50 default.

**As of 2026-08-20, tau_SB=0.42 IS the canonical default** -- see the
top-level `"tau_SB": 0.42` key in `data/results/thresholds.json` and the
2026-08-20 entry in `IMPLEMENTATION_PLAN.md`'s "Recent history". The value
that was "alternate" here is now the default; tau_SB=0.50 is the one that
is no longer active anywhere (its own retired snapshot lives in
`data/results_tau_sb_0.50/`).

## Results (`data/results_tau_sb_0.42/`) -- now a FULL mirror

**Updated 2026-08-20 (was a 9-file subset until this session).** An audit
this session found that the original 0.42 switch (`3427ff8`,
"Change canonical tau_SB default 0.50 -> 0.42...") was incomplete in two
ways: (1) it only hand-edited `data/results/thresholds.json`'s `"tau_SB"`
key without changing the actual source of truth (`Thresholds.tau_SB`'s
dataclass default in `src/classifier.py`, still 0.50), so that JSON edit was
silently wiped by the very next real `--calibrate` run (`calibrate()` in
`src/calibrate.py` never threaded a `tau_SB` kwarg through, so it always
fell back to the 0.50 class default and overwrote the whole file); and
(2) even setting aside that bug, only 9 of the tau_SB-sensitive CSVs had
ever been relabeled -- `C6H6_EMIT.csv` and every other per-molecule
EMIT/normal CSV were untouched, so canonical `data/results/` was internally
inconsistent (`thresholds.json` claimed 0.42, most files still reflected
0.50 or worse, H2O's EMIT CSVs were stale from *before* the 2026-08 binary-
scheme-default switch entirely). Both bugs are now fixed at the source
(`src/classifier.py`'s `Thresholds.tau_SB` default is genuinely 0.42;
`calibrate()` now writes `"tau_SB"` into every `thresholds.json` it
produces) and canonical `data/results/` was regenerated end-to-end via
`py reproduce.py --skip figures` (figures deliberately out of scope this
session), so it and this folder are once again in sync.

This folder is now a **complete, self-contained copy of every top-level
file in `data/results/`** (73 files: every roster molecule's `*_normal.csv`,
`*_full_ped_table.csv`, the EMIT CSVs for C6H6 and H2O
(`*_EMIT.csv`/`*_EMIT_full.csv`/`*_EMIT_full_cartesian.csv`),
`library_scores.csv`, the benzene diagnostic/reference CSVs, transferability
confusion tables, `thresholds.json`, `cpu_time_benchmark.csv`, the
sensitivity-sweep CSVs, `_run_manifest.json`, etc.) -- copied verbatim
(byte-for-byte identical, confirmed via `filecmp`), not just the
tau_SB-sensitive subset as before.

**Not included:** `data/results/archive/`, `data/results/archive_threeway/`,
`data/results/archive_unweighted/` -- these are frozen historical snapshots
along a *different* axis (scheme/weighting changes, not tau_SB) that predate
and are unrelated to this switch; mirroring them here would just duplicate
already-archived, unrelated history. If a tau_SB=0.50 comparison point is
ever needed for one of those, `data/results_tau_sb_0.50/` (the retired
snapshot, frozen at commit `7c18326`) is the place to look, not here.

`thresholds_active.json` (the old provenance note for this folder, back when
it held an "override" not the canonical value) has been removed -- it was
identical in substance to what's now in the mirrored `thresholds.json`, so
keeping both was redundant.

## Why this value

the all-molecules classification-error-minimizing value from data/results/thresholds.json's tau_SB_error_sweep block (optimal_tau_all, min_error_all ~= 2.28%, n=788 -- over all single-centre + test molecules, excluding only benzene/C6H6).

## How this was built (2026-08-20 remirror)

1. Fixed the root cause: `src/classifier.py`'s `Thresholds.tau_SB` dataclass
   default 0.50 -> 0.42 (the actual, sole source of truth for this constant
   -- no sweep ever computes a replacement for it, unlike tau_TR/tau_S/tau_B);
   `src/calibrate.py`'s `calibrate()` now writes `"tau_SB"` into the
   `thresholds.json` it produces (previously omitted entirely).
2. Ran `py reproduce.py --skip figures` (full pipeline: library ingest,
   calibration, per-molecule scoring, EMIT scoring + projection for C6H6 AND
   H2O (H2O added to `reproduce.py`'s `EMIT_MOLECULES` this session --
   previously hardcoded to `["C6H6"]` only, which is exactly how H2O's EMIT
   CSVs went stale), PED merge, benzene diagnostics, CPU benchmark) against
   canonical `data/results/`.
3. Copied every resulting top-level file in `data/results/` into this folder
   verbatim.

No manual relabeling this time -- every number here came from the real
engine, not a V_Stretch-based re-derivation shortcut.

## Figures (`data/figures_tau_sb_0.42/`)

Unchanged this session -- figures regeneration is a separate, still-open
follow-up (explicitly out of scope for the 2026-08-20 results-consistency
fix described above). The figures-folder documentation below is UNCHANGED
from the prior session and still describes that folder's own state, not
this one's.

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

`data/figures/` still reflects the OLD tau_SB=0.50 default for the 6
tau_SB-sensitive figures listed above and has not been regenerated to
match -- still the one open follow-up, unchanged by this session's CSV work.
