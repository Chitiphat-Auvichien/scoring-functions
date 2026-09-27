# Developer / maintainer workflow

`README.md`'s Steps 1-5 (`-m`/`--molecule` + `--mode`) are the whole job for ordinary use. `main.py` also
exposes a second group of flags — diagnostics, cross-checks, calibration, and manuscript-support
pipelines — that most users won't need. Run `python main.py --help` to see both groups. All commands
below were run from the repository root against the data already checked into this repo and confirmed to
produce real output files.

## Reproducing everything

`python reproduce.py` rebuilds every `data/results/*.csv` and `data/figures/*` from the inputs in one
ordered pass (`--dry-run` lists the stages; `--only`/`--skip` select them). It writes
`data/results/_run_manifest.json` recording which scoring definition, commit, and stages produced the
current outputs.

To rebuild under the original unweighted V-score instead, use `python reproduce.py --v-weighting none`.
Scoring and thresholds must come from the same definition: `thresholds.json` records which one it was
calibrated for, and the program refuses to classify with a mismatched pair rather than silently
mislabelling modes. The pre-2026-08-14 unweighted outputs are kept under `data/results/archive_unweighted/`
and `data/figures/archive_unweighted/`; `python scripts/compare_weighting.py` diffs the two.

## Classification thresholds (CLI overrides)

* `--tau-sb TAU_SB` — override the single binary S/B cutoff for this run only (starts from
  `Thresholds.calibrated()`, never written back to `thresholds.json`). Only affects `scheme=binary` (the
  default).
* `--classify-scheme {binary,threeway}` — Step-2 vibrational-label vocabulary. `binary` (default,
  paper-standard) forces every mode to `"S"` or `"B"` via `tau_SB`, never `"SB"`. `threeway` (explicit
  opt-in) may label a mode `"SB"` (mixed stretch-bend) via a `tau_S`/`tau_B` split.
* `--no-identify-tr` — skip Step 3 entirely; every mode's `tr_label`/`tr_score` stay blank and only
  `vib_label` is reported.
* `--tau-tr TAU_TR`, `--tau-s TAU_S`, `--tau-b TAU_B` — override the corresponding threshold for this run
  only. **Diagnostic-only as of the 2026-08-25 restructuring**: none of these gate any classification
  decision any more (Step 3 no longer gates on `tau_TR` at all; `tau_S`/`tau_B` only matter under
  `--classify-scheme threeway`). Kept for sensitivity reporting (`src/calibrate.py`'s sweeps).

## Per-molecule flags (`-m`/`--molecule` required; each needs its prerequisite CSV to already exist)

Project raw EMIT eigenvectors onto the mass-weighted normal-mode reference basis. Requires
`data/results/<mol>_EMIT.csv` to already exist (run `--mode emit` first); merges the grouped
`C2_Tx`..`C2_VMix` fractions into that file **in place** and writes `data/results/<mol>_EMIT_full.csv`
(per-reference-mode detail). Fails loudly with a clear message if a prerequisite is missing.

```bash
python main.py -m water --mode emit
python main.py -m water --emit-projection
```

Merge real VEDA4 PED (Potential Energy Distribution) output into one molecule's `<mol>_normal.csv`,
appending `PED_Stretch_pct`/`PED_Bend_pct` (+ per-bond-type columns) **in place** for validating
`V_Stretch` against ground-truth %stretch/%bend character. Requires `data/results/<mol>_normal.csv`
(from `--mode normal`) and a `data/ved/<mol>.ved` + `data/ved/<mol>.vdf` pair from a real VEDA4 GUI run
— see `ped/README.md` for the full VEDA4 bridge workflow (`ped/build_veda_fmt.py` builds the `.fmt` file
VEDA4 needs from `data/logs/<mol>.log` and, optionally, `data/fchk/<mol>.fchk`).

```bash
python main.py -m CH4 --mode normal
python main.py -m CH4 --ped-merge
```

## Global flags (molecule-independent; `-m` is ignored if supplied alongside them)

Build the hydride-library dataset from `data/logs/`+`data/gjf/` into `data/results/library_scores.csv`.
Grows automatically as more molecule files are added; can be slow once the library is large, since every
molecule is scored from scratch.

```bash
python main.py --library
```

Calibrate `tau_TR`/`tau_S`/`tau_B` against the library and run the threshold-sensitivity sweep. Writes
`data/results/thresholds.json` and `data/results/tau_sensitivity_sweep.csv`. Also runs the advisory
`tau_SB` error-vs-threshold sweep into `tau_sb_sensitivity_sweep.csv` (prints a suggested `tau_SB` but
never overwrites the frozen `tau_SB` default — see `--tau-sb` above). Rebuilds the library itself first,
so this is also slow.

```bash
python main.py --calibrate
```

Regenerate every manuscript figure from `data/results/*.csv` into `data/figures/*.{pdf,png}`.

```bash
python main.py --figures
```

Run `--ped-merge` for every molecule in `data/mol_list_method.csv`'s roster, skipping (with a message)
any missing its `<mol>_normal.csv` or `data/ved/` pair, and write one combined table across all matched
molecules for correlating `V_Stretch` against real PED %stretch (default
`data/results/combined_ped_vs_scores.csv`, override with `--combined-output`). Note: this roster is the
hydride calibration library, not any transferability-test molecule set you may be tracking separately —
use `ped/merge_ped_scores.py --molecules ...` directly (see `ped/README.md`) to merge PED for an explicit
list of molecules instead.

```bash
python main.py --ped-merge-all
```

Every manuscript figure is also a standalone function in `src/figures.py` (e.g.
`plot_benzene_stress_test`, `plot_confusion_matrix`, `plot_bond_scores`, `plot_boxplots`,
`plot_mode_mixing`, `plot_sensitivity`, `plot_benzene_normal_modes`) if you need to regenerate just one
of them from a Python session; `python main.py --figures` (or `python -m src.figures`) regenerates all
of them in one go.

## `scripts/` (maintainer-only, not needed for ordinary use)

* `benchmark_cpu_time.py` — empirical Gaussian-vs-classifier CPU time comparison, feeds a manuscript
  figure.
* `compare_rerun.py` — G09-vs-G16 rerun consistency check (historical; the G16 promotion this validated
  is complete, see below).
* `compare_weighting.py` — diffs `mu`-weighted vs `none`-weighted library scores.
* `regenerate_figures_at_tau_sb.py` — regenerates the `tau_SB`-sensitive figures/CSVs at an alternate
  threshold without touching canonical outputs.

## Archived/retired snapshots

`data/` and `ped/` carry a few intentionally-kept historical snapshots. None of these are read by the
live pipeline — they exist so a specific past state (a superseded basis set, a retired threshold, a
retired scoring scheme) stays reproducible without digging through git history.

* **`data/logs_g09/`, `data/gjf_g09/`, `data/fchk_g09/`, `data/characterised_modes_g09.csv`** — the
  original Gaussian 09 run, superseded when Gaussian 16 (Revision C.01) was promoted to canonical
  (2026-08-10). The canonical `data/logs/`, `data/gjf/`, `data/fchk/`, `data/characterised_modes.csv`
  hold the current G16 rerun.
* **`data/results_tau_sb_0.50/`** — the frozen pre-2026-08-20 output snapshot from when `tau_SB=0.50` was
  the canonical default (it is now `0.42`; see `Thresholds.tau_SB`'s default in `src/classifier.py`).
* **`data/results/archive/`, `data/results/archive_threeway/`, `data/results/archive_unweighted/`** (and
  the matching `data/figures/archive_threeway/`, `data/figures/archive_unweighted/`) — frozen snapshots
  from, respectively, the pre-2026-08-25 two-gate-purity classifier, the pre-binary-scheme (threeway)
  classification default, and the pre-2026-08-14 unweighted `V_Stretch` definition.
* **`ped/archive_python_ped/`** — a superseded from-scratch PED pipeline, retired in favor of the real-VEDA4
  bridge (`ped/build_veda_fmt.py`); see its own `ARCHIVE_NOTE.md` for the full account of why.

## Tests

```bash
pip install -r requirements-dev.txt
pytest
```

No `conftest.py`/`pytest.ini` — plain pytest auto-discovery over `tests/`. All tests should pass.
