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

Opt-in stages (run only with `--only`): `sync_manuscript`, `si_tables`, and `si_analyses`, the
Supporting-Information analyses of the JCC revision (about 1.5 min): `src/si_degenerate.py` (s[V_S]
over all bases of each degenerate set), `src/si_frame.py` (frame spread; on-axis cut-off),
`src/si_water_overlaps.py` (H2O scores vs Cartesian / mass-weighted overlaps) and
`src/si_tau_robustness.py` (robustness of tau_SB and the S/B ambiguity range). They write
`data/results/si_*.csv` and `data/figures/fig_si_*`; each module also runs on its own
(`python -m src.si_degenerate`, ...).

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

## Principal-axis frame (`src/frame.py`)

The T/R scores are evaluated in a frame chosen from the input geometry and the input modes
(`main.build_scorer_and_final`, the single place every scorer is built):

* **Symmetric and spherical tops** — step 1, symmetry (Gaussian's convention): `z` = the `C_n` axis, `y`
  through an atom of the smallest set of symmetry-equivalent off-axis atoms, `x = y × z`; spherical tops:
  every set of three perpendicular `C2` axes (`C4` for O/Oh). Step 2: among these allowed frames (grouped
  by axis lines), the largest `Σ_k (D_T,k + D_R,k)` over the input modes, `D` = largest minus second-largest
  `|s|` over the rotatable axes. Step 3: names and signs (`+z` to the heaviest on-axis atom off the centre,
  else `Σ m z³ > 0`; lowest-index atom in the `yz` plane on `+y`; then the largest coordinates in input
  order). Step 3 never changes a score magnitude, a label or the T/R assignment.
* **Asymmetric tops** — eigenvectors of the inertia tensor in increasing moment order, right-handed, signs
  by the step-3 key.
* **Linear molecules** — `x` = molecular axis towards the heaviest atom off the centre; the turn about `x`
  with the largest `ΣD` from a 0.1° rescoring scan; `y`/`z` named by the first mode with an in-plane score.
* EMIT vectors are stored in `MIT()`'s frame and are returned to the input frame
  (`ModeScorer.mit_rotation_matrix`) before the frame is chosen; `--emit-projection` asserts that the
  normal-mode and EMIT pools end in the same frame.

Flags:
* `--frame-rule {symmetry,eigh}` — `symmetry` (default) is the rule above; `eigh` restores the submitted
  version's `MIT()` frame (eigenvectors with the heaviest-atom sign fix; can be left-handed).
* `--keep-input-frame` — keep the input frame (e.g. Gaussian's standard orientation) when it is one of the
  allowed frames. Without it the rule alone reproduces Gaussian's axis lines for 24 of the 40 roster tops
  with 3N−6 input (all 40 with 3N input); the rest get a symmetry-equivalent frame.
* `--sym-tol TOL` — symmetry tolerance in Å (default `1.5e-5`).

Known limits:
* `SYM_TOL` is calibrated on the roster (tops symmetric to ≤ 4.5e-6 Å; the ten C2v near-tops miss a `C_n`
  axis by ≥ 5e-5 Å and stay asymmetric). Coordinates printed more coarsely than 1e-5 Å need a larger value.
* Accidentally equal moments without a `C_n` (n ≥ 3) axis are not covered (the eigenvector frame is then arbitrary).
* Linear molecules: the turn is fixed to the 0.1° step (re-oriented input: axes within 0.06°, per-axis
  scores within 4e-4). With a central atom, input bent by more than the on-axis cut-off (`ON_AXIS_TOL`)
  makes `s[R_y]`, `s[R_z]` jump with the turn.
* Axis names and signs can depend on the atom numbering where the alternatives are symmetry-equivalent;
  never on the input orientation, and never in `|s|`, labels or the assignment.
* The on-axis cut-off must exceed the geometry's symmetry error (B3N3H6: 1.5e-6 Å < 1e-5 Å).

## Per-molecule flags (`-m`/`--molecule` required; each needs its prerequisite CSV to already exist)

Project raw EMIT eigenvectors onto the mass-weighted normal-mode reference basis. Requires
`data/results/<mol>_EMIT.csv` to already exist (run `--mode emit` first); merges the grouped
`C2_Tx`..`C2_VMix` fractions into that file **in place** and writes `data/results/<mol>_EMIT_full.csv`
(per-reference-mode detail). Fails loudly with a clear message if a prerequisite is missing.

```bash
python main.py -m H2O --mode emit
python main.py -m H2O --emit-projection
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

Build the library dataset (every molecule in `data/mol_list_method.csv`) from `data/logs/`+`data/gjf/`
into `data/results/library_scores.csv`.
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
`data/results/combined_ped_vs_scores.csv`, override with `--combined-output`). Note: this roster is every
molecule in `data/mol_list_method.csv` (the calibration set and the `test` set alike) — use
`ped/merge_ped_scores.py --molecules ...` directly (see `ped/README.md`) to merge PED for an explicit
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
  the canonical default (it is now `0.43`, the midpoint of the single-centre minimum-error plateau rounded to 0.01; see `Thresholds.tau_SB`'s default in `src/classifier.py`).
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
