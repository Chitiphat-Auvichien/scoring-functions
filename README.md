# Scoring Functions for Classifying Modes of Molecular Motion
![Figure abstract](https://github.com/Chitiphat-Auvichien/scoring-functions/blob/main/Abstract.jpg)

## Overview

Standard visualization of molecular vibrations can be subjective. This repository implements a
**reference-free, vector-mathematical framework** that scores every one of a molecule's `3N` modes
of motion for translational, rotational, and vibrational character, and then classifies each mode
into one of six categories — no pre-computed reference/"clean" mode set required.

This is the reference implementation for the manuscript *"A Unified, Reference-Free Framework for
Classifying the 3N Modes of Molecular Motion,"* in preparation for the *Journal of Computational
Chemistry* (JCC).

Concretely, the code covers:

* **Step 1 — scoring.** Every mode gets six per-axis external scores plus a stretch/bend score:
  * **Translational scores (Tx, Ty, Tz):** does the whole molecule move along the X/Y/Z axis?
  * **Rotational scores (Rx, Ry, Rz):** does the whole molecule rotate about the X/Y/Z axis?
  * **Vibrational score (`V_Stretch`):** is the internal motion bond **stretching** (high score) or
    **bending** (low score)? Each bond is weighted by its reduced mass `μ_AB = m_A·m_B/(m_A+m_B)`, so a
    bond counts for as much as the kinetic energy its stretching motion actually carries. Since `μ`
    divides out when every bond is the same (an AB_n molecule like H₂O or CH₄), this only changes
    molecules with mixed bond types. `--v-weighting none` restores the original unweighted definition.
* **Steps 2-4 — classification.** A Hungarian (`scipy.optimize.linear_sum_assignment`) global
  assignment of modes to the `n_T + n_R` external (translation/rotation) slots, a two-gate purity test,
  and a stretch/bend/mixed split for every remaining internal mode. Each mode ends up labeled as one
  of: a **clean** external slot (`Tx`, `Ty`, `Tz`, `Rx`, `Ry`, `Rz`), a **flagged mixed** external+
  vibration slot (same names with a trailing `*`, e.g. `Tx*`), **stretching** (`S`), **bending** (`B`),
  or **mixed stretch-bend** (`SB`). The framework classifies internal/external mixing; it does not
  attempt to quantify it — the `*` flag is the output, by design.
* **EMIT-mode projection.** Raw EMIT eigenvectors can be projected onto a mass-weighted normal-mode
  reference basis (ideal T/R + real vibrational modes) to get fractional T/R/stretch/bend
  contributions per EMIT mode.
* **Library calibration.** The two-gate purity thresholds (`tau_TR`, `tau_S`, `tau_B`) are calibrated
  against a hydride-library dataset built directly from the Gaussian files in `data/logs/`+`data/gjf/`
  (every score is computed by this program's own engine; a reference spreadsheet supplies only the
  literature stretch/bend ground-truth labels, never the scores themselves) and validated with a
  confusion matrix and a threshold-sensitivity sweep.
* **Figures.** All of the manuscript's data-driven figures can be regenerated from the pipeline's own
  CSV outputs.
* **Real-VEDA4 PED cross-check.** `ped/` bridges a molecule's Gaussian Hessian into a real VEDA4
  Potential Energy Distribution (PED) analysis (`ped/build_veda_fmt.py`), and `ped/merge_ped_scores.py`
  (also reachable via `main.py --ped-merge`/`--ped-merge-all`) merges VEDA4's ground-truth %stretch/%bend
  character back into this program's own scores for a genuine, non-reimplemented validation of
  `V_Stretch`. See `ped/README.md` for the full workflow.

## Getting Started

### Prerequisites
You need **Python 3.0** or higher installed on your computer.

### Installation
1.  Download this repository to your computer.
2.  Open your terminal (Command Prompt, PowerShell, or Terminal).
3.  Navigate to the folder and install the required libraries:

    ```bash
    pip install -r requirements.txt
    ```
    This installs `numpy`, `pandas`, `scipy`, and `matplotlib`, which cover the scoring, classification,
    projection, and figure pipelines. The library-ingestion pipeline (`--library` below) additionally
    needs `openpyxl` to read the source `.xlsx` workbook — install it separately (`pip install openpyxl`)
    if you plan to run that step.

## How to Use

The program is designed to be simple. You provide a molecule name, and it looks for files in the `data/` folder.

### Step 1: Prepare Your Files
1.  Run optimization + frequency calculations in Gaussian (e.g., `freq=hpmodes`).
2.  Rename your files to the molecule name (e.g., `benzene.com` and `benzene.log`). The name must be **consistent**.
3.  Place the Gaussian input file (`.com` or `.gjf`) in the **`data/gjf/`** folder and the Gaussian output file (`.log`) in the **`data/logs/`** folder.

*(Optional)* If you have EMIT mode files, place them in `data/EMIT/` named `benzene_EMIT.txt`.

### Step 2: Run the Script
In your terminal, run the following command (replace `benzene` with your molecule name):

```bash
python main.py -m benzene
```
### Step 3: Choose Calculation Mode
The program will ask you to choose a mode:
* **Type** `1`: To calculate scores for **Normal Modes** (standard Gaussian output).
* **Type** `2`: To calculate scores for **EMIT Modes** (advanced user option).

You can skip the interactive prompt with `--mode {normal,emit}`, e.g. `python main.py -m benzene --mode normal`.

### Reproducing everything

`python reproduce.py` rebuilds every `data/results/*.csv` and `data/figures/*` from the inputs in one
ordered pass (`--dry-run` lists the stages; `--only`/`--skip` select them). It writes
`data/results/_run_manifest.json` recording which scoring definition, commit and stages produced the
current outputs.

To rebuild under the original unweighted V-score instead, use `python reproduce.py --v-weighting none`.
Note that scoring and thresholds must come from the same definition: `thresholds.json` records which one
it was calibrated for, and the program refuses to classify with a mismatched pair rather than silently
mislabelling modes. The pre-2026-08-14 unweighted outputs are kept under
`data/results/archive_unweighted/` and `data/figures/archive_unweighted/`;
`python scripts/compare_weighting.py` diffs the two.

### Step 4: Verify Bonds
The program generates a text file in `data/intermediate/` containing the molecule's geometry and the mode vectors.
* **If bonds are missing (e.g., there is no `.gjf` file):** The program will pause and ask you to open this text file.
* **Action:** Open the file, verify the `BONDS` section, and add any missing bonds (e.g., `1 2` for a bond between Atom 1 and Atom 2 or you can copy the connectivity information directly from the Gaussian input file). Save and close the file, then press **Enter** in the terminal.

### Step 5: View Results
The results are saved as one CSV file per mode type in `data/results/` (`<molecule>_normal.csv` or
`<molecule>_EMIT.csv`) — scores, `Mu`/`K`/`Irrep`, and the Steps 2-4 classification (`label`,
`annotation`, `s_AB`) are all in that single file; there is no separate scores-only output. You can open
this file in Excel.

**Note:** This repository is designed to facilitate calculations of scores for vibrational modes from the Gaussian program or for EMIT modes. However, the user can calculate scores for modes of motion obtained from any programs or methods by adapting the output format to suit this program.

## Developer / maintainer workflow (CLI flags)

Steps 1-5 above (`-m`/`--molecule` + `--mode`) are the whole job for ordinary use. `main.py` also exposes
a second group of flags — diagnostics, cross-checks, and manuscript-support pipelines — that most users
won't need. Run `python main.py --help` to see both groups. All commands below were run from the
repository root (`Github/scoring-functions/`) against the data already checked into this repo and
confirmed to produce real output files.

**Per-molecule flags** (`-m`/`--molecule` required; each requires its prerequisite CSV to already exist):

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

**Global flags** (molecule-independent; `-m` is ignored if supplied alongside them):

Build the hydride-library dataset from `data/logs/`+`data/gjf/` into `data/results/library_scores.csv`.
Grows automatically as more molecule files are added; can be slow once the library is large, since every
molecule is scored from scratch.

```bash
python main.py --library
```

Calibrate `tau_TR`/`tau_S`/`tau_B` against the library and run the threshold-sensitivity sweep. Writes
`data/results/thresholds.json` and `data/results/tau_sensitivity_sweep.csv`. Rebuilds the library itself
first, so this is also slow.

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
68-molecule hydride calibration library, not any transferability-test molecule set you may be tracking
separately — use `ped/merge_ped_scores.py --molecules ...` directly (see `ped/README.md`) to merge PED
for an explicit list of molecules instead.

```bash
python main.py --ped-merge-all
```

Every manuscript figure is also a standalone function in `src/figures.py` (e.g.
`plot_benzene_stress_test`, `plot_confusion_matrix`, `plot_bond_scores`, `plot_boxplots`,
`plot_mode_mixing`, `plot_sensitivity`, `plot_benzene_normal_modes`) if you need to regenerate just one
of them from a Python session; `python main.py --figures` (or `python -m src.figures`) regenerates all
of them in one go.

## Citation
If you use this code in your class or research, please cite:

> A full citation for *"A Unified, Reference-Free Framework for Classifying the 3N Modes of Molecular
> Motion"* (in preparation for the *Journal of Computational Chemistry*) will be added here once the
> manuscript is submitted.
