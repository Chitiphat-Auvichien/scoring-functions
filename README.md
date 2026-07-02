# Scoring Functions for Classifying Modes of Molecular Motion
![Figure abstract](https://github.com/Chitiphat-Auvichien/scoring-functions/blob/main/Abstract.jpg)

## Overview

Standard visualization of molecular vibrations can be subjective. This repository implements a
**reference-free, vector-mathematical framework** that scores every one of a molecule's `3N` modes
of motion for translational, rotational, and vibrational character, and then classifies each mode
into one of six categories — no pre-computed reference/"clean" mode set required.

This is a two-paper research codebase, and both papers are active:

* **Paper I** (the original scoring-functions concept — three raw per-axis scores) was submitted to
  the *Journal of Chemical Education*. It is the historical starting point of this repository and its
  citation entry below still reflects "submitted" status.
* **Paper II**, *"A Unified, Reference-Free Framework for Classifying the 3N Modes of Molecular
  Motion,"* is in preparation for the *Journal of Computational Chemistry* (JCC) and substantially
  extends Paper I's scoring functions into a full classification algorithm, an EMIT-mode projection
  method, and a calibrated, validated pipeline. It does not invalidate or supersede Paper I; it builds
  on it.

Concretely, the code now covers:

* **Step 1 — scoring.** Every mode gets six per-axis external scores plus a stretch/bend score:
  * **Translational scores (Tx, Ty, Tz):** does the whole molecule move along the X/Y/Z axis?
  * **Rotational scores (Rx, Ry, Rz):** does the whole molecule rotate about the X/Y/Z axis?
  * **Vibrational score (`V_Stretch`):** is the internal motion bond **stretching** (high score) or
    **bending** (low score)?
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
  against a ~25-molecule hydride-library dataset (ingested from a spreadsheet of precomputed scores,
  never re-derived from scratch) and validated with a confusion matrix and a threshold-sensitivity
  sweep.
* **Figures.** All of the manuscript's data-driven figures can be regenerated from the pipeline's own
  CSV outputs.

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
    projection, and figure pipelines. The library-ingestion pipeline (`src/excel_ingest.py`) additionally
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

### Step 4: Verify Bonds
The program generates a text file in `data/intermediate/` containing the molecule's geometry and the mode vectors.
* **If bonds are missing (e.g., there is no `.gjf` file):** The program will pause and ask you to open this text file.
* **Action:** Open the file, verify the `BONDS` section, and add any missing bonds (e.g., `1 2` for a bond between Atom 1 and Atom 2 or you can copy the connectivity information directly from the Gaussian input file). Save and close the file, then press **Enter** in the terminal.

### Step 5: View Results
The results are saved as a CSV file in `data/results/` (`<molecule>_normal_scores.csv` or
`<molecule>_EMIT_scores.csv`). You can open this file in Excel.

**Note:** This repository is designed to facilitate calculations of scores for vibrational modes from the Gaussian program or for EMIT modes. However, the user can calculate scores for modes of motion obtained from any programs or methods by adapting the output format to suit this program.

## Classification, EMIT Projection, Calibration, and Figures

`main.py`'s command line currently only exposes the Step-1 scoring workflow above (`-m`/`--molecule`
and `--mode`). The classification, EMIT-projection, library-calibration, and figure-generation
pipelines are fully built and tested, but — as of this writing — are reachable as **Python functions**
you import and call yourself, not as `main.py` subcommands. (Wiring them into `main.py` as
`--classify`/`--emit-projection`/`--library`/`--figures` flags, and a `reproduce.py` orchestrator that
chains everything headlessly, is tracked as outstanding work in `IMPLEMENTATION_PLAN.md`.)

Run these from the repository root (`Github/scoring-functions/`):

```python
from main import run_classify_pipeline, run_projection_pipeline
from src.excel_ingest import run_ingest_pipeline
from src.calibrate import run_calibration_pipeline

# Step 1 scores + Steps 2-4 classification labels for every mode of a molecule.
# mode_type is "normal" or "emit"; writes data/results/<mol>_<normal|EMIT>_classified.csv
df, path = run_classify_pipeline("benzene", "normal")

# Project raw EMIT eigenvectors onto the mass-weighted normal-mode reference basis.
# Writes data/results/<mol>_EMIT_contributions.csv (grouped) and
# data/results/<mol>_EMIT_projection_full.csv (per-reference-mode detail).
df_grouped, df_full, (path_grouped, path_full) = run_projection_pipeline("benzene")

# Ingest the hydride-library spreadsheet and attach geometry-backed classifications.
# Writes data/results/library_scores.csv.
df_lib, path, skip_report = run_ingest_pipeline()

# Calibrate tau_TR/tau_S/tau_B against the ingested library and run the
# threshold-sensitivity sweep. Writes data/results/thresholds.json and
# data/results/tau_sensitivity_sweep.csv.
thresholds, result, sweep_df, (path_json, path_sweep) = run_calibration_pipeline()
```

Every manuscript figure is a function in `src/figures.py` (e.g. `plot_benzene_stress_test`,
`plot_confusion_matrix`, `plot_bond_scores`, `plot_boxplots`, `plot_mode_mixing`,
`plot_sensitivity`, `plot_benzene_normal_modes`); each reads already-computed
`data/results/*.csv`/`thresholds.json` and writes a vector PDF + PNG pair to `data/figures/`. Running
the module directly regenerates all of them in one go:

```bash
python -m src.figures
```

## Citation
If you use this code in your class or research, please cite:

> Auvichien, C.; Therdpraisan, N.; Lertmankha, P.; Paiboonvorachat, N. "Scoring functions for
> classifying modes of molecular motion: I. Bridging mathematics and chemistry education" (Citation
> will be updated after submission to the *J. Chem. Educ.*)

The classification framework, EMIT-projection method, and calibration/validation pipeline described
above are the subject of a second, separate manuscript, in preparation for the *Journal of
Computational Chemistry*: *"A Unified, Reference-Free Framework for Classifying the 3N Modes of
Molecular Motion."* A full citation will be added here once that manuscript is submitted.
