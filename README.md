# Scoring Functions for Classifying Modes of Molecular Motion
![Figure abstract](https://github.com/Chitiphat-Auvichien/scoring-functions/blob/main/Abstract.jpg)

## Overview

Standard visualization of molecular vibrations can be subjective. This repository implements a
**displacement-based, vector-mathematical framework** that scores every one of a molecule's `3N` modes
of motion for translational, rotational, and vibrational character, directly from atomic displacement
vectors, atom connectivity, and molecular geometry — no Hessian and no user-supplied reference
coordinates required.

This is the reference implementation for the manuscript *"A Unified, Displacement-Based Framework for
Classifying Modes of Molecular Motion,"* by Chitiphat Auvichien, Nathani Therdpraisan, Phitchanaka
Lertmankha, Thanchon Boonkrong, Viwat Vchirawongkwin, and Nattapong Paiboonvorachat, submitted to the
*Journal of Computational Chemistry* (JCC).

Concretely, the code covers:

* **Step 1 — scoring.** Every mode gets seven characteristic scores:
  * **Translational scores (Tx, Ty, Tz):** does the whole molecule move along the X/Y/Z axis?
  * **Rotational scores (Rx, Ry, Rz):** does the whole molecule rotate about the X/Y/Z axis?
  * **Vibrational score (`V_Stretch`):** is the internal motion bond **stretching** (high score) or
    **bending** (low score)? Each bond is weighted by its reduced mass `μ_AB = m_A·m_B/(m_A+m_B)`, so a
    bond counts for as much as the kinetic energy its stretching motion actually carries. Since `μ`
    divides out when every bond is the same (an AB_n molecule like H₂O or CH₄), this only changes
    molecules with mixed bond types. `--v-weighting none` restores the original unweighted definition.
* **Step 2 — vibrational classification.** Every one of the `3N` candidate modes gets a `vib_label`
  from `V_Stretch` alone: `"S"` (stretching) or `"B"` (bending) under the default binary scheme, via a
  single calibrated cutoff `tau_SB` (default `0.42`). An explicit `--classify-scheme threeway` opt-in can
  instead label a mode `"SB"` (mixed stretch-bend) using a separate `tau_S`/`tau_B` split.
* **Step 3 — T/R identification (optional).** A Hungarian (`scipy.optimize.linear_sum_assignment`) global
  assignment finds the `n_T + n_R` modes that best represent translation/rotation, each getting a bare
  slot name (`Tx`, `Ty`, `Tz`, `Rx`, `Ry`, `Rz`) and its own raw signed `tr_score` — **unconditionally**,
  with no purity gate. This step is independent of Step 2 and mainly useful when the input is an
  arbitrary/non-normal-mode basis (e.g. EMIT modes) where translation/rotation aren't already known by
  construction. Pass `--no-identify-tr` to skip it.
* **EMIT-mode projection.** Raw EMIT eigenvectors can be projected onto a mass-weighted normal-mode
  reference basis (ideal T/R + real vibrational modes) to get fractional T/R/stretch/bend contributions
  per EMIT mode.
* **Library calibration.** `tau_SB` is fixed by author calibration against a hydride-library dataset
  built directly from the Gaussian files in `data/logs/`+`data/gjf/` (every score is computed by this
  program's own engine; a reference spreadsheet supplies only the literature stretch/bend ground-truth
  labels, never the scores themselves), validated with a confusion matrix and a threshold-sensitivity
  sweep. `tau_TR`/`tau_S`/`tau_B` are also calibrated but are diagnostic-only — they no longer gate any
  classification decision.
* **Figures.** All of the manuscript's data-driven figures can be regenerated from the pipeline's own
  CSV outputs.
* **Real-VEDA4 PED cross-check.** `ped/` bridges a molecule's Gaussian Hessian into a real VEDA4
  Potential Energy Distribution (PED) analysis, and merges VEDA4's ground-truth %stretch/%bend character
  back into this program's own scores for a genuine, non-reimplemented validation of `V_Stretch`. See
  `ped/README.md` for the full workflow.

## Getting Started

### Prerequisites
You need **Python 3.9** or higher installed on your computer.

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

### Step 4: Verify Bonds
The program generates a text file in `data/intermediate/` containing the molecule's geometry and the mode vectors.
* **If bonds are missing (e.g., there is no `.gjf` file):** The program will pause and ask you to open this text file.
* **Action:** Open the file, verify the `BONDS` section, and add any missing bonds (e.g., `1 2` for a bond between Atom 1 and Atom 2 or you can copy the connectivity information directly from the Gaussian input file). Save and close the file, then press **Enter** in the terminal.

### Step 5: View Results
The results are saved as one CSV file per mode type in `data/results/` (`<molecule>_normal.csv` or
`<molecule>_EMIT.csv`) — the seven scores, `Mu`/`K`/`Irrep`, and the Steps 2-3 classification
(`vib_label`, `tr_label`, `tr_score`) are all in that single file. You can open this file in Excel.

**Note:** This repository is designed to facilitate calculations of scores for vibrational modes from the Gaussian program or for EMIT modes. However, the user can calculate scores for modes of motion obtained from any programs or methods by adapting the output format to suit this program.

### The web app

The classification framework also has a browser-based front-end, *Vscore*, developed in this work and
available at <https://vscore.vercel.app/>. Its source lives in its own companion repository,
<https://github.com/ThanchonBK/vscore-webapp>, which vendors `src/{scoring,utils,classifier}.py` and
`data/results/thresholds.json` as verbatim copies — a change to the classification scheme here needs
re-syncing there. See that repo's own README for setup and its test suite, which compares every score
and label against this repository's `data/results/*.csv` when checked out alongside it.

## Repository layout

* `src/` — the core package (scoring, classification, projection, calibration, parsing, figures).
* `main.py` / `reproduce.py` — the CLI entry point and the full-pipeline orchestrator.
* `data/` — inputs (`gjf/`, `logs/`, `EMIT/`) and outputs (`results/`, `figures/`, `intermediate/`); see
  `docs/DEVELOPMENT.md` for what the archived/retired snapshot folders inside `data/` are.
* `ped/` — the real-VEDA4 PED cross-check bridge.
* `tests/` — the pytest suite.
* `docs/DEVELOPMENT.md` — calibration, EMIT projection, the PED cross-check, figure regeneration, and
  every other maintainer-facing CLI flag and script not needed for ordinary use.

## License

MIT — see [`LICENSE.md`](LICENSE.md).

## Citation
If you use this code in your class or research, please cite:

> A full citation for *"A Unified, Displacement-Based Framework for Classifying Modes of Molecular
> Motion"* by C. Auvichien, N. Therdpraisan, P. Lertmankha, T. Boonkrong, V. Vchirawongkwin, and
> N. Paiboonvorachat (submitted to the *Journal of Computational Chemistry*) will be added here once the
> manuscript is accepted.
