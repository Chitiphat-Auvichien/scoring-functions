# Scoring Functions for Classifying Modes of Molecular Motion

**Authors:** Chitiphat Auvichien, Nathani Therdpraisan, Phitchanaka Lertmankha, and Nattapong Paiboonvorachat

**Affiliation:** Department of Chemistry, Faculty of Science, Chulalongkorn University

This repository contains the Python implementation of the scoring functions described in our paper submitted to the *Journal of Chemical Education*.

## 1. Overview

Standard visualization of molecular vibrations can be subjective. This tool uses vector mathematics to quantitatively classify normal modes into three categories:

* **Translational Scores (T):** Does the whole molecule move along the X, Y, or Z axis?
* **Rotational Scores (R):** Does the whole molecule rotate around the X, Y, or Z axis?
* **Vibrational Scores (V):** Is the motion a bond **Stretching** (high score) or **Bending** (low score)?

## 2. Getting Started

### Prerequisites
You need **Python 3.0** or higher installed on your computer.

### Installation
1.  Download this repository to your computer.
2.  Open your terminal (Command Prompt, PowerShell, or Terminal).
3.  Navigate to the folder and install the required libraries:

    ```bash
    pip install -r requirements.txt
    ```
    
## 3. How to Use

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

### Step 4: Verify Bonds
The program generates a text file in `data/intermediate/` containing the molecule's geometry and the mode vectors.
* **If bonds are missing (e.g., there is no `.gjf` file):** The program will pause and ask you to open this text file.
* **Action:** Open the file, verify the `BONDS` section, and add any missing bonds (e.g., `1 2` for a bond between Atom 1 and Atom 2 or you can copy the connectivity information directly from the Gaussian input file). Save and close the file, then press **Enter** in the terminal.

### Step 5: View Results
The results are saved as a CSV file in `data/results/`. You can open this file in Excel.

**Note:** This repository is designed to facilitate calculations of scores for vibrational modes from the Gaussian program or for EMIT modes. However, the user can calculate scores for modes of motion obtained from any programs or methods by adapting the output format to suit this program.

## 5. Citation
If you use this code in your class or research, please cite:
> Auvichien, C.; Therdpraisan, N.; Lertmankha, P.; Paiboonvorachat, N. Scoring functions for classifying modes of molecular motion: I. Bridging mathematics and chemistry education" (Citation will be updated after submission to the *J. Chem. Educ.*)
