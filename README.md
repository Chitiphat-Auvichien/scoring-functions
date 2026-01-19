# Scoring Functions for Classifying Modes of Molecular Motion

**Authors:** Chitiphat Auvichien, Nathani Therdpraisan, Phitchanaka Lertmankha, and Nattapong Paiboonvorachat

**Affiliation:** Department of Chemistry, Faculty of Science, Chulalongkorn University

This repository contains the Python implementation of the scoring functions described in our paper submitted to the *Journal of Chemical Education*.

## 1. Overview

Standard visualization of molecular vibrations can be subjective. This tool uses vector mathematics to quantitatively classify normal modes into three categories:

1. **Translational Scores (T):** Motion along principal axes (x, y, z).

2. **Rotational Scores (R):** Rotation around principal axes (x, y, z).

3. **Vibrational Scores (V):** Distinguishes **Stretching** (high score) from **Bending** (low score).

## 2. Installation

1. Ensure you have **Python 3** or higher installed.

2. Install the required dependencies using the provided file:

   ```
   pip install -r requirements.txt
   ```
## 3. Directory Structure
The program expects files to be organized systematically in the data/ folder:

* `data/logs/`: **[INPUT]** Place your Gaussian `.log` or `.out` files here.

* `data/EMIT/`: **[INPUT]** (Optional) Place EMIT mode text files here (e.g., benzene_EMIT.txt).

* `data/gjf/`: **[INPUT]** (Optional) Place Gaussian input files (.com or .gjf) here for automatic connectivity reading.

* `data/intermediate/`: **[EDITABLE]** The program generates readable text files here for you to verify/add bonds.

* `data/results/`: **[OUTPUT]** Final CSV score tables are saved here.

**File Naming Convention**: For a molecule named benzene, name your files consistently:

* Log: `data/logs/benzene.log`

* EMIT: `data/EMIT/benzene_EMIT.txt`

* GJF: `data/gjf/benzene.com`

## 4. Usage Guide
### Step 1: Prepare your files
Copy your Gaussian log file (and optionally EMIT/GJF files) into the correct subfolders inside `data/`.

### Step 2: Run the Script
Run the main program by specifying the molecule name (without extension) using the `-m` flag:

```
python main.py -m benzene
```

### Step 3: Choose Calculation Type
The program will ask:

1. **Normal Modes:** Calculates scores for the standard normal modes found in the Gaussian log.

2. **EMIT Modes:** Calculates scores for the 3*N* EMIT modes found in the corresponding text file in `data/EMIT/`.

### Step 4: Check Bonds
The program will generate an intermediate file (e.g., `data/intermediate/benzene_data.txt`).

* **If bonds are found (from GJF or Log):** It proceeds automatically.

* **If NO bonds are found:** It will **PAUSE** and ask you to open the intermediate text file to manually add bonds (e.g., 1 2 for a bond between Atom 1 and Atom 2 or copy directly from the GJF file).

* Save the file and press ENTER to continue.

### Step 5: Get Results
Results are saved to `data/results/benzene_normal_scores.csv` (or `_EMIT_scores.csv`).

## 5. Citation

> Auvichien, C.; et al. Scoring functions for classifying modes of molecular motion: I. Bridging mathematics and chemistry education" (Citation will be updated after submission to the *J. Chem. Educ.*)
