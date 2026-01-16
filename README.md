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

* `data/logs/`: **[INPUT]** Place your Gaussian `.log` files here.

* `data/intermediate/`: **[EDITABLE]** Extracted geometry files. **Edit this file to add bonds manually if needed.**

* `data/results/`: **[OUTPUT]** Final CSV score tables.

## 4. Usage Guide

1. Place your log file obtained from the Gaussian program (e.g., `h2o.log`) in the `data/logs/` folder.

2. Run the main script:

   ```
   python main.py data/logs/h2o.log
   ```

3. **Check Bonds:** The script will check for connectivity. If none is found, it will **PAUSE** and ask you to edit the generated text file in `data/intermediate/` to add bonds manually (e.g., `1 2` for a bond between Atom 1 and Atom 2). You can copy them from the Gaussian input file.

4. **Results:** Once confirmed, scores are saved to `data/results/h2o_scores.csv`.

## 5. Citation

> Auvichien, C.; et al. Scoring functions for classifying modes of molecular motion: I. Bridging mathematics and chemistry education", *J. Chem. Educ.* **2026**, ...
