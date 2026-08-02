#!/bin/bash
# NOTE (archived 2026-08-02): this pipeline was moved into archive_python_ped/
# -- run this script from INSIDE archive_python_ped/, not from ped/ directly.
set -e
echo "=== Step 1: load Gaussian geometry + normal modes ==="
python3 01_load_gaussian.py
echo
echo "=== Step 2: reconstruct Hessian ==="
python3 02_reconstruct_hessian.py
echo
echo "=== Step 3: internal coordinates / B-matrix ==="
python3 03_build_internal_coords.py
echo
echo "=== Step 4: PED calculation ==="
python3 04_compute_ped.py
echo
echo "=== Step 5: final table ==="
python3 05_final_table.py
