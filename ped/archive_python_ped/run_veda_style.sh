#!/bin/bash
# NOTE (archived 2026-08-02): this pipeline was moved into archive_python_ped/
# -- run this script from INSIDE archive_python_ped/ (it self-locates via
# BASH_SOURCE below, so it also works if invoked from elsewhere).
set -e
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

for f in H_cart.npy L.npy vibfreq.npy B.npy coord_labels.txt benzene_PED_table.csv; do
    if [ ! -f "$f" ]; then
        echo "ERROR: $f not found. Run 'bash run_all.sh' first (steps 01-05) --" >&2
        echo "the VEDA-style cross-check reuses their validated output and does not" >&2
        echo "regenerate it." >&2
        exit 1
    fi
done

echo "=== Step 6: natural (non-redundant) coordinates ==="
python3 06_natural_coords.py
echo
echo "=== Step 7: greedy EPm-maximizing mixing ==="
python3 07_veda_mixing.py
echo
echo "=== Step 8: reduction pass ==="
python3 08_veda_reduction.py
echo
echo "=== Step 9: VEDA-style PED table ==="
python3 09_veda_ped_table.py
echo
echo "=== Step 10: comparison report ==="
python3 10_compare_report.py
