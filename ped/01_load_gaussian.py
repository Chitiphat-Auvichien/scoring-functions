"""
Step 1: Load benzene geometry + all 30 vibrational normal modes directly
from the actual Gaussian log used in the JCC manuscript -- no PySCF, no
independent geometry optimization. This is the same MP2/3-21G, D6h,
freq=hpmodes calculation (data/logs/C6H6_MP2_3-21G_D6h.log) that produces
every benzene figure/table in the paper.

Requires: the parent scoring-functions repo's src/parser.py, src/utils.py
Output: opt_coords_ang.npy (12x3, Angstrom), opt_symbols.txt,
        l_raw.npy (30x12x3, Gaussian's raw printed Cartesian displacement
                   per mode -- normalization convention NOT assumed here,
                   verified empirically in step 2),
        vibfreq.npy (30,), reduced_mass.npy (30,), bonds.txt (0-indexed pairs)
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from src.parser import GaussianParser

LOG_PATH = os.path.join(os.path.dirname(__file__), '..', 'data', 'logs',
                         'C6H6_MP2_3-21G_D6h.log')

parser = GaussianParser(LOG_PATH)
data = parser.parse(parse_modes=True)

symbols = data["atoms"]
coords_ang = data["coords"]
modes = data["modes"]
bonds = data["bonds"]
natm = len(symbols)

n_expected = 3 * natm - 6
if len(modes) != n_expected:
    raise ValueError(f"Expected {n_expected} vibrational modes for benzene "
                      f"(3N-6, N={natm}), parser returned {len(modes)}.")

if any(m["reduced_mass"] is None for m in modes):
    raise ValueError("Some modes are missing a parsed reduced_mass -- "
                      "needed in step 2 to identify Gaussian's displacement "
                      "normalization convention.")

vibfreq = np.array([m["frequency"] for m in modes])
reduced_mass = np.array([m["reduced_mass"] for m in modes])
l_raw = np.array([m["vector"] for m in modes])  # (30, 12, 3)

np.save(os.path.join(os.path.dirname(__file__), 'opt_coords_ang.npy'), coords_ang)
with open(os.path.join(os.path.dirname(__file__), 'opt_symbols.txt'), 'w') as f:
    for s in symbols:
        f.write(s + '\n')
np.save(os.path.join(os.path.dirname(__file__), 'l_raw.npy'), l_raw)
np.save(os.path.join(os.path.dirname(__file__), 'vibfreq.npy'), vibfreq)
np.save(os.path.join(os.path.dirname(__file__), 'reduced_mass.npy'), reduced_mass)
with open(os.path.join(os.path.dirname(__file__), 'bonds.txt'), 'w') as f:
    for a, b in bonds:
        f.write(f"{a} {b}\n")

print(f"Loaded {LOG_PATH}")
print(f"Atoms ({natm}): {' '.join(symbols)}")
print(f"Bonds ({len(bonds)}): {bonds}")
print(f"Vibrational modes: {len(modes)} (expected 3N-6={n_expected})")
print("\nFrequencies (cm^-1):")
print(np.round(np.sort(vibfreq), 2))
print("\nSaved: opt_coords_ang.npy, opt_symbols.txt, l_raw.npy, "
      "vibfreq.npy, reduced_mass.npy, bonds.txt")
