"""
Step 9 (VEDA-style cross-check): final PED table from the optimized/reduced
natural coordinates (B_final.npy), grouped into the SAME five categories as
05_final_table.py's redundant-coordinate table, for direct comparability.

Category assignment: the natural-coordinate construction (step 06) found
NO independent CCC_bend direction for benzene's ring at all -- every
candidate CCC_bend combination turned out to be a linear combination of
CC_stretch/ring_torsion-derived directions already kept (see
nat_coord_labels.txt and 06_natural_coords.py's printed KEEP/DROP list).
This is a real, honestly-derived numerical result, not a bug: it means
this coordinate system has no basis vector whose *dominant* physical
character is "C-C-C angle bending independent of everything else" -- so
the CCC-bend column in the table below will legitimately be ~0% for every
mode. See VEDA_STYLE_METHODOLOGY.md for the full discussion.

Each of the 30 final coordinates (rows of M_final, relative to the natural
basis in nat_coord_labels.txt) is assigned to whichever REPORTED category
carries the largest summed squared weight in its composition -- explicit,
documented convention, not VEDA's own (unpublished) rule.

Requires: B_final.npy (step 08); nat_coord_labels.txt (step 06);
          H_cart.npy, L.npy, vibfreq.npy (step 02)
Output: veda_style_PED_table.csv, veda_style_PED_table.txt (same schema as
        benzene_PED_table.csv/.txt)
"""
import csv
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from _ped_core import compute_PED, freq_to_lambda, EPm

HERE = os.path.dirname(__file__)

# family -> reported category (mirrors 05_final_table.py's merging of
# ring_torsion+CH_wag_oop into "out-of-plane", and CCH_bend_sym+anti into
# "CCH bend")
FAMILY_TO_CATEGORY = {
    'CH_stretch': 'CH stretch',
    'CC_stretch': 'CC stretch',
    'CCC_bend': 'CCC bend',
    'CCH_bend_sym': 'CCH bend',
    'CCH_bend_anti': 'CCH bend',
    'ring_torsion': 'out-of-plane',
    'CH_wag_oop': 'out-of-plane',
}
CATEGORIES = ['CH stretch', 'CC stretch', 'CCC bend', 'CCH bend', 'out-of-plane']

B_final = np.load(os.path.join(HERE, 'B_final.npy'))
M_final = np.load(os.path.join(HERE, 'M_final.npy'))
H = np.load(os.path.join(HERE, 'H_cart.npy'))
L = np.load(os.path.join(HERE, 'L.npy'))
vib = np.load(os.path.join(HERE, 'vibfreq.npy'))
lam = freq_to_lambda(vib)

nat_labels, nat_families = [], []
for line in open(os.path.join(HERE, 'nat_coord_labels.txt')):
    label, family, species = line.rstrip('\n').split('\t')
    nat_labels.append(label)
    nat_families.append(family)

# --- assign each of the 30 FINAL coordinates to a reported category ---
final_category = []
for k in range(30):
    weight = {c: 0.0 for c in CATEGORIES}
    for comp, coeff in enumerate(M_final[k]):
        if abs(coeff) < 1e-9:
            continue
        cat = FAMILY_TO_CATEGORY[nat_families[comp]]
        weight[cat] += coeff ** 2
    final_category.append(max(weight, key=weight.get))

print("Final-coordinate category assignment (30 coordinates):")
for cat in CATEGORIES:
    n = final_category.count(cat)
    print(f"  {cat:12s}: {n}")

PED = compute_PED(B_final, H, L, lam)   # (30, 30): coordinate x mode
colsums = PED.sum(axis=0)
print(f"\nPED column sums: min={colsums.min():.6f}, max={colsums.max():.6f}")
print(f"EPm (this table): {EPm(PED):.6f}")

# Sum raw (signed) PED within each category first, then normalize -- same
# convention as 04_compute_ped.py/05_final_table.py (no abs() on individual
# components before grouping, so a small negative category percentage from
# a non-orthogonal-adjacent coordinate is preserved and visible, not masked).
cat_pct = {c: np.zeros(30) for c in CATEGORIES}
for k, cat in enumerate(final_category):
    cat_pct[cat] += PED[k]
for c in CATEGORIES:
    cat_pct[c] = 100 * cat_pct[c] / colsums

rows = sorted(zip(vib, cat_pct['CH stretch'], cat_pct['CC stretch'],
                   cat_pct['CCC bend'], cat_pct['CCH bend'], cat_pct['out-of-plane']),
              key=lambda r: r[0])

grouped = []
i = 0
while i < len(rows):
    if i + 1 < len(rows) and abs(rows[i][0] - rows[i + 1][0]) < 0.5:
        avg = tuple((a + b) / 2 for a, b in zip(rows[i], rows[i + 1]))
        grouped.append((avg, True))
        i += 2
    else:
        grouped.append((rows[i], False))
        i += 1

hdr = (f"{'#':>3} {'freq(cm-1)':>10} {'deg':>4}  {'CHstr':>6} {'CCstr':>6} "
       f"{'CCCbend':>8} {'CCHbend':>8} {'out-of-plane':>13}  dominant character")
lines = [hdr, "-" * len(hdr)]
csv_rows = []
for n, (row, isdeg) in enumerate(grouped, 1):
    freq, ch, cc, ccc, cch, oop = row
    vals = {'CH stretch': ch, 'CC stretch': cc, 'CCC bend': ccc,
            'CCH bend': cch, 'out-of-plane': oop}
    dom = sorted(vals.items(), key=lambda x: -x[1])
    domstr = " + ".join(f"{v:.0f}% {k}" for k, v in dom if v > 5)
    degstr = "(x2)" if isdeg else ""
    lines.append(f"{n:3d} {freq:10.1f} {degstr:>4}  {ch:6.1f} {cc:6.1f} "
                  f"{ccc:8.1f} {cch:8.1f} {oop:13.1f}  {domstr}")
    csv_rows.append({
        'mode': n, 'freq_cm-1': round(freq, 1), 'degenerate': isdeg,
        'CH_stretch_pct': round(ch, 1), 'CC_stretch_pct': round(cc, 1),
        'CCC_bend_pct': round(ccc, 1), 'CCH_bend_pct': round(cch, 1),
        'out_of_plane_pct': round(oop, 1), 'dominant_character': domstr,
    })

footer = [
    f"\n{len(grouped)} unique lines representing all 30 vibrational modes "
    f"(degenerate e-type pairs shown once, marked (x2)).",
    "\nFrequencies are Gaussian's own computed harmonic values (unscaled), "
    "from data/logs/C6H6.log (MP2/3-21G, D6h, freq=hpmodes) "
    "-- the same benzene calculation used throughout the JCC manuscript.",
    "\nCoordinates: non-redundant natural-coordinate basis (steps 06-08), "
    "NOT the redundant 42-coordinate basis of benzene_PED_table.txt.",
    "\nCCC-bend is legitimately ~0% throughout: the natural-coordinate "
    "construction found no independent CCC-bend direction for this ring "
    "(see VEDA_STYLE_METHODOLOGY.md) -- not a computation error.",
]

for line in lines + footer:
    print(line)

txt_path = os.path.join(HERE, 'veda_style_PED_table.txt')
with open(txt_path, 'w') as f:
    f.write("\n".join(lines + footer) + "\n")

csv_path = os.path.join(HERE, 'veda_style_PED_table.csv')
with open(csv_path, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(csv_rows[0].keys()))
    writer.writeheader()
    writer.writerows(csv_rows)

print(f"\nSaved: {txt_path}")
print(f"Saved: {csv_path}")
