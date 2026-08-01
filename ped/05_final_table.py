"""
Step 5: Merge the out-of-plane categories (ring torsion + C-H wag are
non-orthogonal for a planar hexagonal ring, so their individual PED
values can go negative/>100% -- a known, documented ambiguity, not a
bug), group numerically-degenerate frequency pairs, and print the
final table.

No frequency scale factor is applied: there is no verified published
scale factor for MP2/3-21G (unlike, e.g., Scott & Radom's well-known
0.8929 for HF/6-31G(d)), and PED percentages don't depend on scaling
anyway -- only the displayed cm^-1 column would change. Reported
frequencies are Gaussian's own computed (unscaled) harmonic values.

Requires: vibfreq.npy, PED_group_pct.npy, cats.txt
Output: benzene_PED_table.csv, benzene_PED_table.txt
"""
import csv
import os
import numpy as np

HERE = os.path.dirname(__file__)

PED = np.load(os.path.join(HERE, 'PED_group_pct.npy'))
cats = [l.strip() for l in open(os.path.join(HERE, 'cats.txt'))]
vib = np.load(os.path.join(HERE, 'vibfreq.npy'))
idx = {c: i for i, c in enumerate(cats)}

CC = PED[idx['CC_stretch']]
CH = PED[idx['CH_stretch']]
CCC = PED[idx['CCC_bend']]
CCH = PED[idx['CCH_bend']]
OOP = PED[idx['ring_torsion']] + PED[idx['CH_wag_oop']]

rows = sorted(zip(vib, CC, CH, CCC, CCH, OOP), key=lambda r: r[0])

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

hdr = f"{'#':>3} {'freq(cm-1)':>10} {'deg':>4}  {'CHstr':>6} {'CCstr':>6} {'CCCbend':>8} {'CCHbend':>8} {'out-of-plane':>13}  dominant character"
lines = [hdr, "-" * len(hdr)]
csv_rows = []
for n, (row, isdeg) in enumerate(grouped, 1):
    freq, cc, ch, ccc, cch, oop = row
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
    "\nOut-of-plane category = ring_torsion + CH_wag_oop merged (non-orthogonal "
    "internal coordinates for a planar ring -- see README.md).",
]

for line in lines + footer:
    print(line)

txt_path = os.path.join(HERE, 'benzene_PED_table.txt')
with open(txt_path, 'w') as f:
    f.write("\n".join(lines + footer) + "\n")

csv_path = os.path.join(HERE, 'benzene_PED_table.csv')
with open(csv_path, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(csv_rows[0].keys()))
    writer.writeheader()
    writer.writerows(csv_rows)

print(f"\nSaved: {txt_path}")
print(f"Saved: {csv_path}")
