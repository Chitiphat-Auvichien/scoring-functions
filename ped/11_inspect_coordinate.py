"""
Step 11 (VEDA-style cross-check, diagnostic tool -- reads only, writes only
to ped/inspect_*.png): show what a final optimized coordinate actually IS
-- its composition back through mixing -> natural symmetrization -> raw
local internal coordinates, and its Cartesian atomic-displacement pattern
(the same kind of per-atom arrow picture used for normal-mode figures
elsewhere in this project), for a given mode (by frequency) or a given
final-coordinate index directly.

Usage:
    python 11_inspect_coordinate.py --freq 1056.39
    python 11_inspect_coordinate.py --coord 7
"""
import argparse
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(__file__)
sys.path.insert(0, HERE)
from _ped_core import compute_PED, freq_to_lambda

parser = argparse.ArgumentParser()
parser.add_argument('--freq', type=float, default=None, help='mode frequency (cm-1); dominant coordinate for this mode is inspected')
parser.add_argument('--coord', type=int, default=None, help='final coordinate index (0-29) to inspect directly')
args = parser.parse_args()
if (args.freq is None) == (args.coord is None):
    raise SystemExit("Pass exactly one of --freq or --coord")

B_final = np.load(os.path.join(HERE, 'B_final.npy'))
M_final = np.load(os.path.join(HERE, 'M_final.npy'))
B_nat = np.load(os.path.join(HERE, 'B_nat.npy'))
H = np.load(os.path.join(HERE, 'H_cart.npy'))
L = np.load(os.path.join(HERE, 'L.npy'))
vib = np.load(os.path.join(HERE, 'vibfreq.npy'))
coords = np.load(os.path.join(HERE, 'opt_coords_ang.npy'))
symbols = [l.strip() for l in open(os.path.join(HERE, 'opt_symbols.txt'))]
lam = freq_to_lambda(vib)

nat_labels, nat_families, nat_species = [], [], []
for line in open(os.path.join(HERE, 'nat_coord_labels.txt')):
    label, family, species = line.rstrip('\n').split('\t')
    nat_labels.append(label); nat_families.append(family); nat_species.append(species)

# --- raw B.npy / coord_labels.txt, to trace natural coords back to atom-level bends ---
B_raw = np.load(os.path.join(HERE, 'B.npy'))
raw_labels, raw_defs = [], []
for line in open(os.path.join(HERE, 'coord_labels.txt')):
    lab, idxs = line.rstrip('\n').split('\t')
    raw_labels.append(lab); raw_defs.append(idxs)
IDX_CCH_A = list(range(18, 30, 2))
IDX_CCH_B = list(range(19, 30, 2))


def fourier_basis(n=6):
    j = np.arange(n)
    species, rows = [], []
    species.append('A'); rows.append(np.full(n, 1.0 / np.sqrt(n)))
    for k in range(1, n // 2):
        species.append(f'E{k}_cos'); rows.append(np.sqrt(2.0 / n) * np.cos(2 * np.pi * k * j / n))
        species.append(f'E{k}_sin'); rows.append(np.sqrt(2.0 / n) * np.sin(2 * np.pi * k * j / n))
    species.append('B'); rows.append(np.array([(-1.0) ** jj for jj in j]) / np.sqrt(n))
    return species, np.array(rows)


four_species, four_coeff = fourier_basis(6)

if args.coord is not None:
    k = args.coord
    mode_note = ""
else:
    PED = compute_PED(B_final, H, L, lam)
    mode_idx = int(np.argmin(np.abs(vib - args.freq)))
    k = int(np.argmax(np.abs(PED[:, mode_idx])))
    mode_note = f" (dominant coordinate for the mode nearest {args.freq} cm-1, actual freq {vib[mode_idx]:.2f} cm-1, PED={PED[k, mode_idx]*100:.1f}%)"

print(f"Inspecting final coordinate #{k}{mode_note}")
print("=" * 70)

print("\n1. Composition vs. the natural (30-coordinate) basis (M_final row):")
comps = [(i, c) for i, c in enumerate(M_final[k]) if abs(c) > 1e-6]
comps.sort(key=lambda x: -abs(x[1]))
for i, c in comps:
    print(f"   {c:+.4f} x {nat_labels[i]:28s} (family={nat_families[i]}, species={nat_species[i]})")

print("\n2. Each natural-basis component, traced back to the raw per-carbon")
print("   CCH bend angles (family=CCH_bend_sym/anti; other families omitted here")
print("   if not present above):")
for i, c in comps:
    fam, sp = nat_families[i], nat_species[i]
    if fam not in ('CCH_bend_sym', 'CCH_bend_anti'):
        print(f"   [{fam}] not a CCH-bend family -- skipping atom-level trace for this component")
        continue
    sp_idx = four_species.index(sp)
    four_row = four_coeff[sp_idx]   # (6,) weight per ring-carbon index
    sign = +1 if fam == 'CCH_bend_sym' else -1
    print(f"   {nat_labels[i]} (species {sp}) = per-ring-carbon weights (a=+b combo if sym, a-b if anti):")
    for ring_i in range(6):
        a_idx, b_idx = IDX_CCH_A[ring_i], IDX_CCH_B[ring_i]
        a_def, b_def = raw_defs[a_idx], raw_defs[b_idx]
        print(f"      ring position {ring_i}: weight {four_row[ring_i]:+.4f}  "
              f"-> local bend a={a_def} (C-C-H one side), b={b_def} (C-C-H other side), "
              f"combo = a {'+ ' if sign>0 else '- '}b")

print(f"\n3. Cartesian atomic-displacement pattern for this coordinate")
print("   (B_final row -- how much and which direction each atom moves per")
print("   unit change in this internal coordinate):")
disp = B_final[k].reshape(-1, 3)
mags = np.linalg.norm(disp, axis=1)
order = np.argsort(-mags)
for idx in order:
    if mags[idx] < 0.01 * mags.max():
        continue
    print(f"   atom {idx:2d} ({symbols[idx]}): |d|={mags[idx]:.4f}  "
          f"direction=({disp[idx,0]:+.3f}, {disp[idx,1]:+.3f}, {disp[idx,2]:+.3f})")

# --- visualization: 2D projection (ring plane) with in-plane arrows ---
fig, ax = plt.subplots(figsize=(6, 6))
carbons = [i for i, s in enumerate(symbols) if s == 'C']
hydrogens = [i for i, s in enumerate(symbols) if s == 'H']
xy = coords[:, :2]
ax.scatter(xy[carbons, 0], xy[carbons, 1], s=300, c='#444444', zorder=3)
ax.scatter(xy[hydrogens, 0], xy[hydrogens, 1], s=150, c='#dddddd', edgecolors='#444444', zorder=3)
for i in range(len(symbols)):
    ax.annotate(f"{symbols[i]}{i}", xy[i], ha='center', va='center', fontsize=7, zorder=4)
scale = 1.5
for i in range(len(symbols)):
    dx, dy = disp[i, 0] * scale, disp[i, 1] * scale
    if np.hypot(dx, dy) > 0.02:
        ax.annotate('', xy=(xy[i, 0] + dx, xy[i, 1] + dy), xytext=(xy[i, 0], xy[i, 1]),
                     arrowprops=dict(arrowstyle='->', color='#D55E00', lw=2), zorder=5)
ax.set_aspect('equal')
ax.set_title(f"Final coordinate #{k}{mode_note}\n(in-plane XY component only; see printed z-components above)")
ax.axis('off')
out_png = os.path.join(HERE, f'inspect_coord_{k}.png')
fig.savefig(out_png, dpi=150, bbox_inches='tight')
print(f"\nSaved visualization: {out_png}")
