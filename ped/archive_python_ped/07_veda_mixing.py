"""
Step 7 (VEDA-style cross-check): greedy EPm-maximizing coordinate mixing,
implementing VEDA's documented "mix coordinates, accept if EPm improves"
optimization paradigm (Jamroz, Spectrochim. Acta A 2013 -- objective and
accept-if-improves rule are published; the exact search strategy below is
NOT VEDA's own, since that isn't public -- see VEDA_STYLE_METHODOLOGY.md).

Design choices made here, documented explicitly (not VEDA's own published
algorithm):
  - Moves are pairwise Givens rotations between two rows of the current
    coordinate matrix B: guarantees the working set stays exactly rank-30
    by construction (a rotation cannot change the row space), so there is
    no risk of accidentally reintroducing redundancy during optimization.
  - Candidate pairs are restricted to the SAME family (same original
    coordinate type tag from nat_coord_labels.txt) -- e.g. only mix two
    CH-stretch-derived rows together, never a stretch with a bend. VEDA's
    own description is "superposition of local modes" of compatible type;
    an unrestricted global search risks producing EPm-maximizing but
    chemically uninterpretable composite coordinates that can't be cleanly
    assigned to a PED category afterward.
  - Deterministic sweep order (family by family, i<j within each family),
    full-sweep convergence (a sweep with zero accepted moves stops the
    optimization) or a 50-sweep cap.

Requires: B_nat.npy, nat_coord_labels.txt (step 06); H_cart.npy, L.npy,
          vibfreq.npy (step 02)
Output: B_mixed.npy (30x36), M_mixed.npy (30x30, coefficients vs. the
        natural basis), epm_history.npy/.txt
"""
import os
import sys
import numpy as np
from scipy.optimize import minimize_scalar

sys.path.insert(0, os.path.dirname(__file__))
from _ped_core import compute_PED, freq_to_lambda, EPm

HERE = os.path.dirname(__file__)
MAX_SWEEPS = 50
TOL_IMPROVE = 1e-9

B_nat = np.load(os.path.join(HERE, 'B_nat.npy'))
H = np.load(os.path.join(HERE, 'H_cart.npy'))
L = np.load(os.path.join(HERE, 'L.npy'))
vib = np.load(os.path.join(HERE, 'vibfreq.npy'))
lam = freq_to_lambda(vib)

nat_labels, nat_families = [], []
for line in open(os.path.join(HERE, 'nat_coord_labels.txt')):
    label, family, species = line.rstrip('\n').split('\t')
    nat_labels.append(label)
    nat_families.append(family)

Nnat = B_nat.shape[0]
assert Nnat == 30

# --- same-family candidate pairs only ---
family_groups = {}
for i, fam in enumerate(nat_families):
    family_groups.setdefault(fam, []).append(i)
pairs = []
for fam, idxs in family_groups.items():
    for a in range(len(idxs)):
        for b in range(a + 1, len(idxs)):
            pairs.append((idxs[a], idxs[b]))
print(f"Families and sizes: { {k: len(v) for k, v in family_groups.items()} }")
print(f"Same-family candidate pairs: {len(pairs)}")

B = B_nat.copy()
M = np.eye(Nnat)   # tracks each current row's composition vs. the natural basis

epm0 = EPm(compute_PED(B, H, L, lam))
epm_history = [epm0]
print(f"Initial EPm (natural, unmixed coordinates): {epm0:.6f}")

for sweep in range(1, MAX_SWEEPS + 1):
    improved_this_sweep = False
    for (i, j) in pairs:
        def neg_epm(theta):
            c, s = np.cos(theta), np.sin(theta)
            Bt = B.copy()
            Bt[i], Bt[j] = c * B[i] + s * B[j], -s * B[i] + c * B[j]
            return -EPm(compute_PED(Bt, H, L, lam))

        current_epm = epm_history[-1]
        res = minimize_scalar(neg_epm, bounds=(-np.pi / 2, np.pi / 2), method='bounded')
        if -res.fun > current_epm + TOL_IMPROVE:
            c, s = np.cos(res.x), np.sin(res.x)
            Bi, Bj = B[i].copy(), B[j].copy()
            Mi, Mj = M[i].copy(), M[j].copy()
            B[i], B[j] = c * Bi + s * Bj, -s * Bi + c * Bj
            M[i], M[j] = c * Mi + s * Mj, -s * Mi + c * Mj
            epm_history.append(-res.fun)
            improved_this_sweep = True
    print(f"Sweep {sweep}: EPm = {epm_history[-1]:.6f}  "
          f"({'accepted moves' if improved_this_sweep else 'no improving moves -- converged'})")
    if not improved_this_sweep:
        break
else:
    print(f"WARNING: reached the {MAX_SWEEPS}-sweep cap without full convergence.")

PED_final = compute_PED(B, H, L, lam)
colsums = PED_final.sum(axis=0)
print(f"\nFinal EPm: {epm_history[-1]:.6f}  (initial was {epm0:.6f}, "
      f"delta={epm_history[-1]-epm0:+.6f})")
print(f"PED column sums after mixing: min={colsums.min():.6f}, max={colsums.max():.6f}")

np.save(os.path.join(HERE, 'B_mixed.npy'), B)
np.save(os.path.join(HERE, 'M_mixed.npy'), M)
np.save(os.path.join(HERE, 'epm_history.npy'), np.array(epm_history))
with open(os.path.join(HERE, 'epm_history.txt'), 'w') as f:
    for k, v in enumerate(epm_history):
        f.write(f"{k}\t{v:.8f}\n")

print("\nSaved: B_mixed.npy, M_mixed.npy, epm_history.npy/.txt")
