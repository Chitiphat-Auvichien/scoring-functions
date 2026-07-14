"""
Step 4: Transform the Cartesian Hessian into the redundant internal-
coordinate force constant matrix F_q (Pulay's generalized-inverse
method), project each normal mode onto internal coordinates, and
compute the Potential Energy Distribution (PED).

Theory:
  At a stationary point, H = B^T F_q B exactly.
  For a redundant (non-square) B, use the Moore-Penrose pseudoinverse
  B+ instead of a true inverse:
      F_q = (B+)^T H (B+)
  For each normal mode mu (Cartesian displacement L_mu, eigenvalue
  lambda_mu = omega_mu^2):
      D_mu = B @ L_mu                     (internal-coordinate projection)
      PED[n, mu] = D[n,mu] * (F_q @ D[:,mu])[n] / lambda_mu
  Sum_n PED[n,mu] should equal 1.0 for every mode -- this is the
  correctness check performed below. (lambda's absolute scale/units
  are irrelevant to this ratio -- only its per-mode value matters --
  so amu/Angstrom-based units from steps 1-2 are fine as-is.)

Requires: H_cart.npy, L.npy, vibfreq.npy (step 2); B.npy, coord_labels.txt (step 3)
Output: PED_group_pct.npy (n_categories x 30, in %), cats.txt
"""
import os
import numpy as np

HERE = os.path.dirname(__file__)
C_CM_PER_S = 2.99792458e10  # speed of light, cm/s (exact)

H = np.load(os.path.join(HERE, 'H_cart.npy'))
L = np.load(os.path.join(HERE, 'L.npy'))
vib = np.load(os.path.join(HERE, 'vibfreq.npy'))
B = np.load(os.path.join(HERE, 'B.npy'))
labels = [l.split('\t')[0] for l in open(os.path.join(HERE, 'coord_labels.txt'))]

Nint, Nvib = B.shape[0], L.shape[1]

Bpinv = np.linalg.pinv(B, rcond=1e-8)      # (36, 42)
Fq = Bpinv.T @ H @ Bpinv                   # (42, 42)
D = B @ L                                   # (42, 30)

omega = 2 * np.pi * C_CM_PER_S * vib        # rad/s, same convention as step 2
lam = omega ** 2

PED_raw = np.zeros((Nint, Nvib))
for mu in range(Nvib):
    FqD = Fq @ D[:, mu]
    PED_raw[:, mu] = D[:, mu] * FqD / lam[mu]

colsums = PED_raw.sum(axis=0)
print("PED column sums (correctness check -- should all be ~1.0):")
print(np.round(colsums, 3))
if np.any(np.abs(colsums - 1) > 0.05):
    print("WARNING: some modes deviate from 1.0 -- check coordinate "
          "completeness / redundancy handling.")

cats = sorted(set(labels), key=labels.index)
cat_idx = {c: [i for i, l in enumerate(labels) if l == c] for c in cats}
PED_group = np.array([PED_raw[cat_idx[c], :].sum(axis=0) for c in cats])
PED_group_pct = 100 * PED_group / PED_group.sum(axis=0, keepdims=True)

np.save(os.path.join(HERE, 'PED_group_pct.npy'), PED_group_pct)
with open(os.path.join(HERE, 'cats.txt'), 'w') as f:
    for c in cats:
        f.write(c + '\n')

print("\nCategories:", cats)
print("Saved: PED_group_pct.npy, cats.txt")
