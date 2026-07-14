"""
Step 3: Define a redundant set of internal coordinates for benzene and
build the Wilson B-matrix (B[n,i] = d(q_n)/d(x_i)) by numerical central
differences, in Angstrom (Gaussian's native geometry units).

Ring carbons and their attached hydrogens are detected generically from the
parsed connectivity graph (bonds.txt from step 1) rather than assumed from
a hardcoded atom order -- makes this robust to whatever atom ordering the
source .com/.log file happens to use.

Internal coordinate set (42 total, redundant for 30 vibrational DOF):
  6  C-C stretches
  6  C-H stretches
  6  C-C-C in-plane bond angles
  12 C-C-H in-plane bond angles (2 per ring carbon)
  6  C-C-C-C ring torsions
  6  H-C-C-C out-of-plane wags

Requires: opt_coords_ang.npy, opt_symbols.txt, bonds.txt (from step 1)
Output: B.npy (42x36), coord_labels.txt
"""
import os
import numpy as np

HERE = os.path.dirname(__file__)

coords = np.load(os.path.join(HERE, 'opt_coords_ang.npy'))
symbols = [l.strip() for l in open(os.path.join(HERE, 'opt_symbols.txt'))]
natm = len(symbols)

bonds = []
with open(os.path.join(HERE, 'bonds.txt')) as f:
    for line in f:
        a, b = line.split()
        bonds.append((int(a), int(b)))

# --- Generic ring + substituent detection from the bonds graph ---
carbons = [i for i, s in enumerate(symbols) if s == 'C']
adj = {i: set() for i in range(natm)}
for a, b in bonds:
    adj[a].add(b)
    adj[b].add(a)

carbon_adj = {c: sorted(n for n in adj[c] if n in carbons) for c in carbons}
if any(len(v) != 2 for v in carbon_adj.values()):
    raise ValueError("Expected each ring carbon to have exactly 2 carbon "
                      "neighbors (6-membered ring) -- got: " +
                      str({c: v for c, v in carbon_adj.items() if len(v) != 2}))

# Walk the ring starting from the first carbon
ring = [carbons[0]]
prev = None
cur = carbons[0]
while len(ring) < len(carbons):
    nxt = [n for n in carbon_adj[cur] if n != prev][0]
    ring.append(nxt)
    prev, cur = cur, nxt
if ring[0] not in carbon_adj[ring[-1]]:
    raise ValueError("Ring walk did not close -- carbons are not a single 6-ring.")
C = ring  # C[i] = ring carbon index in cyclic order

Hh = []
for c in C:
    h_neighbors = [n for n in adj[c] if symbols[n] != 'C']
    if len(h_neighbors) != 1:
        raise ValueError(f"Expected exactly 1 non-carbon substituent on ring "
                          f"carbon {c}, found {h_neighbors}")
    Hh.append(h_neighbors[0])

print(f"Detected ring (cyclic order): {C}")
print(f"Attached substituents:        {Hh}")


def bond(a, b, x):
    return np.linalg.norm(x[a] - x[b])


def angle(a, b, c, x):
    """Angle at atom b, between vectors b->a and b->c."""
    v1, v2 = x[a] - x[b], x[c] - x[b]
    v1n, v2n = v1 / np.linalg.norm(v1), v2 / np.linalg.norm(v2)
    return np.arccos(np.clip(np.dot(v1n, v2n), -1, 1))


def dihedral(a, b, c, d, x):
    b1, b2, b3 = x[b] - x[a], x[c] - x[b], x[d] - x[c]
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    return np.arctan2(np.dot(m1, n2), np.dot(n1, n2))


coord_defs = []
for i in range(6):
    coord_defs.append(('CC_stretch', bond, (C[i], C[(i + 1) % 6])))
for i in range(6):
    coord_defs.append(('CH_stretch', bond, (C[i], Hh[i])))
for i in range(6):
    prev_c, this_c, next_c = C[(i - 1) % 6], C[i], C[(i + 1) % 6]
    coord_defs.append(('CCC_bend', angle, (prev_c, this_c, next_c)))
for i in range(6):
    prev_c, this_c, next_c, h_ = C[(i - 1) % 6], C[i], C[(i + 1) % 6], Hh[i]
    coord_defs.append(('CCH_bend', angle, (next_c, this_c, h_)))
    coord_defs.append(('CCH_bend', angle, (h_, this_c, prev_c)))
for i in range(6):
    a, b, c, d = C[(i - 1) % 6], C[i], C[(i + 1) % 6], C[(i + 2) % 6]
    coord_defs.append(('ring_torsion', dihedral, (a, b, c, d)))
for i in range(6):
    prev_c, this_c, next_c, h_ = C[(i - 1) % 6], C[i], C[(i + 1) % 6], Hh[i]
    coord_defs.append(('CH_wag_oop', dihedral, (prev_c, this_c, next_c, h_)))

Nint = len(coord_defs)
print(f"Number of redundant internal coordinates: {Nint} "
      f"(3N-6 = {3*natm-6} vibrational DOF)")

x0 = coords.copy()
disp = 1e-4  # Angstrom
xflat = x0.flatten()
B = np.zeros((Nint, 3 * natm))

for i in range(3 * natm):
    xp, xm = xflat.copy(), xflat.copy()
    xp[i] += disp
    xm[i] -= disp
    Xp, Xm = xp.reshape(natm, 3), xm.reshape(natm, 3)
    qp = np.array([f(*idxs, Xp) for (_, f, idxs) in coord_defs])
    qm = np.array([f(*idxs, Xm) for (_, f, idxs) in coord_defs])
    dq = qp - qm
    dq = (dq + np.pi) % (2 * np.pi) - np.pi   # handle angle wraparound
    B[:, i] = dq / (2 * disp)

np.save(os.path.join(HERE, 'B.npy'), B)
with open(os.path.join(HERE, 'coord_labels.txt'), 'w') as f:
    for (lab, func, idxs) in coord_defs:
        f.write(f"{lab}\t{idxs}\n")

print("B-matrix shape:", B.shape)
print("Saved: B.npy, coord_labels.txt")
