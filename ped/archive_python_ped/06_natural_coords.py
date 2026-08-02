"""
Step 6 (VEDA-style cross-check, additive -- does not touch 01-05):
Build a non-redundant (exactly 3N-6 = 30) "natural internal coordinate" set
for benzene's ring from the already-validated redundant 42x36 B-matrix
(ped/03_build_internal_coords.py's B.npy), following the symmetrized-
linear-combination principle of natural/local-symmetry coordinates
(Pulay, Fogarasi, Pang & Boggs, JACS 1979, 101, 2550 -- 'RN12' in
MolecularMotion.bib; ring-redundancy treatment per Cremer & Pople, JACS
1975, 97, 1354) rather than VEDA's own unpublished exact procedure.

VEDA's documented automatic rule ((N-1) stretches + (N-2) bends +
(N-3) torsions = 3N-6) does not apply directly to a ring: rings introduce
"ring-closure" redundancy that isn't resolved by simply omitting one bond
(VEDA's own documentation flags rings as needing special "ring coordinate"
handling). The approach here:

  1. Symmetrize each of the 7 raw coordinate families (CC-stretch,
     CH-stretch, CCC-bend, CCH-bend split into per-carbon symmetric/
     antisymmetric pairs first, ring-torsion, CH-wag) via the real C6
     cyclic-group (Fourier) basis -- a pure orthonormal basis rotation,
     same row space as the raw 42 coordinates, not a reduction yet.
  2. Rank-revealing greedy selection down to exactly 30 independent rows,
     processing candidates in a priority order (never-redundant
     substituent-type combinations first, totally-symmetric/ring-closure-
     suspect combinations last) but letting numerical rank (SVD-based) be
     the final, authoritative arbiter of what's actually redundant -- not
     an asserted-by-hand count (a by-hand identity count for the bend
     family came up short by one; the numerics settle it, not prose).

Requires: B.npy, coord_labels.txt (step 03)
Output: B_nat.npy (30x36), nat_coord_labels.txt (label, family, species,
        kept/dropped-from index)
"""
import os
import numpy as np

HERE = os.path.dirname(__file__)
TOL = 1e-8
TARGET_RANK = 30

B = np.load(os.path.join(HERE, 'B.npy'))            # (42, 36)
raw_labels = [l.split('\t')[0] for l in open(os.path.join(HERE, 'coord_labels.txt'))]

# --- Index ranges, matching 03_build_internal_coords.py's construction order ---
IDX = {
    'CC_stretch': list(range(0, 6)),
    'CH_stretch': list(range(6, 12)),
    'CCC_bend': list(range(12, 18)),
    'CCH_bend_a': list(range(18, 30, 2)),   # 18,20,22,24,26,28
    'CCH_bend_b': list(range(19, 30, 2)),   # 19,21,23,25,27,29
    'ring_torsion': list(range(30, 36)),
    'CH_wag_oop': list(range(36, 42)),
}
for k, v in IDX.items():
    assert len(v) == 6, f"{k}: expected 6 indices, got {len(v)}"
assert all(raw_labels[i] == 'CC_stretch' for i in IDX['CC_stretch'])
assert all(raw_labels[i] == 'CH_stretch' for i in IDX['CH_stretch'])
assert all(raw_labels[i] == 'CCC_bend' for i in IDX['CCC_bend'])
assert all(raw_labels[i] == 'CCH_bend' for i in IDX['CCH_bend_a'] + IDX['CCH_bend_b'])
assert all(raw_labels[i] == 'ring_torsion' for i in IDX['ring_torsion'])
assert all(raw_labels[i] == 'CH_wag_oop' for i in IDX['CH_wag_oop'])

# --- Local (per-carbon) symmetric/antisymmetric CCH combination, before
#     the ring-level Fourier step -- the standard two-stage natural-
#     coordinate construction (local combination, then ring symmetrization) ---
B_ccha, B_cchb = B[IDX['CCH_bend_a']], B[IDX['CCH_bend_b']]
B_cch_sym = (B_ccha + B_cchb) / np.sqrt(2)
B_cch_anti = (B_ccha - B_cchb) / np.sqrt(2)

families = {
    'CC_stretch': B[IDX['CC_stretch']],
    'CH_stretch': B[IDX['CH_stretch']],
    'CCC_bend': B[IDX['CCC_bend']],
    'CCH_bend_sym': B_cch_sym,
    'CCH_bend_anti': B_cch_anti,
    'ring_torsion': B[IDX['ring_torsion']],
    'CH_wag_oop': B[IDX['CH_wag_oop']],
}


def fourier_basis(n=6):
    """Real, orthonormal C6 cyclic-group basis. Returns (species_labels,
    coeff matrix) with coeff @ coeff.T == I (verified below, not assumed)."""
    j = np.arange(n)
    species, rows = [], []
    # k=0: totally symmetric 'A'
    species.append('A'); rows.append(np.full(n, 1.0 / np.sqrt(n)))
    # 0<k<n/2: doubly-degenerate 'E_k' cos/sin pair
    for k in range(1, n // 2):
        species.append(f'E{k}_cos'); rows.append(np.sqrt(2.0 / n) * np.cos(2 * np.pi * k * j / n))
        species.append(f'E{k}_sin'); rows.append(np.sqrt(2.0 / n) * np.sin(2 * np.pi * k * j / n))
    # k=n/2: alternating 'B' (n even)
    species.append('B'); rows.append(np.array([(-1.0) ** jj for jj in j]) / np.sqrt(n))
    coeff = np.array(rows)
    resid = np.max(np.abs(coeff @ coeff.T - np.eye(n)))
    assert resid < 1e-12, f"Fourier basis not orthonormal, residual={resid:.2e}"
    return species, coeff


species_labels, fourier_coeff = fourier_basis(6)   # 6 species: A, E1_cos, E1_sin, E2_cos, E2_sin, B

sym_rows, sym_labels, sym_family, sym_species = [], [], [], []
for fam_name, fam_B in families.items():
    fam_resid = np.max(np.abs(fourier_coeff @ fourier_coeff.T - np.eye(6)))
    B_sym = fourier_coeff @ fam_B          # (6, 36), pure basis rotation
    # verify same row space (rank-preserving rotation, not a reduction)
    rank_before = np.linalg.matrix_rank(fam_B, tol=TOL)
    rank_after = np.linalg.matrix_rank(B_sym, tol=TOL)
    assert rank_before == rank_after, (
        f"{fam_name}: symmetrization changed rank ({rank_before} -> {rank_after})")
    for sp, row in zip(species_labels, B_sym):
        sym_rows.append(row)
        sym_labels.append(f"{fam_name}_{sp}")
        sym_family.append(fam_name)
        sym_species.append(sp)

B_sym_all = np.array(sym_rows)   # (42, 36), symmetrized, same span as raw B
assert B_sym_all.shape == (42, 36)

# --- Priority order for rank-revealing greedy selection ---
# Tier 1 (never part of ring closure -- solid claim): substituent-type
# combinations -- CH_stretch (all 6), CH_wag_oop (all 6), and the non-
# totally-symmetric components of the local CCH sym/anti combinations.
# Tier 2: non-totally-symmetric CC_stretch / CCC_bend combinations.
# Tier 3 (most likely ring-closure-redundant): totally-symmetric CC_stretch
# and CCC_bend, totally-symmetric CCH sym/anti, and the full ring_torsion
# family. This ordering only determines WHICH 12 of 42 end up dropped when
# redundancy is hit -- the numerical rank check is the actual arbiter.


def indices_of(label_pred):
    return [i for i, lab in enumerate(sym_labels) if label_pred(lab)]


tier1 = (indices_of(lambda l: l.startswith('CH_stretch_'))
         + indices_of(lambda l: l.startswith('CH_wag_oop_'))
         + indices_of(lambda l: l.startswith('CCH_bend_sym_') and not l.endswith('_A'))
         + indices_of(lambda l: l.startswith('CCH_bend_anti_') and not l.endswith('_A')))
tier2 = (indices_of(lambda l: l.startswith('CC_stretch_') and not l.endswith('_A'))
         + indices_of(lambda l: l.startswith('CCC_bend_') and not l.endswith('_A')))
tier3 = (indices_of(lambda l: l == 'CC_stretch_A')
         + indices_of(lambda l: l == 'CCC_bend_A')
         + indices_of(lambda l: l == 'CCH_bend_sym_A')
         + indices_of(lambda l: l == 'CCH_bend_anti_A')
         + indices_of(lambda l: l.startswith('ring_torsion_')))
priority_order = tier1 + tier2 + tier3
assert sorted(priority_order) == list(range(42)), "priority_order must cover all 42 candidates exactly once"
assert len(tier1) == 22 and len(tier2) == 10 and len(tier3) == 10

# --- Rank-revealing greedy selection: max-residual-norm (pivoted-QR-style)
#     pick WITHIN each priority tier, not a naive fixed-order pass/fail
#     against a drifting per-stack tolerance. A first attempt using a fixed
#     priority order with an incremental per-stack SVD check reached 30 rows
#     but left the result ill-conditioned (cond ~1e8, several near-1e-8
#     singular values) -- a real sign that some accepted rows were only
#     marginally independent. Max-residual selection is the standard,
#     numerically robust rank-revealing procedure (equivalent to modified
#     Gram-Schmidt / pivoted QR); restricting it to run tier-by-tier keeps
#     the chemically-motivated preference order (substituent coordinates
#     before ring-closure-suspect ones) while letting the numerics pick the
#     best-conditioned representative within each tier.
REF_SCALE = np.linalg.svd(B_sym_all, compute_uv=False)[0]   # fixed reference, not per-stack
ABS_TOL = TOL * REF_SCALE

kept, dropped = [], []
Q = np.zeros((0, 36))   # orthonormal basis of the accepted-row subspace so far

for tier in (tier1, tier2, tier3):
    remaining = list(tier)
    while remaining and len(kept) < TARGET_RANK:
        best_idx, best_row, best_resid_norm = None, None, -1.0
        for idx in remaining:
            row = B_sym_all[idx]
            resid = row - Q.T @ (Q @ row) if Q.shape[0] else row.copy()
            n = np.linalg.norm(resid)
            if n > best_resid_norm:
                best_idx, best_row, best_resid_norm = idx, resid, n
        if best_resid_norm > ABS_TOL:
            Q = np.vstack([Q, best_row / best_resid_norm])
            kept.append(best_idx)
        else:
            dropped.append(best_idx)
        remaining.remove(best_idx)
    dropped.extend(remaining)   # tier exhausted without reaching target rank, or target already hit

remaining = [i for i in priority_order if i not in kept and i not in dropped]
dropped += remaining
accepted_rows = [B_sym_all[i] for i in kept]

if len(kept) != TARGET_RANK:
    raise SystemExit(f"FATAL: selected {len(kept)} independent coordinates, expected {TARGET_RANK}")

B_nat = np.array(accepted_rows)   # (30, 36)
cond = np.linalg.cond(B_nat)
svals_nat = np.linalg.svd(B_nat, compute_uv=False)

print(f"Kept {len(kept)} / 42 symmetrized coordinates as the natural non-redundant set:")
for i in kept:
    print(f"  KEEP    {sym_labels[i]}")
print(f"\nDropped {len(dropped)} / 42 (ring-closure / redundant directions):")
for i in dropped:
    print(f"  DROP    {sym_labels[i]}")
print(f"\nB_nat condition number: {cond:.6e}")
print(f"B_nat singular values: {np.array2string(svals_nat, precision=4)}")

np.save(os.path.join(HERE, 'B_nat.npy'), B_nat)
with open(os.path.join(HERE, 'nat_coord_labels.txt'), 'w') as f:
    for i in kept:
        f.write(f"{sym_labels[i]}\t{sym_family[i]}\t{sym_species[i]}\n")

print("\nSaved: B_nat.npy, nat_coord_labels.txt")
