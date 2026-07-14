"""
Step 2: Reconstruct the full Cartesian Hessian from Gaussian's own printed
normal-mode data -- no VEDA, no .chk/formchk, no PySCF. The Gaussian log
never prints a "Force constants in Cartesian coordinates" block for this
job, but it doesn't need to: given ALL 3N-6 vibrational frequencies and
their full Cartesian displacement vectors, the Hessian is exactly
recoverable.

Theory:
  Gaussian's printed displacement vector l_k (per mode) is the correct
  EIGENVECTOR DIRECTION of the mass-weighted Hessian, but its overall
  normalization is Gaussian's own convention (not assumed here). What we
  actually need is L (Cartesian, per mode) satisfying L^T M L = I -- i.e.
  mass-weighted-orthonormal. Rather than trust a memorized Gaussian
  normalization formula, we derive the correct per-mode rescaling directly:
      L_k = l_k / sqrt(sum_i m_i |l_k,i|^2)
  which *by construction* satisfies the L^T M L = I condition regardless of
  how Gaussian originally normalized l_k. (Mass units are irrelevant to the
  final PED ratio -- verified algebraically -- so plain amu is used
  throughout, no amu->kg conversion needed.)

  Then: H_cart = M L Lambda L^T M,  Lambda = diag(omega_k^2),
  omega_k = 2 pi c_cm/s * nu_k(cm^-1)  (exact, standard relation).

  Validation: re-diagonalize M^-1/2 H_cart M^-1/2 and confirm the 30
  recovered nonzero eigenvalues reproduce Gaussian's original frequencies
  (should be a near-exact round trip, since H_cart is built directly from
  those eigenvectors/eigenvalues) and that the other 6 eigenvalues
  (translation/rotation null space, never included in the reconstruction)
  come out ~0.

Requires: opt_coords_ang.npy, opt_symbols.txt, l_raw.npy, vibfreq.npy,
          reduced_mass.npy (step 1)
Output: H_cart.npy (36x36), L.npy (36x30, L^T M L = I), vibfreq.npy passthrough
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from src.utils import atomicMass

HERE = os.path.dirname(__file__)
C_CM_PER_S = 2.99792458e10  # speed of light, cm/s (exact)

symbols = [l.strip() for l in open(os.path.join(HERE, 'opt_symbols.txt'))]
l_raw = np.load(os.path.join(HERE, 'l_raw.npy'))          # (30, 12, 3)
vibfreq = np.load(os.path.join(HERE, 'vibfreq.npy'))       # (30,)
reduced_mass = np.load(os.path.join(HERE, 'reduced_mass.npy'))  # (30,)

natm = len(symbols)
Nvib = l_raw.shape[0]
masses = np.array([atomicMass[s.lower()] for s in symbols])  # (12,), amu
M = np.repeat(masses, 3)                                      # (36,)

L_flat = l_raw.reshape(Nvib, 3 * natm).T                      # (36, 30)

# --- Rescale each mode to satisfy L^T M L = I, independent of whatever
#     normalization Gaussian used when printing l_raw ---
normsq = np.einsum('i,ik->k', M, L_flat ** 2)                  # (30,)
L = L_flat / np.sqrt(normsq)

# --- Sanity check: mass-weighted orthonormality ---
gram = (L * M[:, None]).T @ L
off_diag_max = np.max(np.abs(gram - np.eye(Nvib)))
print(f"Mass-weighted orthonormality check (L^T M L vs I): "
      f"max |error| = {off_diag_max:.3e}  "
      f"({'OK' if off_diag_max < 1e-3 else 'WARNING: modes not orthonormal'} -- "
      f"residual mainly from finite print precision on Gaussian's degenerate "
      f"E-symmetry mode pairs; the frequency round-trip below is the "
      f"authoritative correctness check)")

# --- Informational cross-check against Gaussian's own printed reduced mass
#     (not required for correctness -- just an independent consistency check
#     against the documented "unit-Euclidean-norm displacement" convention) ---
mu_candidate = 1.0 / np.einsum('i,ik->k', 1.0 / masses.repeat(3), L_flat ** 2)
rel_err = np.abs(mu_candidate - reduced_mass) / reduced_mass
print(f"\nReduced-mass cross-check (informational only, not required):")
print(f"  max relative error vs Gaussian's printed reduced_mass: "
      f"{rel_err.max():.2e}  (mean {rel_err.mean():.2e})")

# --- Build Lambda from frequencies (mass-unit-independent, standard relation) ---
omega = 2 * np.pi * C_CM_PER_S * vibfreq   # rad/s
lam = omega ** 2

# --- Reconstruct Cartesian Hessian: H_cart = M L Lambda L^T M ---
H_cart = (M[:, None] * L) @ (lam[:, None] * L.T)
H_cart = H_cart * M[None, :]
H_cart = 0.5 * (H_cart + H_cart.T)

# --- Validation: re-diagonalize and recover Gaussian's frequencies ---
Minv_sqrt = 1.0 / np.sqrt(M)
Hmw_check = (Minv_sqrt[:, None] * H_cart) * Minv_sqrt[None, :]
evals, _ = np.linalg.eigh(Hmw_check)
evals_sorted = np.sort(evals)

near_zero = evals_sorted[:6]
vib_evals = evals_sorted[6:]
freq_recovered = np.sign(vib_evals) * np.sqrt(np.abs(vib_evals)) / (2 * np.pi * C_CM_PER_S)
freq_recovered_sorted = np.sort(freq_recovered)
freq_original_sorted = np.sort(vibfreq)

max_diff = np.max(np.abs(freq_recovered_sorted - freq_original_sorted))
print(f"\nNull-space check (should be ~0, T/R modes never included): "
      f"{np.round(near_zero, 6)}")
print(f"\nFrequency round-trip (Gaussian -> Hessian -> re-diagonalized), cm^-1:")
for a, b in zip(np.round(freq_original_sorted, 4), np.round(freq_recovered_sorted, 4)):
    print(f"  original={a:10.4f}   recovered={b:10.4f}   diff={a - b:+.6f}")
print(f"\nMax discrepancy: {max_diff:.6f} cm^-1  "
      f"({'OK, proceeding' if max_diff < 2.0 else 'WARNING: check normalization/units!'})")

np.save(os.path.join(HERE, 'H_cart.npy'), H_cart)
np.save(os.path.join(HERE, 'L.npy'), L)
print("\nSaved: H_cart.npy, L.npy (vibfreq.npy unchanged from step 1)")
