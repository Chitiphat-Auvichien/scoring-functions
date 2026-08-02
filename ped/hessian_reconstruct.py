"""Framework 1: reconstruct the full Cartesian Hessian directly from
Gaussian's own printed normal-mode data -- no VEDA, no .chk/formchk, no
external QM package. A Gaussian freq log doesn't need to print a "Force
constants in Cartesian coordinates" block for this to work: given ALL
3N-6 vibrational frequencies and their full Cartesian displacement vectors,
the Hessian is exactly recoverable.

Theory (unchanged from the archived
ped/archive_python_ped/02_reconstruct_hessian.py):
  Gaussian's printed displacement vector l_k (per mode) is the correct
  EIGENVECTOR DIRECTION of the mass-weighted Hessian, but its overall
  normalization is Gaussian's own convention (not assumed here). What we
  actually need is L (Cartesian, per mode) satisfying L^T M L = I -- i.e.
  mass-weighted-orthonormal. Rather than trust a memorized Gaussian
  normalization formula, the correct per-mode rescaling is derived directly:
      L_k = l_k / sqrt(sum_i m_i |l_k,i|^2)
  which *by construction* satisfies L^T M L = I regardless of how Gaussian
  originally normalized l_k. (Mass units are irrelevant to the final result
  as long as used consistently -- plain amu is used throughout here.)

  Then: H_cart = M L Lambda L^T M,  Lambda = diag(omega_k^2),
  omega_k = 2 pi c_cm/s * nu_k(cm^-1)  (exact, standard relation).
"""
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                        # sibling ped/ modules (blocks.py)
sys.path.insert(0, os.path.join(_HERE, '..'))     # repo root (src/)
from src.parser import GaussianParser  # noqa: E402
from src.utils import atomicMass  # noqa: E402

import blocks  # noqa: E402

C_CM_PER_S = 2.99792458e10  # speed of light, cm/s (exact)


def load_geometry_and_modes(log_path):
    """Wrap GaussianParser(log_path).parse(parse_modes=True) and validate it
    the same way the archived 01_load_gaussian.py did: exactly 3N-6
    vibrational modes (nonlinear molecules only -- this is a known
    limitation carried over from the archived script, not a new one)."""
    parser = GaussianParser(log_path)
    data = parser.parse(parse_modes=True)

    symbols = data["atoms"]
    coords_ang = data["coords"]
    modes = data["modes"]
    bonds = data["bonds"]
    natm = len(symbols)

    n_expected = 3 * natm - 6
    if len(modes) != n_expected:
        raise ValueError(f"Expected {n_expected} vibrational modes (3N-6, "
                          f"N={natm}), parser returned {len(modes)} for "
                          f"{log_path}. (Framework 1/hessian_reconstruct.py "
                          "only supports nonlinear molecules.)")

    if any(m["reduced_mass"] is None for m in modes):
        raise ValueError("Some modes are missing a parsed reduced_mass -- "
                          "needed to identify Gaussian's displacement "
                          "normalization convention.")

    return {
        "atoms": symbols,
        "coords": coords_ang,
        "modes": modes,
        "bonds": bonds,
    }


def reconstruct_hessian_amu(symbols, l_raw, vibfreq, reduced_mass=None):
    """H_cart = M L Lambda L^T M, in plain amu*s^-2 units. `reduced_mass`
    is optional and only used for an informational (non-fatal) cross-check
    print against Gaussian's own printed reduced masses."""
    natm = len(symbols)
    Nvib = l_raw.shape[0]
    masses = np.array([atomicMass[s.lower()] for s in symbols])  # (N,), amu
    M = np.repeat(masses, 3)                                      # (3N,)

    L_flat = l_raw.reshape(Nvib, 3 * natm).T                      # (3N, Nvib)

    # --- Rescale each mode to satisfy L^T M L = I ---
    normsq = np.einsum('i,ik->k', M, L_flat ** 2)
    L = L_flat / np.sqrt(normsq)

    # --- Sanity check: mass-weighted orthonormality ---
    gram = (L * M[:, None]).T @ L
    off_diag_max = np.max(np.abs(gram - np.eye(Nvib)))
    print(f"Mass-weighted orthonormality check (L^T M L vs I): "
          f"max |error| = {off_diag_max:.3e}  "
          f"({'OK' if off_diag_max < 1e-3 else 'WARNING: modes not orthonormal'} -- "
          f"residual mainly from finite print precision on degenerate "
          f"E-symmetry mode pairs; the frequency round-trip is the "
          f"authoritative correctness check)")

    # --- Informational cross-check against Gaussian's own printed reduced
    #     mass (not required for correctness) ---
    if reduced_mass is not None:
        mu_candidate = np.einsum('i,ik->k', masses.repeat(3), L_flat ** 2)
        rel_err = np.abs(mu_candidate - reduced_mass) / reduced_mass
        print(f"\nReduced-mass cross-check (informational only, not required):")
        print(f"  max relative error vs Gaussian's printed reduced_mass: "
              f"{rel_err.max():.2e}  (mean {rel_err.mean():.2e})")

    # --- Build Lambda from frequencies ---
    omega = 2 * np.pi * C_CM_PER_S * vibfreq   # rad/s
    lam = omega ** 2

    # --- Reconstruct Cartesian Hessian: H_cart = M L Lambda L^T M ---
    H_cart = (M[:, None] * L) @ (lam[:, None] * L.T)
    H_cart = H_cart * M[None, :]
    H_cart = 0.5 * (H_cart + H_cart.T)

    return H_cart, L


def hessian_amu_to_hartree_bohr2(H_cart_amu):
    """H_cart_amu carries no length unit (see module docstring), so only a
    mass/time conversion is needed: 1 amu/s^2 = u*a0^2/Eh Hartree/Bohr^2."""
    import scipy.constants as sc

    u = sc.u
    a0 = sc.value('Bohr radius')
    Eh = sc.value('Hartree energy')
    factor = u * a0 ** 2 / Eh

    H_AU = 0.5 * (H_cart_amu + H_cart_amu.T) * factor
    diag = np.diag(H_AU)
    print(f"H_AU diagonal range: [{diag.min():.4f}, {diag.max():.4f}] Hartree/Bohr^2 "
          f"(A1-91.FMT reference range: ~0.03-0.7)")
    return H_AU


def build_hessian(mol, symbols, coords_ang, vibfreq, masses, reduced_mass, l_raw):
    """Orchestrate Framework 1: reconstruct -> convert units -> validate.
    Returns (H_AU, max_diff_cm1). Raises ValueError (via
    blocks.round_trip_validate_au) if the round-trip discrepancy is too
    large -- refusing to hand back an unverified Hessian."""
    print(f"Reconstructing Hessian for {mol.name} from {mol.log_path} "
          "(Framework 1: reconstruct from printed normal modes)")
    H_cart_amu, L = reconstruct_hessian_amu(symbols, l_raw, vibfreq,
                                             reduced_mass=reduced_mass)
    H_AU = hessian_amu_to_hartree_bohr2(H_cart_amu)
    max_diff = blocks.round_trip_validate_au(H_AU, masses, vibfreq)
    return H_AU, max_diff
