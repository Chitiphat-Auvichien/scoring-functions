"""
Shared PED core, factored out of 04_compute_ped.py's formula so the
VEDA-style scripts (06-10) can reuse it without duplicating or touching
04_compute_ped.py itself.

Theory (same as 04_compute_ped.py -- see PED_METHODOLOGY.md):
  F_q = (B+)^T H (B+)           (Pulay generalized-inverse; B+ becomes an
                                  exact inverse-like operator when B has
                                  full row rank, so this formula is valid
                                  for both redundant and non-redundant B)
  D_mu = B @ L_mu
  PED[n, mu] = D[n,mu] * (F_q @ D[:,mu])[n] / lambda_mu
"""
import numpy as np

C_CM_PER_S = 2.99792458e10  # speed of light, cm/s (exact)


def freq_to_lambda(vibfreq_cm1):
    """Convert Gaussian's cm^-1 frequencies to lambda_k = omega_k^2 (rad/s)^2,
    the same convention used throughout ped/02_reconstruct_hessian.py and
    ped/04_compute_ped.py."""
    omega = 2 * np.pi * C_CM_PER_S * vibfreq_cm1
    return omega ** 2


def compute_PED(B, H, L, lam, rcond=1e-8):
    """B: (Nint, 3N) internal-coordinate B-matrix (redundant or not).
    H: (3N, 3N) Cartesian Hessian. L: (3N, Nvib) mass-weighted-orthonormal
    Cartesian mode matrix (L^T M L = I). lam: (Nvib,) eigenvalues omega_k^2.
    Returns PED_raw: (Nint, Nvib), raw (unnormalized-to-100%) PED matrix;
    column sums should equal 1.0 for a correct B/H/L/lam combination."""
    Bpinv = np.linalg.pinv(B, rcond=rcond)
    Fq = Bpinv.T @ H @ Bpinv
    D = B @ L
    Nint, Nvib = B.shape[0], L.shape[1]
    PED_raw = np.zeros((Nint, Nvib))
    for mu in range(Nvib):
        FqD = Fq @ D[:, mu]
        PED_raw[:, mu] = D[:, mu] * FqD / lam[mu]
    return PED_raw


def EPm(PED_raw):
    """VEDA's documented optimization objective: the sum, over all modes,
    of the largest single-coordinate |PED| contribution in that mode.
    abs() is our own explicit convention (PED components can be negative
    for non-orthogonal/redundant coordinate choices -- see
    PED_METHODOLOGY.md Sec 4) -- VEDA's own sign convention inside EPm is
    not published, so this is documented here rather than silently assumed."""
    return float(np.sum(np.max(np.abs(PED_raw), axis=0)))
