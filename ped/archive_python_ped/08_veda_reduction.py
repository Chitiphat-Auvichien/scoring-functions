"""
Step 8 (VEDA-style cross-check): reduction pass, implementing VEDA's
documented final step of eliminating components of multicomponent (mixed)
coordinates that change EPm "in minimal degree" (Jamroz 2013 -- the
existence/purpose of this step is published; the exact numeric threshold
VEDA itself uses is not, so the tolerance below is our own explicit,
documented choice -- see VEDA_STYLE_METHODOLOGY.md and the epm_tol
sensitivity sweep in the validation report).

For each of the 30 mixed coordinates, try dropping its smallest-magnitude
components (relative to the natural basis, from M_mixed) one at a time,
smallest first; accept the drop only if (a) the coefficient matrix stays
full rank and reasonably conditioned (cond<1e8) after re-normalizing, and
(b) the resulting EPm changes by less than epm_tol. This simplifies each
final coordinate toward fewer, more chemically interpretable components
without materially harming PED purity.

Requires: B_nat.npy (step 06); B_mixed.npy, M_mixed.npy (step 07);
          H_cart.npy, L.npy, vibfreq.npy (step 02)
Output: B_final.npy (30x36), M_final.npy (30x30)
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from _ped_core import compute_PED, freq_to_lambda, EPm

HERE = os.path.dirname(__file__)
DEFAULT_EPM_TOL = 0.01
COND_LIMIT = 1e8
ZERO_COEFF = 1e-6


def run_reduction(B_nat, B_mixed, M_mixed, H, L, lam, epm_tol=DEFAULT_EPM_TOL, verbose=True):
    B, M = B_mixed.copy(), M_mixed.copy()
    current = EPm(compute_PED(B, H, L, lam))
    n_dropped = 0
    for k in range(M.shape[0]):
        order = np.argsort(np.abs(M[k]))
        for comp in order:
            if abs(M[k, comp]) < ZERO_COEFF:
                continue
            Mt = M.copy()
            Mt[k, comp] = 0.0
            nrm = np.linalg.norm(Mt[k])
            if nrm < ZERO_COEFF:
                continue
            Mt[k] /= nrm
            if np.linalg.matrix_rank(Mt, tol=1e-10) < M.shape[0]:
                continue
            cond = np.linalg.cond(Mt)
            if cond > COND_LIMIT:
                continue
            Bt = Mt @ B_nat
            trial = EPm(compute_PED(Bt, H, L, lam))
            if abs(current - trial) < epm_tol:
                M[k], B[k], current = Mt[k], Bt[k], trial
                n_dropped += 1
    if verbose:
        print(f"epm_tol={epm_tol}: dropped {n_dropped} components, final EPm={current:.6f}")
    return B, M, current, n_dropped


if __name__ == '__main__':
    epm_tol = float(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_EPM_TOL

    B_nat = np.load(os.path.join(HERE, 'B_nat.npy'))
    B_mixed = np.load(os.path.join(HERE, 'B_mixed.npy'))
    M_mixed = np.load(os.path.join(HERE, 'M_mixed.npy'))
    H = np.load(os.path.join(HERE, 'H_cart.npy'))
    L = np.load(os.path.join(HERE, 'L.npy'))
    vib = np.load(os.path.join(HERE, 'vibfreq.npy'))
    lam = freq_to_lambda(vib)

    epm_before = EPm(compute_PED(B_mixed, H, L, lam))
    print(f"EPm before reduction: {epm_before:.6f}")

    B_final, M_final, epm_after, n_dropped = run_reduction(
        B_nat, B_mixed, M_mixed, H, L, lam, epm_tol=epm_tol)

    PED_final = compute_PED(B_final, H, L, lam)
    colsums = PED_final.sum(axis=0)
    print(f"EPm after reduction:  {epm_after:.6f}  (delta={epm_after-epm_before:+.6f}, "
          f"{n_dropped} components dropped)")
    print(f"PED column sums after reduction: min={colsums.min():.6f}, max={colsums.max():.6f}")
    print(f"Final coefficient-matrix condition number: {np.linalg.cond(M_final):.4e}")

    # average nonzero components per final coordinate -- a concrete
    # "simplification achieved" number for the comparison report
    nnz_per_row = (np.abs(M_final) > ZERO_COEFF).sum(axis=1)
    print(f"Nonzero components per final coordinate: mean={nnz_per_row.mean():.2f}, "
          f"min={nnz_per_row.min()}, max={nnz_per_row.max()} "
          f"(before reduction, mean={(np.abs(M_mixed) > ZERO_COEFF).sum(axis=1).mean():.2f})")

    np.save(os.path.join(HERE, 'B_final.npy'), B_final)
    np.save(os.path.join(HERE, 'M_final.npy'), M_final)
    print("\nSaved: B_final.npy, M_final.npy")
