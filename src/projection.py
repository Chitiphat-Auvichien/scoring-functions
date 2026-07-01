"""EMIT -> normal-mode projection (eq:emitproj): Theta_tilde = Q^T Theta.

Phase-0 locked decision (mass-weighting convention, IMPLEMENTATION_PLAN.md):
------------------------------------------------------------------------------
Scores (s[T]/s[R]/s[V_S], scoring.py) are computed from UNWEIGHTED Cartesian
displacements -- that convention is unchanged and is NOT touched here. This
module is the one place mass-weighting enters, because it is what the
projection reference needs to be a genuine orthonormal basis.

Both Gaussian's printed normal-mode vectors and the raw EMIT eigenvectors are
normalized to unit length under the plain (unweighted) Cartesian inner
product, but they are only mutually ORTHOGONAL under the mass-weighted inner
product <u, v> = sum_A m_A (u_A . v_A) -- this is verified numerically on
benzene: the plain Gram matrix of the 30 real normal modes has off-diagonal
entries up to 0.80, while the mass-weighted Gram matrix's off-diagonals are
<=4e-4 (consistent with Gaussian's ~5-6 significant-figure print precision,
i.e. numerical noise, not a real deviation from orthogonality). The same
holds for the raw EMIT eigenvectors read from data/EMIT/<mol>_EMIT.txt. This
is the standard Eckart-Sayvetz signature: true harmonic normal modes are
orthonormal in mass-weighted coordinates; Gaussian (and, empirically, EMIT)
report them back-transformed to unweighted Cartesian and renormalized to unit
Euclidean length for display.

Convention (LOCKED): for every reference/candidate mode vector v (ideal T/R
built the same geometric way as ModeScorer.construct_T/construct_R, or a real
vibrational normal mode, or a raw EMIT eigenvector), define the mass-weighted
form v_mw[A] = sqrt(m_A) * v[A] (same scalar applied to all 3 Cartesian
components of atom A), then renormalize v_mw to unit Euclidean length. Given
COM + principal-axis alignment (already enforced by ModeScorer.COM/MIT), the
ideal T/R block is analytically exactly orthogonal to itself in this metric,
and the Eckart-Sayvetz theorem makes it (numerically-exactly) orthogonal to
the true vibrational normal modes -- so stacking [ideal T, ideal R, real
vibrational modes] gives a basis Q that is orthonormal to within the input
files' print precision (~1e-4, verified above). Projecting the
similarly-mass-weighted-and-renormalized EMIT eigenvector Theta onto Q,

    Theta_tilde = Q^T Theta                                    (eq:emitproj)

then yields per-column fractional contributions Theta_tilde**2 that sum to
~1 per EMIT mode (Parseval), matching the semantics of the existing
data/results/benzene_EMIT_contributions.csv (validated to ~1e-4 absolute
agreement on EMIT 2/9/34/35/36 against that file during development).

Internal/vibration fractions (C2_VS / C2_VB / C2_VMix): the manuscript's
external split (translation/rotation/vibration fractions) falls straight out
of Q's 6 ideal-T/R columns. The stretch/bend/mixed split of the *vibration*
fraction is NOT itself a projection quantity -- it is obtained by reusing the
Step-4 internal classification (classifier.vib_label) on each REAL
vibrational normal mode's own (unweighted) s[V_S] score, then summing that
normal mode's Theta_tilde**2 into the corresponding bucket. This reproduces
the ground-truth CSV's C2_VS/C2_VB/C2_VMix columns exactly and keeps the
stretch/bend boundary defined in exactly one place (classifier.py), not
duplicated here.
"""

import numpy as np

from .classifier import Thresholds, vib_label

# The 6 possible ideal external reference labels. A linear molecule's pool
# (main.build_scorer_and_final) omits "Rx" (n_R=2); this module handles that
# by simply never finding "Rx" in the reference labels -- no special-casing.
EXTERNAL_LABELS = ("Tx", "Ty", "Tz", "Rx", "Ry", "Rz")

# Internal/vibration buckets a real normal mode's Theta_tilde**2 is summed
# into (mirrors classifier.py's Step-4 labels).
_VIB_GROUPS = ("STRETCHING", "BENDING", "MIXED_STRETCH_BEND")


def mass_weights_from_scorer(scorer):
    """sqrt(mass_A) per atom, repeated x3 (one weight per Cartesian
    component) -- length 3N, in atom order. Reads masses straight off the
    scorer's Atom objects (single source of truth; matches whatever
    atomicMass lookup ModeScorer already did), so this module never
    re-derives masses from symbols independently.
    """
    masses = np.array([atom.rMass for atom in scorer.atoms], dtype=float)
    return np.repeat(np.sqrt(masses), 3)


def _mass_weighted_unit_columns(mode_list, weights):
    """Flatten each mode's (N,3) vector, mass-weight it, and renormalize to
    unit Euclidean length. Returns an ndarray (3N, len(mode_list)) whose
    columns are in the same order as mode_list.
    """
    cols = []
    for mode in mode_list:
        v = np.asarray(mode["vector"], dtype=float).flatten() * weights
        norm = np.linalg.norm(v)
        cols.append(v / norm if norm > 1e-12 else v)
    return np.array(cols).T


def build_reference_basis(scorer, final_normal, thresholds=None):
    """Build the mass-weighted normal-mode reference basis Q.

    Parameters
    ----------
    scorer : ModeScorer
        Already MIT-aligned (as main.build_scorer_and_final(raw, "normal")
        leaves it) -- its current geometry/atoms are the frame final_normal's
        vectors are in, and the source of the per-atom masses.
    final_normal : list of mode dicts
        The 'normal' candidate pool: 3 ideal T + (2 or 3) ideal R + the real
        3N-6 vibrational normal modes (exactly what main.build_scorer_and_final
        (raw, "normal") returns) -- a full 3N-dimensional set.
    thresholds : classifier.Thresholds, optional (defaults to the provisional
        constants also used by classify_all_modes, so the internal
        stretch/bend/mixed split stays consistent with Step 4 elsewhere).

    Returns
    -------
    dict with:
      "Q"       : ndarray (3N, 3N), mass-weighted unit-normalized reference
                  columns, in final_normal's order.
      "labels"  : list of str, the mode labels (Tx..Rz, Vib 1..Vib (3N-6)).
      "groups"  : dict label -> "EXTERNAL" | "STRETCHING" | "BENDING" |
                  "MIXED_STRETCH_BEND".
      "weights" : the mass_weights_from_scorer(scorer) vector (reused so the
                  EMIT side mass-weights identically).
    """
    thresholds = thresholds or Thresholds()
    weights = mass_weights_from_scorer(scorer)
    Q = _mass_weighted_unit_columns(final_normal, weights)

    labels = []
    groups = {}
    for mode in final_normal:
        label = mode.get("label")
        labels.append(label)
        if label in EXTERNAL_LABELS:
            groups[label] = "EXTERNAL"
        else:
            sc = scorer.calculate_scores(mode["vector"])
            groups[label] = vib_label(sc["V"], thresholds)
    return {"Q": Q, "labels": labels, "groups": groups, "weights": weights}


def project_emit(ref, final_emit, sum_tol=0.05):
    """Project EMIT eigenvectors onto the reference basis (eq:emitproj).

    Parameters
    ----------
    ref : dict, output of build_reference_basis().
    final_emit : list of mode dicts -- the raw EMIT candidate pool (exactly
        main.build_scorer_and_final(raw, "emit")'s output: all 3N raw EMIT
        eigenvectors, unrotated, already assumed to be in the same
        principal-axis frame the 'normal' reference geometry was rotated
        into -- see main.build_scorer_and_final's docstring).
    sum_tol : float
        Soft sanity check: each EMIT mode's total fractional contribution
        (sum of C2_* columns) should be ~1 (Parseval, since Q is
        near-orthonormal). Raises if any mode deviates by more than this
        (catches a gross convention/shape bug; small ~1e-4 numerical
        residual from Q's imperfect orthonormality is expected and NOT an
        error).

    Returns
    -------
    (rows, full_rows) : two lists of dicts (row-per-EMIT-mode).
      rows      : {"Mode", "Eigenvalue", "C2_Tx".."C2_Rz", "C2_VS", "C2_VB",
                   "C2_VMix"} -- matches the columns/semantics of the existing
                   data/results/benzene_EMIT_contributions.csv ground truth.
      full_rows : {"Mode", "Eigenvalue", <one column per reference label>}
                  -- the per-individual-normal-mode Theta_tilde**2 detail (the
                  richer "projection-coefficients data file" Phase 2 calls
                  for, not just the grouped external/internal fractions).
    """
    Q, labels, groups, weights = ref["Q"], ref["labels"], ref["groups"], ref["weights"]
    Theta = _mass_weighted_unit_columns(final_emit, weights)

    Proj = Q.T @ Theta          # (n_ref, n_emit)
    Frac = Proj ** 2            # fractional contributions (eq:emitproj coefficients, squared)

    idx_of = {lbl: i for i, lbl in enumerate(labels)}

    rows = []
    full_rows = []
    for j, mode in enumerate(final_emit):
        name = mode.get("label", f"EMIT {j + 1}")
        eig = mode["frequency"]

        row = {"Mode": name, "Eigenvalue": eig}
        for slot in EXTERNAL_LABELS:
            row[f"C2_{slot}"] = float(Frac[idx_of[slot], j]) if slot in idx_of else 0.0

        totals = {g: 0.0 for g in _VIB_GROUPS}
        for lbl in labels:
            if lbl in EXTERNAL_LABELS:
                continue
            totals[groups[lbl]] += Frac[idx_of[lbl], j]
        row["C2_VS"] = totals["STRETCHING"]
        row["C2_VB"] = totals["BENDING"]
        row["C2_VMix"] = totals["MIXED_STRETCH_BEND"]

        total = sum(v for k, v in row.items() if k not in ("Mode", "Eigenvalue"))
        assert abs(total - 1.0) <= sum_tol, (
            f"{name}: projected fractions sum to {total:.4f}, expected ~1.0 "
            f"(Q may not be a near-orthonormal basis -- check inputs)."
        )
        rows.append(row)

        full_row = {"Mode": name, "Eigenvalue": eig}
        for lbl in labels:
            full_row[lbl] = float(Frac[idx_of[lbl], j])
        full_rows.append(full_row)

    return rows, full_rows
