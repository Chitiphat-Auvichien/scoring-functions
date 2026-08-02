"""EMIT -> normal-mode projection (eq:emitproj): Theta_tilde = Q^T Theta.

Scores elsewhere (scoring.py) use unweighted Cartesian displacements; this
module is the one place mass-weighting enters, because a genuine orthonormal
reference basis requires it. Both Gaussian's normal-mode vectors and raw EMIT
eigenvectors are unit-length under the plain Cartesian inner product but only
mutually orthogonal under the mass-weighted inner product <u,v> = sum_A m_A
(u_A . v_A) (Eckart-Sayvetz). Verified on benzene's 30 real normal modes: the
unweighted Gram matrix has off-diagonals up to 0.80 (not orthogonal); the
mass-weighted Gram matrix's off-diagonals are <=4e-4 (numerical noise only).

Convention (locked): every mode vector v is mass-weighted as
v_mw[A] = sqrt(m_A) * v[A], then renormalized to unit length, before use in Q
or Theta. Stacking [ideal T, ideal R, real vibrational modes] this way gives
a near-orthonormal basis Q, so Theta_tilde = Q^T Theta (eq:emitproj) yields
per-column fractions Theta_tilde**2 summing to ~1 per EMIT mode (Parseval).

Internal/vibration fractions (C2_VS/C2_VB/C2_VMix): the stretch/bend/mixed
split of a real normal mode's contribution reuses Step-4's classification
(classifier.vib_label) on that mode's own s[V_S] score, so the boundary is
defined in exactly one place (classifier.py), not duplicated here.
"""

import numpy as np

from .classifier import Thresholds, vib_label, STRETCHING, BENDING, MIXED_STRETCH_BEND

# The 6 possible ideal external reference labels. A linear molecule's pool
# omits "Rx" (n_R=2); handled by simply never finding "Rx" in the labels.
EXTERNAL_LABELS = ("Tx", "Ty", "Tz", "Rx", "Ry", "Rz")

# Internal/vibration buckets a real normal mode's Theta_tilde**2 is summed
# into; imported from classifier.py rather than hardcoded to stay in sync.
_VIB_GROUPS = (STRETCHING, BENDING, MIXED_STRETCH_BEND)


def mass_weights_from_scorer(scorer):
    """sqrt(mass_A) per atom, repeated x3 (length 3N, atom order), read from
    the scorer's own masses array so masses stay a single source of truth."""
    return np.repeat(np.sqrt(scorer.masses), 3)


def _mass_weighted_unit_columns(mode_list, weights):
    """Flatten, mass-weight, and unit-renormalize each mode's (N,3) vector.
    Returns ndarray (3N, len(mode_list)), columns in mode_list order."""
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
        row["C2_VS"] = totals[STRETCHING]
        row["C2_VB"] = totals[BENDING]
        row["C2_VMix"] = totals[MIXED_STRETCH_BEND]

        total = sum(v for k, v in row.items() if k not in ("Mode", "Eigenvalue"))
        if abs(total - 1.0) > sum_tol:
            raise ValueError(
                f"{name}: projected fractions sum to {total:.4f}, expected ~1.0 "
                f"(Q may not be a near-orthonormal basis -- check inputs)."
            )
        rows.append(row)

        full_row = {"Mode": name, "Eigenvalue": eig}
        for lbl in labels:
            full_row[lbl] = float(Frac[idx_of[lbl], j])
        full_rows.append(full_row)

    return rows, full_rows
