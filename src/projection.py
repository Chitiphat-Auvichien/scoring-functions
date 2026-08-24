"""EMIT -> normal-mode projection (eq:emitproj): Theta_tilde = Q^T Theta.

Scores elsewhere (scoring.py) act on unweighted Cartesian displacements; this
module is the one place mass-weighted *coordinates* enter, because a genuine
orthonormal reference basis requires them. (Distinct from the V-score's own
per-bond reduced-mass weight under the 'mu' variant -- that scales whole-bond
contributions, it does not rescale the displacement vectors themselves; see
ModeScorer._bond_weights.) Both Gaussian's normal-mode vectors and raw EMIT
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

Second, explicitly non-orthonormal pathway (build_reference_basis_cartesian /
project_emit_cartesian, "Ocart_*" columns): kept in plain Cartesian (NOT
mass-weighted) coordinates deliberately, so it stays directly comparable to
the scores themselves (Tscore/Rscore/Vscore in scoring.py all act on raw,
unweighted Cartesian displacements) rather than to the mass-weighted C2_*
Parseval fractions. The basis Q is NOT mutually orthogonal (off-diagonals up
to 0.80), so "Ocart_Tx".."Ocart_Rz" are left as the raw signed overlap
<q_ref, theta_emit> itself (a cosine similarity between mode shapes,
unsquared -- squaring would only mean "fraction of character" if Q were
orthonormal, which it isn't) -- the same signed, unnormalized form as
Tscore/Rscore.

"Ocart_VS"/"Ocart_VB", by contrast, ARE normalized -- same mechanism Vscore
itself uses (magnitude-weighted total, not Parseval): each is the sum of
|overlap| over that group's reference modes, divided by the S+B total
|overlap| (Mix excluded from the denominator -- classify_all_modes()'s
default scheme="binary" never produces a MIXED_STRETCH_BEND reference mode,
so under normal use totals[Mix] is 0 anyway). This forces Ocart_VS +
Ocart_VB == 1 exactly and each term into [0,1], the same range as Vscore's
V_Stretch, without touching Q's geometry (no orthogonalization/rotation) --
consistent with Vscore's own non-orthogonal per-bond basis never being
orthogonalized either, just magnitude-normalized. "Ocart_VMix" is the same
group's |overlap| total divided by that same S+B-only denominator (so it is
NOT itself bounded by the VS/VB normalization -- it's 0 under the binary
scheme, and only nonzero at all if project_emit_cartesian is ever called
against a scheme="threeway" reference basis). "Ocart_Sum" reports the sum of
Ocart_Tx..Rz (raw signed) plus Ocart_VS/VB/VMix (now S+B-normalized) as a
diagnostic with no expected target value -- not a Parseval fraction, and not
asserted to sum to 1 (contrast C2_*'s own Sum, which is).
"""

import numpy as np

from .classifier import Thresholds, vib_label, STRETCHING, BENDING, MIXED_STRETCH_BEND
from .scoring import EPS_DENOM

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
    # v_weighting="*": this only buckets reference modes by vib_label(), it
    # never goes through classify_all_modes' guard, and the provisional
    # tau_S/tau_B here are deliberately definition-agnostic.
    thresholds = thresholds or Thresholds(v_weighting="*")
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


def _cartesian_unit_columns(mode_list):
    """Flatten and unit-renormalize each mode's raw Cartesian displacement
    vector -- no mass weighting (contrast _mass_weighted_unit_columns).
    Returns ndarray (3N, len(mode_list)), columns in mode_list order."""
    cols = []
    for mode in mode_list:
        v = np.asarray(mode["vector"], dtype=float).flatten()
        norm = np.linalg.norm(v)
        cols.append(v / norm if norm > 1e-12 else v)
    return np.array(cols).T


def build_reference_basis_cartesian(scorer, final_normal, thresholds=None):
    """Cartesian-overlap counterpart of build_reference_basis(): the same T/R/
    vibrational reference set, unit-normalized in the plain Cartesian inner
    product instead of the mass-weighted one. See module docstring -- this
    basis is NOT orthonormal in general, so it exists only as an explicit
    point of comparison, not a replacement for the mass-weighted pathway.

    Parameters / Returns: mirror build_reference_basis(), minus "weights"
    (there are none to reuse on the EMIT side under this convention).
    """
    thresholds = thresholds or Thresholds(v_weighting="*")
    Q = _cartesian_unit_columns(final_normal)

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
    return {"Q": Q, "labels": labels, "groups": groups}


def project_emit_cartesian(ref, final_emit):
    """Cartesian-overlap counterpart of project_emit() -- deliberately NOT
    squared. project_emit() squares Q.T @ Theta because Q is (near-)
    orthonormal there, so the squared coefficients have a Parseval
    ("fraction of mode character") interpretation. That basis property does
    not hold here (module docstring: off-diagonals up to 0.80), so squaring
    would not carry the same meaning -- "Ocart_Tx".."Ocart_Rz" instead report
    the raw signed overlap <q_ref, theta_emit> itself (each unit-normalized
    in plain Cartesian coordinates), same as a cosine similarity between the
    two mode shapes -- and the same signed, unnormalized form as Tscore/Rscore.

    "Ocart_VS"/"Ocart_VB" ARE normalized, though, via the same mechanism
    Vscore itself uses (magnitude-weighted total, not Parseval/squaring):
    |overlap| is summed per group, then divided by the S+B total |overlap|
    (Mix excluded from the denominator -- classify_all_modes()'s default
    scheme="binary" never produces a Mix reference mode, so totals[Mix] is 0
    under normal use). This forces Ocart_VS + Ocart_VB == 1 exactly, each
    in [0,1] -- the same range as V_Stretch, directly comparable to it,
    without cancellation and without touching Q's geometry. "Ocart_VMix" is
    that same group's |overlap| total divided by the same S+B-only
    denominator, so it is NOT covered by the VS/VB normalization (only
    nonzero if project_emit_cartesian is ever called against a
    scheme="threeway" reference basis). "Ocart_Sum" is the sum of all 9
    group columns (Tx..Rz raw signed + VS/VB/VMix now S+B-normalized),
    still just a diagnostic with no expected target value (contrast
    project_emit()'s Parseval-motivated sum~1 check).

    Parameters
    ----------
    ref : dict, output of build_reference_basis_cartesian().
    final_emit : list of mode dicts -- see project_emit().

    Returns
    -------
    (rows, full_rows) : same shape as project_emit(), with "Ocart_Tx".."Ocart_Rz"
        (raw signed overlaps) and "Ocart_VS"/"Ocart_VB"/"Ocart_VMix" (S+B-
        normalized fractions) plus "Ocart_Sum".
    """
    Q, labels, groups = ref["Q"], ref["labels"], ref["groups"]
    Theta = _cartesian_unit_columns(final_emit)

    Overlap = Q.T @ Theta  # raw signed overlap, NOT squared -- see docstring

    idx_of = {lbl: i for i, lbl in enumerate(labels)}

    rows = []
    full_rows = []
    for j, mode in enumerate(final_emit):
        name = mode.get("label", f"EMIT {j + 1}")
        eig = mode["frequency"]

        row = {"Mode": name, "Eigenvalue": eig}
        for slot in EXTERNAL_LABELS:
            row[f"Ocart_{slot}"] = float(Overlap[idx_of[slot], j]) if slot in idx_of else 0.0

        # VS/VB/VMix: |overlap| summed per group, normalized by the S+B
        # total only (Mix excluded from the denominator -- see docstring)
        # so Ocart_VS + Ocart_VB == 1 exactly, same [0,1] range as Vscore.
        totals = {g: 0.0 for g in _VIB_GROUPS}
        for lbl in labels:
            if lbl in EXTERNAL_LABELS:
                continue
            totals[groups[lbl]] += abs(Overlap[idx_of[lbl], j])
        denom = totals[STRETCHING] + totals[BENDING]
        if denom > EPS_DENOM:
            row["Ocart_VS"] = totals[STRETCHING] / denom
            row["Ocart_VB"] = totals[BENDING] / denom
            row["Ocart_VMix"] = totals[MIXED_STRETCH_BEND] / denom
        else:
            row["Ocart_VS"] = 0.0
            row["Ocart_VB"] = 0.0
            row["Ocart_VMix"] = 0.0
        row["Ocart_Sum"] = sum(v for k, v in row.items() if k not in ("Mode", "Eigenvalue"))
        rows.append(row)

        full_row = {"Mode": name, "Eigenvalue": eig}
        for lbl in labels:
            full_row[lbl] = float(Overlap[idx_of[lbl], j])
        full_rows.append(full_row)

    return rows, full_rows
