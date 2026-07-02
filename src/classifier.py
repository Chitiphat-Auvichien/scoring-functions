"""Unified mode classifier -- Algorithm 1 ``classify_all_modes``.

Implements the authoritative spec: IMPLEMENTATION_PLAN.md "Authoritative spec"
section (canonical, includes the two-gate purity refinement) and the
pseudocode in ``JCC_manuscript_structure_scoped.md`` Section 7.2 (structure of
Steps 1-4). Consumes a ``(scorer, final)`` pair of the same shape
``main.build_scorer_and_final`` produces, so classification always runs on the
identical mode pool that ``score_modes()`` scores (ideal T/R references
prepended for normal modes; raw EMIT eigenvectors only for EMIT). This module
is self-contained (no dependency on ``main.py``); ``main.py`` orchestrates by
calling into it -- see ``main.run_classify_pipeline``.

Pipeline
--------
Step 1  score every mode: {s[Tx..Tz], s[Rx..Rz], s[V_S]}, per-bond {s_AB}.
Step 2  global external-mode assignment: one-to-one ``linear_sum_assignment``
        (scipy Hungarian solver) of the n_T+n_R external slots against ALL
        modes in the pool, maximizing sum |score|. Plain assignment only --
        no degenerate-axis-block special-casing (retracted 2026-07-01,
        Decision 8: axis choice within a degenerate inertia tensor is a
        labeling convention fixed by the eigensolver, not an assignment
        ambiguity; normal-mode T/R references are built directly from
        geometry, never searched for).
Step 3  two-gate purity test on each assigned (slot, mode) pair: clean iff
        |score_for_slot| >= tau_TR AND s[V_S] <= tau_B -> CLEAN_TRANSLATION /
        CLEAN_ROTATION; else MIXED_EXTERNAL_WITH_VIBRATION, annotated with the
        dominant external slot and vib_label(s[V_S]).
Step 4  every mode NOT assigned an external slot in Step 2: vib_label(s[V_S])
        -> STRETCHING / BENDING / MIXED_STRETCH_BEND. Per-bond s_AB attached
        for STRETCHING / MIXED_STRETCH_BEND.

n_T = 3; n_R = 2 if linear else 3 (linear: smallest principal moment ~= 0).
"""

import json
import os
from dataclasses import dataclass

import numpy as np
from scipy.optimize import linear_sum_assignment

# Phase-3 calibration output (src/calibrate.py). Thresholds.calibrated() reads
# this if present; Thresholds() below remains the explicit, hardcoded
# provisional fallback for a fresh clone / before Phase 3 has been run.
DEFAULT_CALIBRATION_PATH = os.path.join("data", "results", "thresholds.json")


# --- Classification labels ---------------------------------------------
CLEAN_TRANSLATION = "CLEAN_TRANSLATION"
CLEAN_ROTATION = "CLEAN_ROTATION"
MIXED_EXTERNAL_WITH_VIBRATION = "MIXED_EXTERNAL_WITH_VIBRATION"
STRETCHING = "STRETCHING"
BENDING = "BENDING"
MIXED_STRETCH_BEND = "MIXED_STRETCH_BEND"

_T_SLOTS = ("Tx", "Ty", "Tz")
_R_SLOTS_FULL = ("Rx", "Ry", "Rz")

# Tolerance for the "smallest principal moment ~= 0" linear-molecule test,
# relative to the largest moment (mirrors DEGEN_TOL-style relative guards
# elsewhere in scoring.py; this is a separate, purpose-specific constant).
LINEAR_TOL = 1e-6


@dataclass
class Thresholds:
    """Step 2/3/4 thresholds.

    tau_TR : purity bar for a clean external (near 1; PDF example 0.95).
    tau_S  : stretching bar on s[V_S] (>= -> STRETCHING).
    tau_B  : bending bar on s[V_S] (<= -> BENDING); also gate 2 of Step 3.

    The defaults below (0.95 / 0.9 / 0.2) are the PROVISIONAL constants used
    before Phase-3 calibration existed -- kept here, unchanged, as the
    explicit hardcoded fallback (`Thresholds()`), not overwritten by
    calibration. Phase-3's calibrated values (derived in src/calibrate.py from
    the hydride library's literature stretch/bend labels + a tau_TR
    sensitivity sweep; see IMPLEMENTATION_PLAN.md) are loaded on demand via
    `Thresholds.calibrated()`, which classify_all_modes() now uses as its
    default when `thresholds=None` is passed and data/results/thresholds.json
    exists (falling back to these same provisional numbers if it does not --
    e.g. a fresh clone before Phase 3 has been run). This split is
    deliberate: tests/test_classifier.py pins `Thresholds()` explicitly, so
    those regression goldens stay fixed to these exact numbers even if a
    future recalibration (new library data) changes thresholds.json; the
    calibrated behavior has its own dedicated tests
    (tests/test_calibrate.py) that load Thresholds.calibrated() explicitly.
    """
    tau_TR: float = 0.95
    tau_S: float = 0.9
    tau_B: float = 0.2

    @classmethod
    def calibrated(cls, path=DEFAULT_CALIBRATION_PATH):
        """Load the Phase-3 calibrated thresholds from `path` if it exists;
        otherwise fall back to the provisional class defaults above (no
        warning -- running before Phase 3 has produced thresholds.json is an
        expected, supported state, not an error)."""
        if os.path.exists(path):
            with open(path) as f:
                data = json.load(f)
            return cls(tau_TR=data["tau_TR"], tau_S=data["tau_S"], tau_B=data["tau_B"])
        return cls()


def is_linear(scorer, tol=LINEAR_TOL):
    """True if the smallest principal moment of inertia is ~0 (linear molecule).

    Must be called on a scorer whose geometry reflects the frame the modes
    were scored in (i.e. after MIT() has rotated it into principal axes, as
    build_scorer_and_final always does) -- consistent with construct_R()'s
    n_R=2 special case for the on-axis rotation.
    """
    moments, _ = scorer.principal_axes()
    scale = max(float(np.max(np.abs(moments))), 1e-12)
    return moments[0] <= tol * scale


def external_slots(scorer):
    """Return the (n_T, n_R, slot_labels) external slots per the spec.

    n_T = 3 always. n_R = 2 if linear else 3. For a linear molecule, MIT()
    rotates the molecular (smallest-moment) axis onto the new X axis (its
    rotation matrix places the ascending-eigenvalue eigenvector in column 0),
    so the ill-defined on-axis rotation is always 'Rx' post-alignment; it is
    excluded, leaving Ry/Rz.
    """
    linear = is_linear(scorer)
    r_slots = ("Ry", "Rz") if linear else _R_SLOTS_FULL
    slots = _T_SLOTS + r_slots
    return 3, len(r_slots), slots


def vib_label(v, thresholds):
    """Step-4 internal sub-classification from s[V_S]."""
    if v >= thresholds.tau_S:
        return STRETCHING
    if v <= thresholds.tau_B:
        return BENDING
    return MIXED_STRETCH_BEND


def _score_slot(scores, slot):
    axis_type = "T" if slot[0] == "T" else "R"
    axis = slot[1].lower()
    return scores[axis_type][axis]


def classify_all_modes(scorer, final, thresholds=None):
    """Algorithm 1: score, globally assign externals, apply two-gate purity,
    then classify remaining internal modes.

    Parameters
    ----------
    scorer : ModeScorer
        Already constructed and MIT-aligned (as build_scorer_and_final leaves
        it) -- its current geometry is the frame 'final' vectors are in.
    final : list of mode dicts {frequency, vector, label?, is_emit?}
        The full candidate pool (ideal T/R + vibrational modes for 'normal';
        raw EMIT eigenvectors for 'emit') -- exactly what score_modes() scores.
    thresholds : Thresholds, optional. Defaults to Thresholds.calibrated() --
        the Phase-3 calibrated values if data/results/thresholds.json exists,
        else the same provisional constants as Thresholds().

    Returns
    -------
    list of dicts, one per mode in 'final' (same order), each:
        {name, frequency, is_emit, T, R, V, classification, annotation, bonds}
    'bonds' is a list of {'i','j','s_AB'} (0-based atom indices), populated
    only for STRETCHING / MIXED_STRETCH_BEND classifications; [] otherwise.
    'annotation' is "dominant_external=<slot>; vibration=<vib_label>" for
    MIXED_EXTERNAL_WITH_VIBRATION modes; "" otherwise.
    """
    thresholds = thresholds or Thresholds.calibrated()
    n_T, n_R, slots = external_slots(scorer)

    # ---- Step 1: score every mode (+ per-bond s_AB) ----
    scored = []
    for i, mode in enumerate(final):
        sc = scorer.calculate_scores(mode["vector"])
        bonds = scorer.score_bonds()  # reflects the just-loaded displacement
        scored.append({
            "name": mode.get("label", f"Mode {i+1}"),
            "frequency": mode["frequency"],
            "is_emit": mode.get("is_emit", False),
            "T": sc["T"],
            "R": sc["R"],
            "V": sc["V"],
            "_bonds_all": bonds,   # every bond's s_AB; filtered later by label
            "classification": None,
            "annotation": "",
        })

    n_modes = len(scored)

    # ---- Step 2: global external-mode assignment (plain one-to-one) ----
    cost = np.zeros((len(slots), n_modes))
    for si, slot in enumerate(slots):
        for mi, m in enumerate(scored):
            cost[si, mi] = abs(_score_slot(m, slot))

    try:
        row_ind, col_ind = linear_sum_assignment(cost, maximize=True)
    except TypeError:
        # Fallback for older scipy without the maximize kwarg: negate and minimize.
        row_ind, col_ind = linear_sum_assignment(-cost)

    assignment = {}  # mode index -> (slot, signed score value)
    for si, mi in zip(row_ind, col_ind):
        slot = slots[si]
        assignment[mi] = (slot, _score_slot(scored[mi], slot))

    # ---- Step 3: two-gate purity test ----
    for mi, (slot, score_value) in assignment.items():
        v = scored[mi]["V"]
        if abs(score_value) >= thresholds.tau_TR and v <= thresholds.tau_B:
            scored[mi]["classification"] = (
                CLEAN_TRANSLATION if slot in _T_SLOTS else CLEAN_ROTATION
            )
        else:
            scored[mi]["classification"] = MIXED_EXTERNAL_WITH_VIBRATION
            scored[mi]["annotation"] = (
                f"dominant_external={slot}; vibration={vib_label(v, thresholds)}"
            )

    # ---- Step 4: classify remaining (unassigned) internal modes ----
    for mi in range(n_modes):
        if scored[mi]["classification"] is None:
            scored[mi]["classification"] = vib_label(scored[mi]["V"], thresholds)

    # Attach per-bond s_AB only for STRETCHING / MIXED_STRETCH_BEND, per spec.
    for m in scored:
        if m["classification"] in (STRETCHING, MIXED_STRETCH_BEND):
            m["bonds"] = m.pop("_bonds_all")
        else:
            m.pop("_bonds_all")
            m["bonds"] = []

    return scored


def classify_to_rows(scored):
    """Flatten classify_all_modes() output into CSV-row dicts.

    Columns: Mode, Freq/Eigenvalue, Tx,Ty,Tz,Rx,Ry,Rz,V_Stretch, label,
    annotation, s_AB (1-based 'i-j:value' list, semicolon-joined; blank if
    bonds not attached for this mode's classification).
    """
    rows = []
    for m in scored:
        is_emit = m["is_emit"]
        bonds_str = ";".join(
            f"{b['i']+1}-{b['j']+1}:{b['s_AB']:.4f}" for b in m["bonds"]
        )
        rows.append({
            "Mode": m["name"],
            ("Eigenvalue" if is_emit else "Freq"): m["frequency"],
            "Tx": m["T"]["x"], "Ty": m["T"]["y"], "Tz": m["T"]["z"],
            "Rx": m["R"]["x"], "Ry": m["R"]["y"], "Rz": m["R"]["z"],
            "V_Stretch": m["V"],
            "label": m["classification"],
            "annotation": m["annotation"],
            "s_AB": bonds_str,
        })
    return rows
