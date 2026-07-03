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

Label vocabulary (author-approved rename, 2026-07-02; short, axis-specific
symbolic scheme, replacing the earlier CLEAN_TRANSLATION/CLEAN_ROTATION/
MIXED_EXTERNAL_WITH_VIBRATION/STRETCHING/BENDING/MIXED_STRETCH_BEND strings)
--------------------------------------------------------------------------
Clean external (Step 3, gate pass)   -> the specific Step-2 slot name:
    "Tx", "Ty", "Tz", "Rx", "Ry", "Rz".
Mixed external + vibration (Step 3, gate fail) -> the same slot name with a
    trailing "*": "Tx*", "Ty*", "Tz*", "Rx*", "Ry*", "Rz*".
Stretching / bending / mixed stretch-bend (Step 4) -> "S" / "B" / "SB"
    (Python constant names STRETCHING/BENDING/MIXED_STRETCH_BEND unchanged,
    only their string VALUES changed, to minimize import-site churn).
See ``is_external_label``/``external_axis``/``is_clean_external``/
``is_mixed_external``/``is_translation``/``is_rotation`` below for the
reusable predicates downstream code should use instead of hand-rolling regex
against these strings.

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
        |score_for_slot| >= tau_TR AND s[V_S] <= tau_B -> the bare slot name
        (e.g. "Tx"); else the slot name with a trailing "*" (e.g. "Tx*"),
        annotated with vib_label(s[V_S]).
Step 4  every mode NOT assigned an external slot in Step 2: vib_label(s[V_S])
        -> STRETCHING ("S") / BENDING ("B") / MIXED_STRETCH_BEND ("SB").
        Per-bond s_AB attached for STRETCHING / MIXED_STRETCH_BEND.

n_T = 3; n_R = 2 if linear else 3 (linear: smallest principal moment ~= 0).
"""

import json
import os
import re
from dataclasses import dataclass

import numpy as np
from scipy.optimize import linear_sum_assignment

# Phase-3 calibration output (src/calibrate.py). Thresholds.calibrated() reads
# this if present; Thresholds() below remains the explicit, hardcoded
# provisional fallback for a fresh clone / before Phase 3 has been run.
DEFAULT_CALIBRATION_PATH = os.path.join("data", "results", "thresholds.json")


# --- Classification labels ---------------------------------------------
# Clean-external / mixed-external labels are no longer fixed constants -- they
# are axis-specific strings ("Tx".."Rz", optionally with a trailing "*")
# assigned dynamically in classify_all_modes() (Step 3) from whichever slot
# Step 2 assigned. See the module docstring above and the predicates below.
STRETCHING = "S"
BENDING = "B"
MIXED_STRETCH_BEND = "SB"

_T_SLOTS = ("Tx", "Ty", "Tz")
_R_SLOTS_FULL = ("Rx", "Ry", "Rz")

_EXTERNAL_LABEL_RE = re.compile(r"^[TR][xyz]\*?$")


def is_external_label(label):
    """True if `label` matches the external-slot pattern ^[TR][xyz]\\*?$
    (clean, e.g. "Tx", or mixed, e.g. "Tx*")."""
    return isinstance(label, str) and bool(_EXTERNAL_LABEL_RE.match(label))


def external_axis(label):
    """The bare slot name for an external label ("Tx*" -> "Tx"; "Tx" ->
    "Tx"), or None if `label` is not an external label at all."""
    if not is_external_label(label):
        return None
    return label[:-1] if label.endswith("*") else label


def is_clean_external(label):
    """True iff `label` is an external label with NO trailing "*" (passed
    both Step-3 purity gates)."""
    return is_external_label(label) and not label.endswith("*")


def is_mixed_external(label):
    """True iff `label` is an external label WITH a trailing "*" (failed at
    least one Step-3 purity gate)."""
    return is_external_label(label) and label.endswith("*")


def is_translation(label):
    """True iff `label` is an external label (clean or mixed) whose slot is
    a translation (Tx/Ty/Tz)."""
    axis = external_axis(label)
    return axis is not None and axis[0] == "T"


def is_rotation(label):
    """True iff `label` is an external label (clean or mixed) whose slot is
    a rotation (Rx/Ry/Rz)."""
    axis = external_axis(label)
    return axis is not None and axis[0] == "R"


def classification_bucket(label):
    """Map any classify_all_modes() classification label to its semantic
    bucket name -- "translation"/"rotation"/"mixed_external" for external
    slots (clean vs. mixed distinguished by the trailing "*"), or
    "stretch"/"bend"/"mixed" for the Step-4 internal vibration labels.
    Falls back to returning `label` unchanged for anything else (e.g. a
    string that is already a bucket name), mirroring the ``dict.get(x, x)``
    fallback pattern used by the old per-label dict lookups this replaces
    (src/calibrate.py's ``_BUCKET``, src/benzene_validation.py's
    ``_PRED_TO_BUCKET``) -- a single shared home for that logic instead of
    two near-duplicate dicts that would otherwise need 6 axis-specific keys
    each for clean and 6 more for mixed.
    """
    if is_mixed_external(label):
        return "mixed_external"
    if is_translation(label):
        return "translation"
    if is_rotation(label):
        return "rotation"
    if label == STRETCHING:
        return "stretch"
    if label == BENDING:
        return "bend"
    if label == MIXED_STRETCH_BEND:
        return "mixed"
    return label

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
    'bonds' is a list of {'i','j','s_AB','rel_db','i_label','j_label'} (see
    ModeScorer.score_bonds()), populated only for STRETCHING /
    MIXED_STRETCH_BEND classifications; [] otherwise.
    'classification' is the bare Step-2 slot name ("Tx".."Rz") for a clean
    external, that same slot name with a trailing "*" (e.g. "Tx*") for a
    mixed external+vibration mode, or "S"/"B"/"SB" (STRETCHING/BENDING/
    MIXED_STRETCH_BEND) for a Step-4 internal mode.
    'annotation' is "vibration=<vib_label>" for mixed-external ("*"-suffixed)
    modes -- the axis is intentionally NOT repeated here since the top-level
    classification string already names it; "" otherwise.
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
            scored[mi]["classification"] = slot
        else:
            scored[mi]["classification"] = slot + "*"
            # The axis is already encoded in the classification string above
            # (e.g. "Tx*"), so the annotation only adds the one piece of
            # information it doesn't already carry: the vibration sub-label.
            scored[mi]["annotation"] = f"vibration={vib_label(v, thresholds)}"

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
    annotation, s_AB (atom-symbol 'Elem#-Elem#:value' list, e.g.
    "C1-C2:0.0342", semicolon-joined; blank if bonds not attached for this
    mode's classification).
    """
    rows = []
    for m in scored:
        is_emit = m["is_emit"]
        bonds_str = ";".join(
            f"{b['i_label']}-{b['j_label']}:{b['s_AB']:.4f}" for b in m["bonds"]
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
