"""Unified mode classifier -- Algorithm 1 ``classify_all_modes``.

Implements the spec in IMPLEMENTATION_PLAN.md's "Authoritative spec" section
(includes the two-gate purity refinement) and the Steps 1-4 pseudocode in
``JCC_manuscript_structure_scoped.md`` Section 7.2. Consumes a
``(scorer, final)`` pair as produced by ``main.build_scorer_and_final``, so
classification runs on the identical mode pool ``score_modes()`` scores.
Self-contained; ``main.py`` orchestrates via ``main.run_classify_pipeline``.

Label vocabulary: a clean external (Step 3 gate pass) is the Step-2 slot name
("Tx".."Rz"); a mixed external+vibration (gate fail) is the slot name with a
trailing "*" (e.g. "Tx*"); internal modes (Step 4) are "S"/"B"/"SB"
(stretching/bending/mixed). Use ``is_external_label``/``external_axis``/
``is_clean_external``/``is_mixed_external``/``is_translation``/``is_rotation``
below rather than hand-rolling regex against these strings.

Pipeline
--------
Step 1  score every mode: {s[Tx..Tz], s[Rx..Rz], s[V_S]}, per-bond {s_AB}
        (s_AB is signed: positive = stretching, negative = compressing;
        s[V_S] itself is not signed).
Step 2  global external-mode assignment: one-to-one ``linear_sum_assignment``
        (scipy Hungarian solver) of the n_T+n_R external slots against ALL
        modes in the pool, maximizing sum |score|. Plain assignment only --
        no degenerate-axis-block special-casing (Decision 8: axis choice
        within a degenerate inertia tensor is a labeling convention fixed by
        the eigensolver, not an assignment ambiguity).
Step 3  two-gate purity test on each assigned (slot, mode) pair: clean iff
        |score_for_slot| >= tau_TR AND s[V_S] <= tau_B -> the bare slot name
        (e.g. "Tx"); else the slot name with a trailing "*" (e.g. "Tx*"),
        annotated with vib_label(s[V_S]).
Step 4  every mode NOT assigned an external slot in Step 2: vib_label(s[V_S])
        -> STRETCHING ("S") / BENDING ("B") / MIXED_STRETCH_BEND ("SB").
        Per-bond s_AB (signed) attached for STRETCHING / MIXED_STRETCH_BEND.

n_T = 3; n_R = 2 if linear else 3 (linear: smallest principal moment ~= 0).
"""

import json
import os
import re
from dataclasses import dataclass

import numpy as np
from scipy.optimize import linear_sum_assignment

from src.scoring import format_bond_map

# Phase-3 calibration output (src/calibrate.py); Thresholds() below is the
# explicit hardcoded fallback used if this file doesn't exist yet.
DEFAULT_CALIBRATION_PATH = os.path.join("data", "results", "thresholds.json")


# --- Classification labels ---------------------------------------------
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
    """Map a classify_all_modes() label to its semantic bucket:
    "translation"/"rotation"/"mixed_external" for external slots, or
    "stretch"/"bend"/"mixed" for Step-4 internal labels. Single shared home
    for this logic (replaces near-duplicate per-label dicts elsewhere);
    falls back to returning `label` unchanged."""
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

# Relative tolerance for the "smallest principal moment ~= 0" linear-molecule test.
LINEAR_TOL = 1e-6


@dataclass
class Thresholds:
    """Step 2/3/4 thresholds.

    tau_TR : purity bar for a clean external (near 1; PDF example 0.95).
    tau_S  : stretching bar on s[V_S] (>= -> STRETCHING).
    tau_B  : bending bar on s[V_S] (<= -> BENDING); also gate 2 of Step 3.

    Defaults (0.95/0.9/0.2) are the provisional pre-calibration constants,
    deliberately NOT auto-overwritten by Phase-3 calibration -- use
    `Thresholds.calibrated()` for the calibrated values instead.
    tests/test_classifier.py pins `Thresholds()` explicitly so its regression
    goldens stay fixed even if thresholds.json is later recalibrated;
    calibrated behavior has its own tests (tests/test_calibrate.py).
    """
    tau_TR: float = 0.95
    tau_S: float = 0.9
    tau_B: float = 0.2

    @classmethod
    def calibrated(cls, path=DEFAULT_CALIBRATION_PATH):
        """Load calibrated thresholds from `path` if it exists, else the class defaults."""
        if os.path.exists(path):
            with open(path) as f:
                data = json.load(f)
            return cls(tau_TR=data["tau_TR"], tau_S=data["tau_S"], tau_B=data["tau_B"])
        return cls()


def is_linear(scorer, tol=LINEAR_TOL):
    """True if the smallest principal moment of inertia is ~0 (linear molecule).
    Must be called post-MIT() (as build_scorer_and_final always does) -- see
    main.py::build_scorer_and_final for the n_R=2 / Rx-placeholder invariant this feeds."""
    moments, _ = scorer.principal_axes()
    scale = max(float(np.max(np.abs(moments))), 1e-12)
    return moments[0] <= tol * scale


def external_slots(scorer):
    """Return the (n_T, n_R, slot_labels) external slots per the spec:
    n_T=3 always, n_R=2 if linear else 3 (see main.py::build_scorer_and_final
    for why the excluded axis is always 'Rx')."""
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
    """Algorithm 1: score every mode, globally assign externals (Step 2:
    plain Hungarian assignment, no degenerate-axis special-casing -- Decision
    8, deliberately not reintroduced), apply two-gate purity (Step 3), then
    classify remaining internal modes (Step 4).

    scorer must already be MIT-aligned (as build_scorer_and_final leaves it);
    final is its candidate mode pool. thresholds defaults to
    Thresholds.calibrated(). Returns a list of dicts, one per mode in `final`
    (same order): {name, frequency, is_emit, T, R, V, classification,
    annotation, bonds, bonds_all}. 'bonds' (see ModeScorer.score_bonds()) is
    populated only for STRETCHING/MIXED_STRETCH_BEND, per spec; 'bonds_all'
    carries the same per-bond s_AB list unfiltered, for every mode regardless
    of classification (consumed by ped/merge_ped_scores.py's per-bond-type
    breakdown; classify_to_rows() ignores it, so the on-disk classified.csv
    schema is unchanged). 'classification' is the bare Step-2 slot name for a
    clean external, that slot name with a trailing "*" for mixed
    external+vibration, or "S"/"B"/"SB" for a Step-4 internal mode.
    'annotation' is "vibration=<vib_label>" for mixed-external modes, "" otherwise.
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
            # mu/k/irrep: only real Gaussian normal modes carry these; see
            # main.py::score_modes for the .get()-defaulting rationale.
            "reduced_mass": mode.get("reduced_mass"),
            "force_constant": mode.get("force_constant"),
            "irrep": mode.get("irrep"),
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

    # Attach per-bond s_AB to 'bonds' only for STRETCHING / MIXED_STRETCH_BEND,
    # per spec -- but keep the full unfiltered per-bond list on 'bonds_all'
    # for every mode (used by ped/merge_ped_scores.py's per-bond-type
    # breakdown, which needs s_AB regardless of a mode's S/B/SB label).
    for m in scored:
        m["bonds_all"] = m.pop("_bonds_all")
        m["bonds"] = m["bonds_all"] if m["classification"] in (STRETCHING, MIXED_STRETCH_BEND) else []

    return scored


def classify_to_rows(scored):
    """Flatten classify_all_modes() output into CSV-row dicts -- this is the
    single, complete per-mode result row (scores + Mu/K/Irrep + classification);
    there is no separate scores-only row shape. s_AB is a semicolon-joined
    'Elem#-Elem#:value' list, e.g. "C1-C2:0.0342". Mu/K/Irrep are None for
    EMIT modes and the synthetic ideal T/R references (only real Gaussian
    normal modes carry them)."""
    rows = []
    for m in scored:
        is_emit = m["is_emit"]
        bonds_str = format_bond_map(m["bonds"], "s_AB")
        rows.append({
            "Mode": m["name"],
            ("Eigenvalue" if is_emit else "Freq"): m["frequency"],
            "Tx": m["T"]["x"], "Ty": m["T"]["y"], "Tz": m["T"]["z"],
            "Rx": m["R"]["x"], "Ry": m["R"]["y"], "Rz": m["R"]["z"],
            "V_Stretch": m["V"],
            "Mu": m["reduced_mass"], "K": m["force_constant"], "Irrep": m["irrep"],
            "label": m["classification"],
            "annotation": m["annotation"],
            "s_AB": bonds_str,
        })
    return rows
