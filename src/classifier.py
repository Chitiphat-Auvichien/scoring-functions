"""Unified mode classifier -- Algorithm 1 ``classify_all_modes``.

Implements the spec in IMPLEMENTATION_PLAN.md's "Authoritative spec" section
(includes the two-gate purity refinement) and the Steps 1-4 pseudocode in
``JCC_manuscript_structure_scoped.md`` Section 7.2. Consumes a
``(scorer, final)`` pair as produced by ``main.build_scorer_and_final``, so
classification runs on the identical mode pool ``score_modes()`` scores.
Self-contained; ``main.py`` orchestrates via ``main.run_classify_pipeline``.

Label vocabulary: a clean external (Step 3 gate pass) is the Step-2 slot name
("Tx".."Rz"); a mixed external+vibration (gate fail) is the slot name with a
trailing "*" (e.g. "Tx*"); internal modes (Step 4) are "S"/"B" under the
default BINARY scheme (``scheme="binary"``, single cutoff ``tau_SB``), which
forces every internal mode to STRETCHING or BENDING and never produces "SB"
-- this is the paper-standard scheme as of the 2026-08 binary-classification
decision. A second, THREE-WAY scheme (``scheme="threeway"``, explicit
opt-in) splits on ``tau_S``/``tau_B`` and may additionally produce "SB"
(mixed stretch/bend) -- kept fully functional for comparison/on-demand use,
just no longer the default. See ``vib_label``/``vib_label_binary`` below.
Use ``is_external_label``/``external_axis``/
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
        |score_for_slot| >= tau_TR AND s[V_S] <= gate2_bar -> the bare slot
        name (e.g. "Tx"); else the slot name with a trailing "*" (e.g.
        "Tx*"), annotated with vib_label(s[V_S]). gate2_bar is
        SCHEME-DEPENDENT: tau_purity (fixed, 0.05) under scheme="binary";
        tau_B (calibrated, ~0.17) under scheme="threeway", unchanged from
        before this split existed. Gate 1 (tau_TR) and Step 2's assignment
        itself remain scheme-independent.
Step 4  every mode NOT assigned an external slot in Step 2: vib_label(s[V_S])
        -> STRETCHING ("S") / BENDING ("B") / MIXED_STRETCH_BEND ("SB").
        Per-bond s_AB (signed) attached for every mode, regardless of label.

n_T = 3; n_R = 2 if linear else 3 (linear: smallest principal moment ~= 0).
"""

import json
import os
import re
from dataclasses import dataclass

import numpy as np
from scipy.optimize import linear_sum_assignment

from src.scoring import DEFAULT_V_WEIGHTING

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
    tau_S  : stretching bar on s[V_S] (>= -> STRETCHING). scheme="threeway" only.
    tau_B  : bending bar on s[V_S] (<= -> BENDING), and Step 3's gate 2 under
        scheme="threeway" ONLY (unchanged from before tau_purity existed).
        Under scheme="binary", gate 2 uses tau_purity instead -- see below.
    tau_SB : single-cutoff S/B split used by the *binary* classification
        scheme (scheme="binary"); unrelated to tau_S/tau_B's three-way split
        (a different threshold, not a synonym -- distinct name deliberately
        chosen to avoid collision). v >= tau_SB -> STRETCHING, else BENDING;
        never produces "SB".
    tau_purity : Step-3 gate 2's vibrational-leakage bound for a clean
        external, scheme="binary" only -- analogous to tau_TR's fixed
        purity bar, decoupled from tau_B/tau_S's three-way split (which
        conceptually doesn't exist under binary). Fixed default 0.05 ("allow
        5% vibrational leakage", symmetric to tau_TR's "allow 5%
        uncertainty" on the T/R-likeness side) -- a lighter-weight,
        exploratory constant like tau_SB's own initial rollout, NOT swept/
        calibrated here. scheme="threeway" continues to use tau_B for gate
        2, completely unchanged.
    v_weighting : which eq:vscore bond weighting these thresholds were
        calibrated against ('mu' or 'none'), or '*' to match any. tau_S and
        tau_B are read off an s[V_S] distribution, so they are only meaningful
        against the definition that produced it -- classify_all_modes() refuses
        a mismatch rather than silently mislabelling modes.

    Defaults (0.95/0.9/0.2/0.42/0.05) are the provisional pre-calibration
    constants, deliberately NOT auto-overwritten by Phase-3 calibration --
    use `Thresholds.calibrated()` for the calibrated values instead. tau_SB
    is a partial exception to that "provisional" framing: unlike tau_TR/
    tau_S/tau_B, no sweep ever computes a replacement for it (Phase-3
    calibration leaves it untouched -- see run_calibration_pipeline in
    src/calibrate.py), so THIS field default is the actual, sole source of
    truth for the canonical tau_SB value; there is no separate "calibrated"
    tau_SB to defer to. **2026-08-20: canonical default changed 0.50 -> 0.42**
    (the all-molecules classification-error-minimizing value from the
    tau_SB_error_sweep -- see IMPLEMENTATION_PLAN.md's Recent history). An
    earlier same-day attempt at this change only hand-edited
    data/results/thresholds.json's "tau_SB" key without changing this field;
    that was silently wiped by the next real `--calibrate` run because
    calibrate() (src/calibrate.py) never threads a tau_SB kwarg through to
    Thresholds(...) -- it relies entirely on this default. Fixed here (this
    field IS now 0.42) and in calibrate() (now writes "tau_SB" into the
    thresholds.json it produces, sourced from this same default, so the
    value is visible in the JSON without being a second, driftable source
    of truth).
    tests/test_classifier.py pins `Thresholds()` explicitly so its regression
    goldens stay fixed even if thresholds.json is later recalibrated;
    calibrated behavior has its own tests (tests/test_calibrate.py).
    """
    tau_TR: float = 0.95
    tau_S: float = 0.9
    tau_B: float = 0.2
    tau_SB: float = 0.42
    tau_purity: float = 0.05
    v_weighting: str = DEFAULT_V_WEIGHTING

    @classmethod
    def calibrated(cls, path=DEFAULT_CALIBRATION_PATH):
        """Load calibrated thresholds from `path` if it exists, else the class defaults."""
        if os.path.exists(path):
            with open(path) as f:
                data = json.load(f)
            # A thresholds.json written before the weighting variant existed
            # carries no stamp, and was by definition calibrated unweighted.
            # A thresholds.json written before tau_SB/tau_purity existed
            # carries no such key either -- fall back to the class default
            # (0.42 / 0.05 respectively) rather than KeyError. tau_purity is
            # a fixed exploratory constant, not swept by --calibrate, so it
            # is not expected to ever appear in thresholds.json; the
            # data.get() fallback is future-proofing, not the normal path.
            # tau_SB IS now written by every real --calibrate run (see
            # calibrate() in src/calibrate.py) -- the fallback here only
            # matters for a thresholds.json frozen before 2026-08-20.
            return cls(tau_TR=data["tau_TR"], tau_S=data["tau_S"], tau_B=data["tau_B"],
                       tau_SB=data.get("tau_SB", 0.42),
                       tau_purity=data.get("tau_purity", 0.05),
                       v_weighting=data.get("v_weighting", "none"))
        return cls()

    @classmethod
    def bootstrap(cls):
        """Provisional thresholds that match ANY weighting (sentinel '*').

        Needed for the chicken-and-egg first pass after switching variants:
        building library_scores.csv requires thresholds, but the calibrated
        ones are still stamped for the old definition. Safe because tau_S/tau_B
        derivation reads V_Stretch only and is itself threshold-independent --
        only the provisional pass's predicted_label column is affected, which
        is why --library is run again after --calibrate.
        """
        return cls(v_weighting="*")


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


def vib_label_binary(v, tau_SB):
    """Binary-scheme Step-4 internal sub-classification from s[V_S]: forces
    every mode to STRETCHING or BENDING via a single cutoff, never MIXED_STRETCH_BEND."""
    return STRETCHING if v >= tau_SB else BENDING


def vib_label(v, thresholds, scheme="binary"):
    """Step-4 internal sub-classification from s[V_S].

    scheme="binary" (default, paper-standard as of the 2026-08 binary-
    classification decision): delegates to vib_label_binary(v,
    thresholds.tau_SB) -- forces S or B, never SB. scheme="threeway"
    (explicit opt-in, kept fully functional for comparison): the three-way
    tau_S/tau_B split (may produce MIXED_STRETCH_BEND). One entry point
    keeps Step 3's mixed-external annotation and Step 4's label in
    agreement on which scheme is active.
    """
    if scheme == "binary":
        return vib_label_binary(v, thresholds.tau_SB)
    if v >= thresholds.tau_S:
        return STRETCHING
    if v <= thresholds.tau_B:
        return BENDING
    return MIXED_STRETCH_BEND


def gate2_bar(thresholds, scheme):
    """Step-3 gate 2's comparison bar for `scheme`: thresholds.tau_purity
    (fixed, 0.05) under scheme="binary"; thresholds.tau_B (calibrated) under
    scheme="threeway", unchanged from before this split existed. Single
    source of truth used by both classify_all_modes()'s own Step 3 and
    rescheme_external_label() below, so the two can never drift apart.
    """
    return thresholds.tau_purity if scheme == "binary" else thresholds.tau_B


def rescheme_external_label(predicted_label, score_value, v_stretch, thresholds, scheme):
    """Cheaply re-derive a mode's Step-3 clean-vs-mixed-external status under
    a DIFFERENT `scheme`/`thresholds` than the one it was originally
    classified with, from its own already-computed `score_value` (the
    assigned slot's own T/R score, e.g. a library_scores.csv row's "Tx"
    column value for a mode assigned the "Tx" slot) and `v_stretch`
    (s[V_S], Step 1's score) alone -- both scheme-independent Step-1/2
    outputs -- with NO re-assignment (Step 2) and no re-parsing. Mirrors
    `rescheme_internal_label`'s contract for Step 4.

    `predicted_label` must already be an external label (clean, e.g. "Tx",
    or mixed, e.g. "Tx*") for the recompute to apply; any other (Step-4
    internal) label is returned unchanged, since gate 2's threshold never
    touches Step 4. Gate 1 (tau_TR) is unaffected by `scheme` -- only gate
    2's bar (see `gate2_bar`) is scheme-dependent, so only rows that pass
    gate 1 can ever change clean<->mixed status here.
    """
    axis = external_axis(predicted_label)
    if axis is None:
        return predicted_label
    if abs(score_value) >= thresholds.tau_TR and v_stretch <= gate2_bar(thresholds, scheme):
        return axis
    return axis + "*"


def rescheme_internal_label(predicted_label, v_stretch, thresholds, scheme):
    """Cheaply re-derive a mode's Step-4 internal label under a DIFFERENT
    `scheme` than the one it was originally classified with, from its own
    `v_stretch` (s[V_S], Step 1's score -- always scheme-independent) alone,
    with NO re-run of Step 2/3 (external assignment/purity, also scheme-
    independent). `predicted_label` is the mode's already-computed
    classify_all_modes() label (under whichever scheme produced it, e.g. a
    library_scores.csv row); if its bucket is a vibrational one
    (stretch/bend/mixed) it is recomputed via vib_label(v_stretch,
    thresholds, scheme); a Step-2 external-slot label (clean or
    mixed-external, e.g. "Tx"/"Tx*") is returned unchanged, since scheme
    never touches Step 2/3.

    This is the general form of the per-row re-labeling pattern already used
    ad hoc in a few places (e.g. figures.py's plot_confusion_matrix_binary/
    plot_transferability_confusion_binary/plot_benzene_internal_confusion_binary's
    own inline `_binarize` closures, which special-case scheme="binary"
    specifically) -- lets a caller obtain ANY scheme's labels from an
    already-scored table without re-parsing/re-running the whole classifier
    (e.g. src.benzene_validation's threeway-only SB diagnostics, needed once
    library_scores.csv's own predicted_label became binary-by-default).
    """
    bucket = classification_bucket(predicted_label)
    if bucket in ("stretch", "bend", "mixed"):
        return vib_label(v_stretch, thresholds, scheme=scheme)
    return predicted_label


def _score_slot(scores, slot):
    axis_type = "T" if slot[0] == "T" else "R"
    axis = slot[1].lower()
    return scores[axis_type][axis]


def _assert_weighting_match(scorer, thresholds):
    """Refuse to classify when the thresholds were calibrated against a
    different V-score definition than the scorer is producing.

    tau_S/tau_B are cut points on an s[V_S] distribution, so pairing them with
    the other definition shifts every stretch/bend boundary silently -- the
    modes would still be labelled, just wrongly. Every pipeline funnels through
    classify_all_modes(), so this one check covers them all.
    """
    if thresholds.v_weighting in ("*", scorer.v_weighting):
        return
    raise ValueError(
        f"V-score definition mismatch: thresholds were calibrated with "
        f"v_weighting={thresholds.v_weighting!r} but this scorer is using "
        f"{scorer.v_weighting!r}. tau_S/tau_B are cut points on the s[V_S] "
        f"distribution, so mixing the two silently mislabels stretches and "
        f"bends. Either re-run `python main.py --library --calibrate` under "
        f"the current weighting, or score with "
        f"`--v-weighting {thresholds.v_weighting}`.")


def classify_all_modes(scorer, final, thresholds=None, scheme="binary"):
    """Algorithm 1: score every mode, globally assign externals (Step 2:
    plain Hungarian assignment, no degenerate-axis special-casing -- Decision
    8, deliberately not reintroduced), apply two-gate purity (Step 3), then
    classify remaining internal modes (Step 4).

    scorer must already be MIT-aligned (as build_scorer_and_final leaves it);
    final is its candidate mode pool. thresholds defaults to
    Thresholds.calibrated(). `scheme` is "binary" (default, paper-standard:
    single tau_SB cutoff, never produces "SB") or "threeway" (explicit
    opt-in, may produce "SB"; kept fully functional for comparison/on-demand
    use) -- passed through to every vib_label() call site below (Step 3's
    mixed-external annotation AND Step 4's internal label), so the two stay
    in agreement. It does NOT affect Step 2's assignment itself, or Step 3's
    gate 1 (tau_TR) -- but it DOES select Step 3's gate 2 bar (see
    `gate2_bar`): tau_purity under "binary", tau_B under "threeway"
    (unchanged historical behavior). Returns a list of dicts, one
    per mode in `final`
    (same order): {name, frequency, is_emit, T, R, V, classification,
    annotation, bonds, bonds_all}. 'bonds' (see ModeScorer.score_bonds()) carries
    the per-bond s_AB list for every mode regardless of classification; 'bonds_all'
    is the same list under the name ped/merge_ped_scores.py's per-bond-type
    breakdown reads. 'classification' is the bare Step-2 slot name for a
    clean external, that slot name with a trailing "*" for mixed
    external+vibration, or "S"/"B"/"SB" (or just "S"/"B" under scheme="binary")
    for a Step-4 internal mode.
    'annotation' is "vibration=<vib_label>" for mixed-external modes, "" otherwise.
    """
    if scheme not in ("threeway", "binary"):
        raise ValueError(f"scheme must be 'threeway' or 'binary', got {scheme!r}")
    thresholds = thresholds or Thresholds.calibrated()
    _assert_weighting_match(scorer, thresholds)
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
    # gate2_bar is scheme-dependent: tau_purity (fixed) under "binary",
    # tau_B (calibrated) under "threeway" -- see gate2_bar()'s docstring.
    gate2 = gate2_bar(thresholds, scheme)
    for mi, (slot, score_value) in assignment.items():
        v = scored[mi]["V"]
        if abs(score_value) >= thresholds.tau_TR and v <= gate2:
            scored[mi]["classification"] = slot
        else:
            scored[mi]["classification"] = slot + "*"
            # The axis is already encoded in the classification string above
            # (e.g. "Tx*"), so the annotation only adds the one piece of
            # information it doesn't already carry: the vibration sub-label.
            scored[mi]["annotation"] = f"vibration={vib_label(v, thresholds, scheme)}"

    # ---- Step 4: classify remaining (unassigned) internal modes ----
    for mi in range(n_modes):
        if scored[mi]["classification"] is None:
            scored[mi]["classification"] = vib_label(scored[mi]["V"], thresholds, scheme)

    # 'bonds' carries the per-bond s_AB list for every mode regardless of
    # classification; 'bonds_all' is kept as an alias (same list) for
    # ped/merge_ped_scores.py, which reads it by that name.
    for m in scored:
        m["bonds_all"] = m.pop("_bonds_all")
        m["bonds"] = m["bonds_all"]

    return scored


def classify_to_rows(scored):
    """Flatten classify_all_modes() output into CSV-row dicts -- this is the
    single, complete per-mode result row (scores + Mu/K/Irrep + classification);
    there is no separate scores-only row shape. Each bond gets its own
    's_AB[Elem#-Elem#]' column (e.g. "s_AB[C1-C2]") rather than one
    semicolon-joined string, so the CSV is a plain rectangular table --
    every row lists the same molecule's bonds in the same order (see
    ModeScorer.bList), so the column set is identical across rows and easy
    to sort/filter/plot in Excel. Mu/K/Irrep are None for EMIT modes and the
    synthetic ideal T/R references (only real Gaussian normal modes carry
    them)."""
    rows = []
    for m in scored:
        is_emit = m["is_emit"]
        row = {
            "Mode": m["name"],
            ("Eigenvalue" if is_emit else "Freq"): m["frequency"],
            "Tx": m["T"]["x"], "Ty": m["T"]["y"], "Tz": m["T"]["z"],
            "Rx": m["R"]["x"], "Ry": m["R"]["y"], "Rz": m["R"]["z"],
            "V_Stretch": m["V"],
            "Mu": m["reduced_mass"], "K": m["force_constant"], "Irrep": m["irrep"],
            "label": m["classification"],
            "annotation": m["annotation"],
        }
        for b in m["bonds"]:
            row[f"s_AB[{b['i_label']}-{b['j_label']}]"] = b["s_AB"]
        rows.append(row)
    return rows
