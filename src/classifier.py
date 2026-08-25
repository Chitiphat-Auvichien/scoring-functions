"""Unified mode classifier -- Algorithm 1 ``classify_all_modes``.

Implements the flowchart referenced by the JCC manuscript's ``fig:flowchart``
(``JCC/JCC_man_scoring/images/classifier_algorithm.pdf``), restructured
2026-08-25 to separate two previously-entangled questions:

Step 1  score every mode: {s[Tx..Tz], s[Rx..Rz], s[V_S]}, per-bond {s_AB}
        (s_AB is signed: positive = stretching, negative = compressing;
        s[V_S] itself is not signed).
Step 2  UNCONDITIONAL, every one of the 3N candidate modes (including the
        synthetic ideal T/R reference rows -- their V_Stretch=0, so they
        trivially get vib_label="B" under the binary scheme; expected, not a
        bug, since downstream consumers must treat tr_label as authoritative
        for T/R identity, never vib_label): a Step-4-style vibrational label
        ("S"/"B" under the default binary scheme, or "S"/"B"/"SB" under
        scheme="threeway") from s[V_S] alone, via ``vib_label()``. Never
        depends on whether the mode is a good T/R match.
Step 3  OPTIONAL (``identify_tr=True`` by default), independent of Step 2:
        find the best n_T+n_R modes for the T/R slots via the same one-to-one
        ``linear_sum_assignment(cost, maximize=True)`` Hungarian solver as
        before (Decision 8: plain assignment, no degenerate-axis-block
        special-casing), but with NO purity gate at all. Every winner gets
        the slot name (Tx..Rz) and the assignment's own raw signed score --
        no asterisk, no "impure" variant, no annotation string. If
        ``identify_tr=False``, every mode's tr_label/tr_score stay None.

n_T = 3; n_R = 2 if linear else 3 (linear: smallest principal moment ~= 0).

RETIRED 2026-08-25 (do not reintroduce without re-reading
IMPLEMENTATION_PLAN.md's dated entry): the old two-gate purity test
(``gate2_bar``, ``Thresholds.tau_purity``, the starred "Tx*" mixed-external
label + its "vibration=<label>" annotation, ``rescheme_external_label``).
That design made Step 3 (T/R identity) depend on Step 2's vibrational
character, which was almost always true (hence "impure") for non-normal-mode
input (e.g. EMIT) since there is no real translation/rotation in that basis
at all -- the flag was nearly useless outside the normal-mode case. The
manuscript prose describing the old two-gate behavior (including the benzene
EMIT 34-36 "flagged with contaminating SB character" narrative) has NOT yet
been updated to match this file -- a separate ``/revise-section`` pass is
needed once this code change lands (intentionally out of scope here).
``Thresholds.tau_TR`` itself is NOT deleted -- it survives as a diagnostic-
only lens (e.g. for a human to eyeball "how confident is this T/R
assignment"), just no longer used to gate anything in the pipeline.

Label vocabulary: ``vib_label`` (Step 2, every mode) is one of "S"/"B" under
the default BINARY scheme (``scheme="binary"``, single cutoff ``tau_SB``),
which forces every mode to STRETCHING or BENDING and never produces "SB" --
this is the paper-standard scheme as of the 2026-08 binary-classification
decision. A second, THREE-WAY scheme (``scheme="threeway"``, explicit
opt-in) splits on ``tau_S``/``tau_B`` and may additionally produce "SB"
(mixed stretch/bend) -- kept fully functional for comparison/on-demand use,
just no longer the default. ``tr_label`` (Step 3, optional) is either None
(no external slot won, or ``identify_tr=False``) or one of the 6 axis-
specific slot names ("Tx".."Rz") -- never starred. Use
``is_external_label``/``external_axis``/``is_clean_external``/
``is_mixed_external``/``is_translation``/``is_rotation`` below rather than
hand-rolling regex against a label string; these are unchanged generic
utilities -- ``is_mixed_external`` will simply never match a
freshly-produced ``tr_label`` now (no more asterisk), which is correct, not
a regression.
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
    (clean, e.g. "Tx", or -- from legacy pre-2026-08-25 data only -- mixed,
    e.g. "Tx*"). A freshly-produced ``tr_label`` is never starred."""
    return isinstance(label, str) and bool(_EXTERNAL_LABEL_RE.match(label))


def external_axis(label):
    """The bare slot name for an external label ("Tx*" -> "Tx"; "Tx" ->
    "Tx"), or None if `label` is not an external label at all."""
    if not is_external_label(label):
        return None
    return label[:-1] if label.endswith("*") else label


def is_clean_external(label):
    """True iff `label` is an external label with NO trailing "*". Every
    freshly-produced ``tr_label`` satisfies this trivially (Step 3 no longer
    produces a starred variant at all)."""
    return is_external_label(label) and not label.endswith("*")


def is_mixed_external(label):
    """True iff `label` is an external label WITH a trailing "*" -- a
    legacy (pre-2026-08-25) concept. Always False on freshly-produced
    ``tr_label`` values, since Step 3 no longer gates purity."""
    return is_external_label(label) and label.endswith("*")


def is_translation(label):
    """True iff `label` is an external label (clean or, legacy, mixed) whose
    slot is a translation (Tx/Ty/Tz)."""
    axis = external_axis(label)
    return axis is not None and axis[0] == "T"


def is_rotation(label):
    """True iff `label` is an external label (clean or, legacy, mixed) whose
    slot is a rotation (Rx/Ry/Rz)."""
    axis = external_axis(label)
    return axis is not None and axis[0] == "R"


def classification_bucket(label):
    """Map a label string ("Tx".."Rz", or "S"/"B"/"SB") to its semantic
    bucket: "translation"/"rotation" for external slots, "stretch"/"bend"/
    "mixed" for Step-2 vibrational labels. Legacy starred labels ("Tx*")
    still map to "mixed_external" for backward compatibility reading old
    data, but that case is never produced fresh anymore. Falls back to
    returning `label` unchanged."""
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


def predicted_category(tr_label, vib_label_value):
    """The canonical predicted-category rule for one already-scored mode,
    now that Step 3 no longer gates: the Step-3 winning slot (e.g. "Tx") if
    this mode was assigned one, else its own already-computed Step-2
    vib_label ("S"/"B"/"SB"). Both inputs are themselves already-derived
    labels -- this function does no scoring of its own, it only picks
    between two columns/fields a caller supplies (e.g.
    library_scores.csv's ``predicted_tr_label``/``predicted_vib_label``,
    ``classify_all_modes()``'s own ``m["tr_label"]``/``m["vib_label"]``, or
    a freshly re-derived vib_label under a different scheme/thresholds).
    Reused by ``src/calibrate.py``'s ``confusion_matrix_stats`` and
    ``src/figures.py``'s joint-confusion-table helpers so the rule never
    drifts between call sites.
    """
    if isinstance(tr_label, str) and tr_label.strip():
        return tr_label
    return vib_label_value


def predicted_category_column(tr_col, vib_col):
    """Vectorized ``predicted_category()`` over two pandas Series (a
    tr_label-like column and a vib_label-like column of the same length/
    index) -- handles pandas' NaN-for-empty-cell CSV round-trip."""
    has_tr = tr_col.notna() & (tr_col.astype(str).str.strip() != "")
    return tr_col.where(has_tr, vib_col)

# Relative tolerance for the "smallest principal moment ~= 0" linear-molecule test.
LINEAR_TOL = 1e-6


@dataclass
class Thresholds:
    """Step 2/4 thresholds (Step 3, the optional T/R assignment, uses no
    threshold at all -- see module docstring).

    tau_TR : DIAGNOSTIC-ONLY as of 2026-08-25 (no longer gates anything in
        the pipeline) -- a purity bar near 1 a human/diagnostic can compare
        a tr_score against by eye. Kept, not deleted, since
        src/calibrate.py's tau_TR sensitivity sweep still reports it as a
        descriptive statistic of the assignment's own score distribution.
    tau_S  : stretching bar on s[V_S] (>= -> STRETCHING). scheme="threeway" only.
    tau_B  : bending bar on s[V_S] (<= -> BENDING). scheme="threeway" only.
    tau_SB : single-cutoff S/B split used by the *binary* classification
        scheme (scheme="binary"); unrelated to tau_S/tau_B's three-way split
        (a different threshold, not a synonym -- distinct name deliberately
        chosen to avoid collision). v >= tau_SB -> STRETCHING, else BENDING;
        never produces "SB".
    v_weighting : which eq:vscore bond weighting these thresholds were
        calibrated against ('mu' or 'none'), or '*' to match any. tau_S and
        tau_B are read off an s[V_S] distribution, so they are only meaningful
        against the definition that produced it -- classify_all_modes() refuses
        a mismatch rather than silently mislabelling modes.

    Defaults (0.95/0.9/0.2/0.42) are the provisional pre-calibration
    constants, deliberately NOT auto-overwritten by Phase-3 calibration --
    use `Thresholds.calibrated()` for the calibrated values instead. tau_SB
    is a partial exception to that "provisional" framing: unlike tau_TR/
    tau_S/tau_B, no sweep ever computes a replacement for it (Phase-3
    calibration leaves it untouched -- see run_calibration_pipeline in
    src/calibrate.py), so THIS field default is the actual, sole source of
    truth for the canonical tau_SB value; there is no separate "calibrated"
    tau_SB to defer to. Canonical value 0.42 (the all-molecules
    classification-error-minimizing value from the tau_SB_error_sweep -- see
    IMPLEMENTATION_PLAN.md's Recent history).
    tests/test_classifier.py pins `Thresholds()` explicitly so its regression
    goldens stay fixed even if thresholds.json is later recalibrated;
    calibrated behavior has its own tests (tests/test_calibrate.py).
    """
    tau_TR: float = 0.95
    tau_S: float = 0.9
    tau_B: float = 0.2
    tau_SB: float = 0.42
    v_weighting: str = DEFAULT_V_WEIGHTING

    @classmethod
    def calibrated(cls, path=DEFAULT_CALIBRATION_PATH):
        """Load calibrated thresholds from `path` if it exists, else the class defaults."""
        if os.path.exists(path):
            with open(path) as f:
                data = json.load(f)
            # A thresholds.json written before the weighting variant existed
            # carries no stamp, and was by definition calibrated unweighted.
            # A thresholds.json written before tau_SB existed carries no such
            # key either -- fall back to the class default (0.42) rather than
            # KeyError.
            return cls(tau_TR=data["tau_TR"], tau_S=data["tau_S"], tau_B=data["tau_B"],
                       tau_SB=data.get("tau_SB", 0.42),
                       v_weighting=data.get("v_weighting", "none"))
        return cls()

    @classmethod
    def bootstrap(cls):
        """Provisional thresholds that match ANY weighting (sentinel '*').

        Needed for the chicken-and-egg first pass after switching variants:
        building library_scores.csv requires thresholds, but the calibrated
        ones are still stamped for the old definition. Safe because tau_S/tau_B
        derivation reads V_Stretch only and is itself threshold-independent --
        only the vib_label column is affected, which is why --library is run
        again after --calibrate.
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
    """Binary-scheme Step-2 vibrational sub-classification from s[V_S]: forces
    every mode to STRETCHING or BENDING via a single cutoff, never MIXED_STRETCH_BEND."""
    return STRETCHING if v >= tau_SB else BENDING


def vib_label(v, thresholds, scheme="binary"):
    """Step-2 vibrational sub-classification from s[V_S] -- runs
    UNCONDITIONALLY on every mode (see module docstring), independent of
    whether that mode also wins a Step-3 T/R slot.

    scheme="binary" (default, paper-standard as of the 2026-08 binary-
    classification decision): delegates to vib_label_binary(v,
    thresholds.tau_SB) -- forces S or B, never SB. scheme="threeway"
    (explicit opt-in, kept fully functional for comparison): the three-way
    tau_S/tau_B split (may produce MIXED_STRETCH_BEND).
    """
    if scheme == "binary":
        return vib_label_binary(v, thresholds.tau_SB)
    if v >= thresholds.tau_S:
        return STRETCHING
    if v <= thresholds.tau_B:
        return BENDING
    return MIXED_STRETCH_BEND


def rescheme_internal_label(predicted_label, v_stretch, thresholds, scheme):
    """Cheaply re-derive a mode's vibrational label under a DIFFERENT
    `scheme` than the one it was originally classified with, from its own
    `v_stretch` (s[V_S], Step 1's score -- always scheme-independent) alone,
    with NO re-run of Step 1/3. `predicted_label` is the mode's already-
    computed label (under whichever scheme produced it, e.g. a
    library_scores.csv row's combined predicted-category string); if its
    bucket is a vibrational one (stretch/bend/mixed) it is recomputed via
    vib_label(v_stretch, thresholds, scheme); a Step-3 external-slot label
    (e.g. "Tx") is returned unchanged, since scheme never touches Step 3.

    Unchanged since before the 2026-08-25 Step-3-degating restructuring --
    still the general form of the per-row re-labeling pattern used ad hoc in
    a few places (e.g. figures.py's plot_confusion_matrix_binary/
    plot_transferability_confusion_binary/plot_benzene_internal_confusion_binary's
    own inline `_binarize` closures) -- lets a caller obtain ANY scheme's
    vib_label from an already-scored table without re-parsing/re-running the
    whole classifier.
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


def classify_all_modes(scorer, final, thresholds=None, scheme="binary", identify_tr=True):
    """Algorithm 1: score every mode (Step 1), classify every mode's
    vibrational character unconditionally (Step 2), then OPTIONALLY assign
    the best n_T+n_R modes to T/R slots with no purity gate (Step 3).

    scorer must already be MIT-aligned (as build_scorer_and_final leaves it);
    final is its candidate mode pool. thresholds defaults to
    Thresholds.calibrated(). `scheme` is "binary" (default, paper-standard:
    single tau_SB cutoff, never produces "SB") or "threeway" (explicit
    opt-in, may produce "SB"; kept fully functional for comparison/on-demand
    use) -- affects Step 2's vib_label vocabulary ONLY; it does not affect
    Step 3's assignment (which uses no threshold or scheme at all).
    `identify_tr` (default True) toggles Step 3 -- when False, every mode's
    tr_label/tr_score stay None and only vib_label is populated; every
    reproduce.py/headless default keeps this True.

    Returns a list of dicts, one per mode in `final` (same order): {name,
    frequency, is_emit, T, R, V, vib_label, tr_label, tr_score, bonds,
    bonds_all}. 'bonds' (see ModeScorer.score_bonds()) carries the per-bond
    s_AB list for every mode regardless of label; 'bonds_all' is the same
    list under the name ped/merge_ped_scores.py's per-bond-type breakdown
    reads. 'vib_label' is "S"/"B" (or "S"/"B"/"SB" under scheme="threeway"),
    always populated. 'tr_label' is one of "Tx".."Rz" if this mode won that
    slot in Step 3, else None (including whenever identify_tr=False).
    'tr_score' is that slot's own signed score value, or None to match.
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
            "vib_label": None,
            "tr_label": None,
            "tr_score": None,
        })

    n_modes = len(scored)

    # ---- Step 2: vibrational label, UNCONDITIONAL, every mode ----
    for m in scored:
        m["vib_label"] = vib_label(m["V"], thresholds, scheme)

    # ---- Step 3: OPTIONAL, ungated T/R identification ----
    if identify_tr:
        cost = np.zeros((len(slots), n_modes))
        for si, slot in enumerate(slots):
            for mi, m in enumerate(scored):
                cost[si, mi] = abs(_score_slot(m, slot))

        try:
            row_ind, col_ind = linear_sum_assignment(cost, maximize=True)
        except TypeError:
            # Fallback for older scipy without the maximize kwarg: negate and minimize.
            row_ind, col_ind = linear_sum_assignment(-cost)

        for si, mi in zip(row_ind, col_ind):
            slot = slots[si]
            scored[mi]["tr_label"] = slot
            scored[mi]["tr_score"] = _score_slot(scored[mi], slot)

    # 'bonds' carries the per-bond s_AB list for every mode regardless of
    # classification; 'bonds_all' is kept as an alias (same list) for
    # ped/merge_ped_scores.py, which reads it by that name.
    for m in scored:
        m["bonds_all"] = m.pop("_bonds_all")
        m["bonds"] = m["bonds_all"]

    return scored


def classify_to_rows(scored):
    """Flatten classify_all_modes() output into CSV-row dicts -- this is the
    single, complete per-mode result row (scores + Mu/K/Irrep + Step 2/3
    labels); there is no separate scores-only row shape. Each bond gets its
    own 's_AB[Elem#-Elem#]' column (e.g. "s_AB[C1-C2]") rather than one
    semicolon-joined string, so the CSV is a plain rectangular table --
    every row lists the same molecule's bonds in the same order (see
    ModeScorer.bList), so the column set is identical across rows and easy
    to sort/filter/plot in Excel. Mu/K/Irrep are None for EMIT modes and the
    synthetic ideal T/R references (only real Gaussian normal modes carry
    them). 'tr_label' is written as "" (not NaN) when no slot was assigned,
    matching the pre-existing empty-string-for-absent convention; 'tr_score'
    is left as Python None, which pandas renders as an empty cell."""
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
            "vib_label": m["vib_label"],
            "tr_label": m["tr_label"] or "",
            "tr_score": m["tr_score"],
        }
        for b in m["bonds"]:
            row[f"s_AB[{b['i_label']}-{b['j_label']}]"] = b["s_AB"]
        rows.append(row)
    return rows
