"""pipeline.py -- parse -> align -> score -> classify, in one call.

Mirrors ``main.build_scorer_and_final`` + ``classify_all_modes`` from the
reference implementation so this app and the CLI agree by construction. The
whole thing runs in ~11 ms for naphthalene (18 atoms, 54 modes), which is why
there is no job queue, no polling and no server-side run state anywhere in
this app.

Everything the browser needs -- table rows, the 3Dmol payload, provenance --
comes out of :func:`analyse` as one JSON-safe dict.
"""

from __future__ import annotations

import numpy as np

from .classifier import Thresholds, classify_all_modes, is_linear
from .parsers import ParseError, displacement_precision
from .scoring import ModeScorer

SCORE_KEYS = ("Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "V_S")

# Labels produced by Step 4 of the classifier, with the plain-language reading
# the UI shows. External slots (Tx..Rz) are described in place.
LABEL_TEXT = {
    "S": "stretching",
    "B": "bending",
    "SB": "mixed stretch/bend",
}


def _rotation_between(a, b):
    """Angle in degrees of the rigid rotation taking geometry `a` onto `b`.

    Both are COM-centred, so the optimal rotation is the Kabsch one; the angle
    comes from its trace. Used only to report how far the principal-axis
    alignment turned the molecule.
    """
    if len(a) < 2:
        return 0.0
    u, _, vt = np.linalg.svd(a.T @ b)
    d = np.sign(np.linalg.det(u @ vt))
    r = u @ np.diag([1.0, 1.0, d]) @ vt
    cos = (np.trace(r) - 1.0) / 2.0
    return round(float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))), 1)


def _detect_mode_set(n_modes, n_atoms, n_ext, linear):
    """Work out from the counts alone whether a file holds 3N or 3N-6 modes.

    The two are never ambiguous -- 3N and 3N-n_ext differ by 5 or 6 -- so the
    user does not need to declare it. Anything else is refused rather than
    guessed at: a count that matches neither means the file is incomplete or
    holds something other than one molecule's modes, and either way the six
    external slots would be filled from the wrong pool.
    """
    full = 3 * n_atoms
    vib_only = full - n_ext
    if n_modes == full:
        return "3n"
    if n_modes == vib_only:
        return "3n-6"
    shape = "3N-5" if linear else "3N-6"
    raise ParseError(
        f"The file has {n_modes} modes, which is neither 3N = {full} (every mode, "
        f"translations and rotations included) nor {shape} = {vib_only} "
        f"(vibrations only) for {n_atoms} atoms. Modes cannot be classified from "
        "a partial set: the six external slots would be filled from whatever was "
        "supplied.")


def analyse(atoms, coords, bonds, modes, title="", source="", warnings=None,
            mode_set="auto"):
    """Score and classify one molecule. Returns a JSON-safe payload.

    ``coords``/``modes`` are taken in whatever frame they arrived in; the
    scorer rotates BOTH into the principal-axis frame together (MIT), which is
    the frame the T/R scores are defined in and the frame the viewer shows.

    ``mode_set`` says what the supplied modes are, which decides whether the
    ideal T/R references are added to the pool:

    "3n-6"  vibrations only (3N-6, or 3N-5 linear). The six external slots have
            no genuine candidates among the supplied modes, so the constructed
            ideal T/R modes are prepended to fill them. Without that, Step 2's
            assignment claims real vibrations instead -- on H2O it labelled a
            0.997 O-H stretch "Tx*".

    "3n"    the complete 3N set, T/R included. The external slots are filled
            from the supplied modes themselves and nothing is constructed, so
            the six T/R modes are identified by the algorithm rather than
            supplied alongside it.
    """
    warnings = list(warnings or [])
    notes = []                      # informational; not a problem with the input
    if mode_set not in ("auto", "3n", "3n-6"):
        raise ValueError(
            f"mode_set must be 'auto', '3n' or '3n-6', got {mode_set!r}")

    if not bonds:
        # Belt and braces: every reader already refuses this, because an empty
        # bond list makes Vscore() return exactly 0.0 for every mode, and
        # 0.0 <= tau_B, so the whole table would come back labelled "B".
        raise ParseError(
            "No bond connectivity available. The V-score is a sum over bonds; "
            "without them every mode would score 0.000 and be labelled bending.")

    thresholds = Thresholds.calibrated()
    scorer = ModeScorer(atoms, coords, bonds)
    _as_given = np.array(scorer.coords, dtype=float)   # COM-shifted, not yet aligned

    rotated = scorer.MIT([dict(m) for m in modes], rotate_modes=True)
    frame_rotation = _rotation_between(_as_given, np.array(scorer.coords, dtype=float))
    for i, m in enumerate(rotated):
        m["label"] = f"Vib {i + 1}"

    linear = is_linear(scorer)
    ideal_R = scorer.construct_R()
    if linear:
        # MIT always places the linear axis on x, so construct_R()'s "Rx" is an
        # all-zero placeholder rather than a real external mode. Left in, it
        # falls through to Step 4 and is mislabelled BENDING (V=0 <= tau_B).
        ideal_R = [m for m in ideal_R if m["label"] != "Rx"]

    n_ext = 3 + len(ideal_R)                 # external slots: 3 T + 3 R (2 linear)
    detected = _detect_mode_set(len(modes), len(atoms), n_ext, linear)
    if mode_set == "auto":
        mode_set = detected
    elif mode_set != detected:
        raise ParseError(
            f"This file holds {len(modes)} modes, i.e. {detected}, "
            f"but {mode_set} was requested.")

    if mode_set == "3n":
        final = rotated                      # nothing constructed; T/R must be in here
        n_ref = 0
    else:
        final = scorer.construct_T() + ideal_R + rotated
        n_ref = len(scorer.construct_T()) + len(ideal_R)

    scored = classify_all_modes(scorer, final, thresholds)
    # classify_all_modes returns one entry per mode in `final`, in the same
    # order, but drops the displacement vector -- the viewer needs it, so pair
    # them back up here rather than re-deriving.
    rows = [_row(m, vec["vector"], i, is_reference=(i < n_ref))
            for i, (m, vec) in enumerate(zip(scored, final))]
    references, vibrations = rows[:n_ref], rows[n_ref:]

    if mode_set == "3n":
        starred = [r for r in vibrations if r["label"].endswith("*")]
        assigned = [r for r in vibrations if r["label"][0] in "TR"]
        if starred:
            # Not a fault: a starred label is the algorithm reporting that the
            # mode is genuinely part rigid-body motion and part vibration.
            # Routine for EMIT, where T/R character is spread across the basis.
            notes.append(
                f"{len(starred)} mixed external"
                + ("s" if len(starred) > 1 else "") + ": "
                + ", ".join(f"{r['name']} {r['label']}" for r in starred)
                + ". Each is part rigid-body motion, part vibration — both its "
                  "V_S and its largest T/R score are highlighted.")
        if len(assigned) < n_ext:
            warnings.append(
                f"Only {len(assigned)} of the {n_ext} external slots were "
                "assigned, so the file may not hold a complete 3N mode set.")

    dp = displacement_precision(modes)
    if dp <= 2:
        warnings.append(
            f"Displacements carry only {dp} decimal places, which is standard "
            "Gaussian output rather than freq=hpmodes. The V-score and the "
            "S/B/SB labels are robust to this, but individual T-scores can shift "
            "by up to ~0.3 because Tscore() unit-normalises each atom's "
            "displacement. Re-run with freq=hpmodes for publication numbers.")

    return {
        "title": title or "molecule",
        "source": source,
        "n_atoms": len(atoms),
        "n_modes": len(vibrations),
        "linear": bool(linear),
        "atoms": list(atoms),
        "bonds": [[int(a), int(b)] for a, b in bonds],
        "bond_labels": [f"{atoms[a]}{a+1}-{atoms[b]}{b+1}" for a, b in bonds],
        "geometry": np.asarray(scorer.coords).round(6).tolist(),
        # How far MIT turned the molecule from the frame it arrived in. The
        # scores -- and the axes the viewer draws -- are defined in the rotated
        # frame, so a large angle here is why the picture looks turned relative
        # to the input deck. It is not a disagreement between the two.
        "frame_rotation_deg": frame_rotation,
        "vibrations": vibrations,
        "references": references,
        "thresholds": {
            "tau_TR": thresholds.tau_TR,
            "tau_S": thresholds.tau_S,
            "tau_B": thresholds.tau_B,
            "v_weighting": thresholds.v_weighting,
        },
        "precision_dp": dp,
        "mode_set": mode_set,
        "warnings": warnings,
        "notes": notes,
        "freq_range": _freq_range(vibrations),
    }


def _row(m, vector, index, is_reference=False):
    """One table row: the seven scores, the label, and which score dominates."""
    scores = {
        "Tx": m["T"]["x"], "Ty": m["T"]["y"], "Tz": m["T"]["z"],
        "Rx": m["R"]["x"], "Ry": m["R"]["y"], "Rz": m["R"]["z"],
        "V_S": m["V"],
    }
    # Largest |value|, not largest value: s[V_S] alone cannot be negative, so a
    # raw comparison would structurally favour it and would pick +0.087 over
    # -0.283. A mode translating along -x is as translational as one along +x.
    dominant = max(SCORE_KEYS, key=lambda k: abs(scores[k]))

    label = m["classification"]

    # What the UI highlights follows the LABEL, not the raw abs-max, because
    # the label is what the algorithm actually concluded:
    #
    #   clean external (Tx..Rz)  -> that slot's score. The mode is a rigid-body
    #                               motion; its s[V_S] is ~0 by the purity gate.
    #   internal (S/B/SB)        -> s[V_S], the score the label is derived from.
    #                               The T/R columns are diagnostics here, and
    #                               highlighting an abs-max that lands there
    #                               (43% of real modes) would assert external
    #                               character the algorithm never assigned.
    #   mixed external (Tx*..)   -> BOTH: s[V_S] and the largest |T/R|. The star
    #                               means the mode failed the purity gate, i.e.
    #                               it is genuinely part external and part
    #                               vibration, so one number cannot describe it.
    _TR = tuple(k for k in SCORE_KEYS if k != "V_S")
    is_external = label[0] in "TR"
    if is_external and label.endswith("*"):
        highlights = ["V_S", max(_TR, key=lambda k: abs(scores[k]))]
    elif is_external:
        highlights = [label]
    else:
        highlights = ["V_S"]
    highlight = highlights[0]
    return {
        "index": index,
        "name": m["name"],
        "is_reference": bool(is_reference),
        "highlight": highlight,
        "highlights": highlights,
        "frequency": _num(m["frequency"]),
        "scores": {k: round(float(v), 4) for k, v in scores.items()},
        "dominant": dominant,
        "dominant_value": round(float(scores[dominant]), 4),
        "label": label,
        "label_text": _label_text(label),
        "annotation": m["annotation"],
        "irrep": m["irrep"],
        "reduced_mass": _num(m["reduced_mass"]),
        "force_constant": _num(m["force_constant"]),
        "bonds": [{"pair": f"{b['i_label']}-{b['j_label']}",
                   "s_AB": round(float(b["s_AB"]), 4),
                   "rel_db": round(float(b["rel_db"]), 5)} for b in m["bonds"]],
        "vector": np.asarray(vector).round(5).tolist(),
    }


def _label_text(label):
    if label in LABEL_TEXT:
        return LABEL_TEXT[label]
    if label.endswith("*"):
        return f"mixed external ({label[:-1]}) + vibration"
    return f"external: {label}"


def _num(v):
    if v is None:
        return None
    f = float(v)
    return f if np.isfinite(f) else None


def _freq_range(rows):
    freqs = [r["frequency"] for r in rows if r["frequency"] is not None]
    if not freqs:
        return [0.0, 0.0]
    return [float(min(freqs)), float(max(freqs))]


def to_csv_rows(payload, include_references=False):
    """Flatten to CSV rows in the same shape as the reference implementation's
    ``classify_to_rows`` output, so a download drops straight into the existing
    analysis (one ``s_AB[X#-Y#]`` column per bond)."""
    rows = []
    source = (payload["references"] + payload["vibrations"]
              if include_references else payload["vibrations"])
    for r in source:
        row = {
            "Mode": r["name"],
            "Freq": r["frequency"],
            "Tx": r["scores"]["Tx"], "Ty": r["scores"]["Ty"], "Tz": r["scores"]["Tz"],
            "Rx": r["scores"]["Rx"], "Ry": r["scores"]["Ry"], "Rz": r["scores"]["Rz"],
            "V_Stretch": r["scores"]["V_S"],
            "Mu": r["reduced_mass"], "K": r["force_constant"], "Irrep": r["irrep"],
            "label": r["label"],
            "annotation": r["annotation"],
            # highlighted_score is what the table marks: the dominant score on a
            # constructed T/R reference, and always V_S on a real vibration.
            # dominant_score is the raw abs-max across all seven, kept because it
            # is a real diagnostic -- but it is NOT a claim of external character
            # (no real vibration is ever assigned a T/R slot).
            "highlighted_score": "+".join(r["highlights"]),
            "dominant_score": r["dominant"],
            "is_reference_mode": r["is_reference"],
        }
        for b in r["bonds"]:
            row[f"s_AB[{b['pair']}]"] = b["s_AB"]
        rows.append(row)
    return rows
