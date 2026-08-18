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


_EXTERNAL_SLOTS = {"Tx", "Ty", "Tz", "Rx", "Ry", "Rz"}


def analyse(atoms, coords, bonds, modes, title="", source="", warnings=None):
    """Score and classify one molecule. Returns a JSON-safe payload.

    ``coords``/``modes`` are taken in whatever frame they arrived in; the
    scorer rotates BOTH into the principal-axis frame together (MIT), which is
    the frame the T/R scores are defined in and the frame the viewer shows.

    Mirrors ``main.build_scorer_and_final``'s two branches:
      - normal modes: MIT rotates the molecule AND the mode vectors
        (rotate_modes=True); 3 ideal T + 3 ideal R references are prepended,
        and (by construction, verified 0/463 counterexamples in the reference
        set) always win their own slot in Step 2's assignment -- so a mode's
        POSITION in `final` reliably tells reference from real vibration.
      - EMIT modes: already live in principal axes, so MIT rotates only the
        molecule (rotate_modes=False); no ideal references are added, all 3N
        raw EMIT eigenvectors are the candidate pool, and REAL eigenvectors
        do win external slots (that's the whole point of the EMIT test) -- so
        here the CLASSIFICATION, not position, says which modes are external.
    """
    warnings = list(warnings or [])

    if not bonds:
        # Belt and braces: every reader already refuses this, because an empty
        # bond list makes Vscore() return exactly 0.0 for every mode, and
        # 0.0 <= tau_B, so the whole table would come back labelled "B".
        raise ParseError(
            "No bond connectivity available. The V-score is a sum over bonds; "
            "without them every mode would score 0.000 and be labelled bending.")

    thresholds = Thresholds.calibrated()
    scorer = ModeScorer(atoms, coords, bonds)

    is_emit = any(m.get("is_emit") for m in modes)
    rotated = scorer.MIT([dict(m) for m in modes], rotate_modes=not is_emit)

    linear = is_linear(scorer)

    if is_emit:
        final = rotated
        scored = classify_all_modes(scorer, final, thresholds)
        rows = [_row(m, vec["vector"], i,
                     is_reference=(m["classification"] in _EXTERNAL_SLOTS))
                for i, (m, vec) in enumerate(zip(scored, final))]
        references, vibrations = [], rows
    else:
        for i, m in enumerate(rotated):
            m["label"] = f"Vib {i + 1}"

        ideal_R = scorer.construct_R()
        if linear:
            # MIT always places the linear axis on x, so construct_R()'s "Rx"
            # is an all-zero placeholder rather than a real external mode.
            # Left in, it falls through to Step 4 and is mislabelled BENDING
            # (V=0 <= tau_B).
            ideal_R = [m for m in ideal_R if m["label"] != "Rx"]

        final = scorer.construct_T() + ideal_R + rotated
        scored = classify_all_modes(scorer, final, thresholds)

        n_ref = len(scorer.construct_T()) + len(ideal_R)
        # classify_all_modes returns one entry per mode in `final`, in the
        # same order, but drops the displacement vector -- the viewer needs
        # it, so pair them back up here rather than re-deriving.
        rows = [_row(m, vec["vector"], i, is_reference=(i < n_ref))
                for i, (m, vec) in enumerate(zip(scored, final))]
        references, vibrations = rows[:n_ref], rows[n_ref:]

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
        "vibrations": vibrations,
        "references": references,
        "thresholds": {
            "tau_TR": thresholds.tau_TR,
            "tau_S": thresholds.tau_S,
            "tau_B": thresholds.tau_B,
            "v_weighting": thresholds.v_weighting,
        },
        "precision_dp": dp,
        "warnings": warnings,
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

    # What the UI highlights, which is NOT always the abs-max. Only the six
    # constructed T/R reference modes are candidates for external character --
    # Step 2 assigns the external slots to them every time, and no real
    # vibration has ever taken one (0 of 463 across the reference set). So on a
    # real mode the T/R columns are diagnostics, not a claim about its
    # character, and highlighting an abs-max that lands there (43% of modes)
    # would assert translation the algorithm never assigned. Real modes
    # therefore highlight s[V_S], the score their label actually comes from.
    # (`is_reference` means "position-based ideal T/R row" for normal modes
    # but "classification-based clean external" for EMIT -- see analyse().)
    highlight = dominant if is_reference else "V_S"

    label = m["classification"]
    return {
        "index": index,
        "name": m["name"],
        "is_reference": bool(is_reference),
        "is_emit": bool(m["is_emit"]),
        "highlight": highlight,
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
            ("Eigenvalue" if r["is_emit"] else "Freq"): r["frequency"],
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
            "highlighted_score": r["highlight"],
            "dominant_score": r["dominant"],
            "is_reference_mode": r["is_reference"],
        }
        for b in r["bonds"]:
            row[f"s_AB[{b['pair']}]"] = b["s_AB"]
        rows.append(row)
    return rows
