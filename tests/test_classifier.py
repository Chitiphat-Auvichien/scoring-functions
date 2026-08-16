"""Golden-reference regression tests for src/classifier.py (Algorithm 1).

Run from ``Github/scoring-functions/``:
    py -m pytest tests/                (with pytest)
    py tests/test_classifier.py       (standalone; no pytest needed)

Pins the label-level validation targets from IMPLEMENTATION_PLAN.md: water's
6 constructed T/R references classify clean, water's 3 real vibrations split
BENDING/STRETCHING/STRETCHING, and benzene EMIT 34/35 (stretching-mixed
externals) vs EMIT 36 (the documented out-of-plane-bending blind spot, which
is INTENTIONALLY reported CLEAN_TRANSLATION -- see IMPLEMENTATION_PLAN.md
Decision X) resolve as designed.

Threshold policy (decided Phase 3, 2026-07): this file explicitly pins
Thresholds() -- the hardcoded PROVISIONAL constants (tau_TR=0.95, tau_S=0.9,
tau_B=0.2) -- rather than relying on classify_all_modes()'s default (which
now auto-loads the Phase-3 CALIBRATED values via Thresholds.calibrated() when
data/results/thresholds.json exists). This keeps these regression goldens
fixed and reproducible even if a future recalibration (e.g. an expanded
library) changes thresholds.json's numbers. The calibrated-threshold behavior
is independently re-verified against the SAME targets in
tests/test_calibrate.py, which loads Thresholds.calibrated() explicitly --
see that file for the confirmation that calibration did not change any of
these labels (tau_S 0.9->0.9037, tau_B 0.2->0.1733, tau_TR unchanged at
0.95; every target below is robust to that shift).

**2026-07-07 basename update:** "water"/"benzene" -> the finalized roster
basenames `H2O`/`C6H6` (data/mol_list_method.csv).
Every assertion in this file is CLASSIFICATION-LEVEL (bucket labels like
BENDING/STRETCHING/"Tx", bond counts), not a pinned V_Stretch/T/R number, so
none needed updating -- water's new file's V_Stretch values differ from the
old file's (see tests/test_scores.py's module docstring) but land on the
same side of tau_S/tau_B either way.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from main import load_inputs, build_scorer_and_final              # noqa: E402
from src.classifier import (                                       # noqa: E402
    classify_all_modes, classify_to_rows, Thresholds,
    is_mixed_external,
    STRETCHING, BENDING, MIXED_STRETCH_BEND,
)

TOL = 5e-4  # 3 decimal places

# Explicit provisional thresholds -- see module docstring for why this file
# does not rely on classify_all_modes()'s calibrated-by-default resolution.
_PROVISIONAL = Thresholds()


def _classify(mol, mode_type):
    raw, _ = load_inputs(mol, mode_type, os.path.join(ROOT, "data"))
    scorer, final = build_scorer_and_final(raw, mode_type)
    scored = classify_all_modes(scorer, final, _PROVISIONAL)
    return {m["name"]: m for m in scored}


def test_water_normal_externals_clean():
    """The 6 ideal T/R reference modes classify clean (score>=tau_TR, V<=tau_B) --
    i.e. the bare slot name itself ("Tx".."Rz"), per the 2026-07-02 axis-specific
    label rename (no more generic CLEAN_TRANSLATION/CLEAN_ROTATION constants)."""
    t = _classify("H2O", "normal")
    for lbl in ("Tx", "Ty", "Tz", "Rx", "Ry", "Rz"):
        assert t[lbl]["classification"] == lbl, t[lbl]["classification"]


def test_water_normal_vibrations():
    """Vib1 = bend -> BENDING; Vib2/Vib3 = stretches -> STRETCHING, with |s_AB| summing to V
    (s_AB is signed; V sums magnitudes)."""
    t = _classify("H2O", "normal")
    assert t["Vib 1"]["classification"] == BENDING
    assert len(t["Vib 1"]["bonds"]) == 2  # bonds are attached for every mode, not just S/SB

    for name in ("Vib 2", "Vib 3"):
        assert t[name]["classification"] == STRETCHING
        bonds = t[name]["bonds"]
        assert len(bonds) == 2  # water has 2 O-H bonds
        assert abs(sum(abs(b["s_AB"]) for b in bonds) - t[name]["V"]) < 1e-6


def test_benzene_emit_34_35_mixed_external():
    """EMIT 34/35 (E1u, |s[T]|=1 but V_Stretch=0.667/0.577 > tau_B) fail gate 2 -> flagged
    with the axis-specific mixed-external label ("Tx*"/"Ty*", 2026-07-02 rename); the axis
    is now IN the classification itself, so the annotation only carries the vibration
    sub-label (no more redundant "dominant_external=..." text)."""
    t = _classify("C6H6", "emit")
    assert t["EMIT 34"]["classification"] == "Tx*"
    assert is_mixed_external(t["EMIT 34"]["classification"])
    assert t["EMIT 34"]["annotation"] == f"vibration={MIXED_STRETCH_BEND}"
    assert t["EMIT 35"]["classification"] == "Ty*"
    assert t["EMIT 35"]["annotation"] == f"vibration={MIXED_STRETCH_BEND}"


def test_benzene_emit_36_clean_translation_blind_spot():
    """EMIT 36 (A2u, |s[Tz]|=1, V_Stretch=0) passes BOTH gates -> bare "Tz" (clean).

    This is the documented, INTENTIONAL blind spot (IMPLEMENTATION_PLAN.md Decision X):
    the two-gate purity test cannot see EMIT 36's out-of-plane bending residual because
    s[V_S]=0 for it too (bending, not stretching). Reproducing clean "Tz" here is
    correct -- do not "fix" this.
    """
    t = _classify("C6H6", "emit")
    assert t["EMIT 36"]["classification"] == "Tz"
    assert abs(t["EMIT 36"]["V"]) < TOL


def test_co2_linear_no_spurious_onaxis_mode():
    """Linear-molecule candidate pool must have exactly 3N modes, not 3N+1.

    Regression for a bug found by formula-auditor: construct_R() always builds
    3 ideal rotation references, but for a linear molecule MIT() places the
    molecular (smallest-moment) axis on the new X axis, so the ideal "Rx"
    reference is an all-zero vector (n_R=2 excludes it, spec). Left in the
    pool, it fell through Step 2 (no slot claims an all-zero mode) into Step 4
    and was mislabeled BENDING (V=0 <= tau_B). build_scorer_and_final() must
    drop it before scoring, and the on-axis slot must simply be absent (not
    present-and-clean, not present-and-mislabeled).
    """
    t = _classify("CO2", "normal")
    assert len(t) == 9  # 3N for a 3-atom linear molecule, not 3N+1
    assert "Rx" not in t
    assert t["Ry"]["classification"] == "Ry"
    assert t["Rz"]["classification"] == "Rz"
    # the 2 bends and 2 stretches must be the REAL vibrational modes, not a
    # spurious all-zero placeholder
    vib_labels = {name: m["classification"] for name, m in t.items() if name.startswith("Vib")}
    assert sorted(vib_labels.values()) == [BENDING, BENDING, STRETCHING, STRETCHING]


def test_classify_to_rows_shape():
    """classify_to_rows() produces the documented CSV columns for both mol/mode_type combos,
    with one 's_AB[iLabel-jLabel]' column per bond (not a joined string) -- and every row
    (every mode) carries the same bond columns, since a molecule's bond list doesn't
    depend on which mode is being scored."""
    base_cols = {"Mode", "Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "V_Stretch",
                 "Mu", "K", "Irrep", "label", "annotation"}
    for mol, mode_type, freq_col in [("H2O", "normal", "Freq"),
                                      ("C6H6", "emit", "Eigenvalue")]:
        raw, _ = load_inputs(mol, mode_type, os.path.join(ROOT, "data"))
        scorer, final = build_scorer_and_final(raw, mode_type)
        scored = classify_all_modes(scorer, final)
        rows = classify_to_rows(scored)
        assert rows, "no rows produced"
        bond_cols = {f"s_AB[{b['i_label']}-{b['j_label']}]" for b in scored[0]["bonds"]}
        assert bond_cols, "no bond columns produced"
        expected_cols = base_cols | {freq_col} | bond_cols
        for row in rows:
            assert set(row.keys()) == expected_cols, set(row.keys())


def test_rejects_thresholds_calibrated_for_the_other_weighting():
    """tau_S/tau_B are cut points on an s[V_S] distribution, so they only mean
    anything against the definition that produced it. Pairing them with the
    other definition would still label every mode -- just wrongly, and with no
    error to notice. classify_all_modes() refuses instead.
    """
    raw, _ = load_inputs("H2O", "normal", os.path.join(ROOT, "data"))
    scorer, final = build_scorer_and_final(raw, "normal", v_weighting="mu")

    raised = False
    try:
        classify_all_modes(scorer, final, Thresholds(v_weighting="none"))
    except ValueError as e:
        raised = "mismatch" in str(e)
    assert raised, "mu scorer + 'none' thresholds should raise"

    # Matching, and the '*' bootstrap sentinel, both go through.
    assert classify_all_modes(scorer, final, Thresholds(v_weighting="mu"))
    assert classify_all_modes(scorer, final, Thresholds.bootstrap())


if __name__ == "__main__":
    # Standalone runner so the suite works even without pytest installed.
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    passed = 0
    for fn in tests:
        try:
            fn()
            print(f"PASS  {fn.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"FAIL  {fn.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            print(f"ERROR {fn.__name__}: {type(e).__name__}: {e}")
    print(f"\n{passed}/{len(tests)} passed")
    sys.exit(0 if passed == len(tests) else 1)
