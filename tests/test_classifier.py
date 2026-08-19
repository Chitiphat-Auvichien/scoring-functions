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
    is_mixed_external, vib_label, vib_label_binary,
    STRETCHING, BENDING, MIXED_STRETCH_BEND,
    gate2_bar, rescheme_external_label,
)

TOL = 5e-4  # 3 decimal places

# Explicit provisional thresholds -- see module docstring for why this file
# does not rely on classify_all_modes()'s calibrated-by-default resolution.
_PROVISIONAL = Thresholds()


def _classify(mol, mode_type, scheme="threeway"):
    """`scheme` pinned to "threeway" by default: every golden target below
    (water B/S, benzene EMIT 34/35/36 SB annotation, CO2) was derived/pinned
    under the three-way scheme, so this keeps them fixed and reproducible
    even though classify_all_modes()'s own default flipped to "binary" (the
    2026-08 paper-standard switch) -- see
    test_classify_all_modes_default_scheme_is_now_binary below for the
    dedicated check that the new global default is genuinely binary."""
    raw, _ = load_inputs(mol, mode_type, os.path.join(ROOT, "data"))
    scorer, final = build_scorer_and_final(raw, mode_type)
    scored = classify_all_modes(scorer, final, _PROVISIONAL, scheme=scheme)
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


def test_vib_label_binary_boundary_and_no_third_outcome():
    """vib_label_binary forces every value to S or B -- there is no third
    outcome, unlike the threeway split's MIXED_STRETCH_BEND band."""
    tau_SB = 0.5
    assert vib_label_binary(0.0, tau_SB) == BENDING
    assert vib_label_binary(0.499999, tau_SB) == BENDING
    assert vib_label_binary(0.5, tau_SB) == STRETCHING  # boundary: >= -> STRETCHING
    assert vib_label_binary(1.0, tau_SB) == STRETCHING
    # A value that would be MIXED_STRETCH_BEND under the threeway split
    # (strictly between tau_B and tau_S) still lands cleanly on one side here.
    assert vib_label_binary(0.6, tau_SB) in (STRETCHING, BENDING)


def test_vib_label_scheme_param_delegates_to_binary():
    """vib_label(..., scheme="binary") must delegate to vib_label_binary
    (same tau_SB cutoff, same S/B-only vocabulary) rather than reimplementing it."""
    thresholds = Thresholds(tau_S=0.9, tau_B=0.2, tau_SB=0.5)
    for v in (0.0, 0.2, 0.45, 0.5, 0.6, 0.9, 1.0):
        assert vib_label(v, thresholds, scheme="binary") == vib_label_binary(v, thresholds.tau_SB)
    # Threeway is unaffected by tau_SB entirely -- a value in the SB band
    # under threeway must still come out MIXED_STRETCH_BEND there, even
    # though the binary scheme (same thresholds object) would call it S or B.
    assert vib_label(0.45, thresholds, scheme="threeway") == MIXED_STRETCH_BEND
    assert vib_label(0.45, thresholds, scheme="binary") == BENDING


def test_classify_all_modes_binary_scheme_never_produces_sb():
    """On a real mode pool (benzene's normal modes) with known V_Stretch
    values straddling tau_SB=0.5 and both strictly inside the threeway
    SB band (tau_B=0.2 < V < tau_S=0.9), scheme='threeway' must produce
    MIXED_STRETCH_BEND for both, while scheme='binary' must produce S or B
    for EVERY internal mode in the pool (never 'SB') -- and specifically
    diverge on these two modes (Vib 13 -> B, Vib 14 -> S), confirming the
    two schemes are genuinely different, not just synonyms.
    """
    raw, _ = load_inputs("C6H6", "normal", os.path.join(ROOT, "data"))
    scorer, final = build_scorer_and_final(raw, "normal")
    thresholds = Thresholds()  # provisional tau_S=0.9, tau_B=0.2, tau_SB=0.5

    scored_3 = {m["name"]: m for m in classify_all_modes(scorer, final, thresholds, scheme="threeway")}
    scored_b = {m["name"]: m for m in classify_all_modes(scorer, final, thresholds, scheme="binary")}

    # Both real modes' V_Stretch sit strictly inside the threeway SB band.
    assert 0.2 < scored_3["Vib 13"]["V"] < 0.9
    assert 0.2 < scored_3["Vib 14"]["V"] < 0.9
    assert scored_3["Vib 13"]["classification"] == MIXED_STRETCH_BEND
    assert scored_3["Vib 14"]["classification"] == MIXED_STRETCH_BEND

    # Binary scheme diverges: same V_Stretch, different verdict either side
    # of tau_SB=0.5 (Vib 13's V~0.451 < 0.5 -> BENDING; Vib 14's V~0.509 >= 0.5 -> STRETCHING).
    assert scored_b["Vib 13"]["classification"] == BENDING
    assert scored_b["Vib 14"]["classification"] == STRETCHING

    # No internal (non-external-slot) mode in the WHOLE pool is ever labeled
    # "SB" under the binary scheme.
    for m in scored_b.values():
        assert m["classification"] != MIXED_STRETCH_BEND

    # Step 2 (assignment) and gate 1 (tau_TR) are untouched by scheme; gate 2
    # is scheme-dependent as of the tau_purity split (see gate2_bar), so this
    # equality only holds here because benzene's normal-mode external slots
    # are always the trivially-perfect synthetic ideal T/R references
    # (V_Stretch=0 exactly -- see build_scorer_and_final), which pass gate 2
    # under EITHER bar. See test_gate2_scheme_divergence_synthetic_case below
    # for a case constructed to actually diverge on gate 2 itself.
    for name in scored_3:
        if scored_3[name]["classification"] != MIXED_STRETCH_BEND:
            assert scored_3[name]["classification"] == scored_b[name]["classification"], name


def test_classify_all_modes_default_scheme_is_now_binary():
    """2026-08 paper-standard switch: classify_all_modes()'s own default
    (no `scheme` kwarg at all, unlike `_classify()` above which pins
    "threeway" explicitly) must now be BINARY, not threeway. Benzene EMIT
    34/35 (V_Stretch 0.667/0.577, both >= tau_SB=0.5) are STRETCHING under
    binary, never MIXED_STRETCH_BEND -- the mirror-image check of
    test_benzene_emit_34_35_mixed_external's threeway-pinned "SB" result."""
    raw, _ = load_inputs("C6H6", "emit", os.path.join(ROOT, "data"))
    scorer, final = build_scorer_and_final(raw, "emit")
    scored = {m["name"]: m for m in classify_all_modes(scorer, final, _PROVISIONAL)}  # no scheme=
    assert scored["EMIT 34"]["classification"] == "Tx*"
    assert scored["EMIT 34"]["annotation"] == f"vibration={STRETCHING}"
    assert scored["EMIT 35"]["classification"] == "Ty*"
    assert scored["EMIT 35"]["annotation"] == f"vibration={STRETCHING}"
    # Confirms this genuinely differs from the threeway-pinned result above,
    # not just a coincidentally-identical annotation string.
    assert vib_label(0.667, _PROVISIONAL) != vib_label(0.667, _PROVISIONAL, scheme="threeway")


# --------------------------------------------------------------------------
# Step-3 gate 2 decoupling (tau_purity, 2026-08): gate 2's bar is now
# scheme-dependent -- tau_purity (fixed 0.05) under "binary", tau_B
# (calibrated) under "threeway", unchanged from before this split existed.
# Gate 1 (tau_TR) and Step 2's assignment remain scheme-independent.
# --------------------------------------------------------------------------

def test_gate2_bar_scheme_dependent():
    """gate2_bar() returns tau_B for threeway (unchanged historical
    behavior) and tau_purity for binary (the new decoupled constant) -- and
    with Thresholds()'s provisional defaults the two genuinely differ
    (0.2 vs 0.05), not just in name."""
    th = _PROVISIONAL  # tau_B=0.2, tau_purity=0.05 (provisional defaults)
    assert gate2_bar(th, "threeway") == th.tau_B == 0.2
    assert gate2_bar(th, "binary") == th.tau_purity == 0.05
    assert gate2_bar(th, "threeway") != gate2_bar(th, "binary")


def test_gate2_bar_threeway_matches_pre_tau_purity_behavior():
    """Overriding tau_B (the OLD gate-2 bar, still used unconditionally by
    every pipeline before tau_purity existed) must move threeway's gate 2 --
    proving scheme="threeway" is wired to tau_B specifically, not some other
    threshold, and that tau_purity plays no role in this scheme at all."""
    th = Thresholds(tau_B=0.37, tau_purity=0.05)
    assert gate2_bar(th, "threeway") == 0.37
    # Changing tau_B must never move binary's gate 2 -- the whole point of
    # the decoupling.
    assert gate2_bar(th, "binary") == 0.05


def test_rescheme_external_label_gate2_scheme_divergence():
    """The exact case the tau_purity decoupling was introduced for: a mode
    that passes gate 1 (|score|=1.0 >= tau_TR=0.95) with V_Stretch=0.1,
    strictly between tau_purity=0.05 and tau_B=0.2 -- clean under threeway
    (0.1 <= tau_B=0.2), mixed under binary (0.1 > tau_purity=0.05). Proves
    the two schemes now genuinely diverge on STEP 3 output, not just Step 4
    internal vocabulary."""
    th = _PROVISIONAL  # tau_TR=0.95, tau_B=0.2, tau_purity=0.05
    assert th.tau_purity < 0.1 <= th.tau_B  # sanity: 0.1 is in the sensitive band
    assert rescheme_external_label("Tx", 1.0, 0.1, th, scheme="threeway") == "Tx"
    assert rescheme_external_label("Tx", 1.0, 0.1, th, scheme="binary") == "Tx*"
    # Gate 1 still governs regardless of scheme: a score below tau_TR is
    # mixed under EITHER scheme, even at V=0 (trivially passes gate 2 both
    # ways) -- gate 2's scheme-dependence never overrides gate 1.
    assert rescheme_external_label("Tx", 0.5, 0.0, th, scheme="threeway") == "Tx*"
    assert rescheme_external_label("Tx", 0.5, 0.0, th, scheme="binary") == "Tx*"


def test_rescheme_external_label_internal_label_passthrough():
    """A Step-4 internal label (S/B/SB) is returned unchanged -- gate 2's
    threshold never touches Step 4, mirroring rescheme_internal_label's
    symmetric contract for external labels."""
    th = _PROVISIONAL
    for label in (STRETCHING, BENDING, MIXED_STRETCH_BEND):
        assert rescheme_external_label(label, 1.0, 0.1, th, scheme="binary") == label
        assert rescheme_external_label(label, 1.0, 0.1, th, scheme="threeway") == label


class _StubScorer:
    """Minimal classify_all_modes()-compatible stub -- calculate_scores()
    returns the caller-supplied {"T":.., "R":.., "V":..} dict verbatim
    (the mode's own "vector" IS that dict, an opaque payload as far as
    classify_all_modes is concerned), ignoring geometry entirely. Lets Step
    3 be exercised with an EXACT, hand-picked V_Stretch (0.1) without
    needing a real molecule whose external-slot V happens to land in the
    (tau_purity, tau_B] band -- none of this repo's example molecules do:
    normal-mode T/R rows are always the trivially-perfect synthetic ideal
    references (V=0 exactly, see build_scorer_and_final), and every real
    EMIT external-slot mode either fails gate 1 outright or has V well
    above tau_B=0.2 (empirically audited when tau_purity was introduced --
    zero rows in library_scores.csv/C6H6_EMIT.csv/H2O_EMIT.csv flip under
    the new gate 2 bar)."""
    def __init__(self, v_weighting):
        self.v_weighting = v_weighting

    def calculate_scores(self, vector):
        return vector

    def score_bonds(self):
        return []

    def principal_axes(self):
        # Non-linear (smallest moment far from 0) -- n_R=3, all 6 slots present.
        import numpy as np
        return np.array([1.0, 2.0, 3.0]), None


def test_gate2_scheme_divergence_synthetic_case():
    """classify_all_modes()-level integration test (not just the pure-
    function unit tests above): a synthetic 6-mode pool, one mode per
    external slot, each perfectly aligned with its own axis (score=1.0,
    trivial unambiguous Step-2 assignment) so Step 3 is the only thing
    under test. "Ext_Tx" is hand-picked with V_Stretch=0.1 -- passes gate 2
    at tau_B=0.2 (threeway) but fails at tau_purity=0.05 (binary); every
    other slot has V=0 (clean under either bar), isolating the divergence
    to exactly the one mode constructed to exercise it.
    """
    th = _PROVISIONAL  # tau_TR=0.95, tau_B=0.2, tau_purity=0.05
    mode_specs = [
        ("Ext_Tx", {"x": 1.0, "y": 0.0, "z": 0.0}, {"x": 0.0, "y": 0.0, "z": 0.0}, 0.1),
        ("Ext_Ty", {"x": 0.0, "y": 1.0, "z": 0.0}, {"x": 0.0, "y": 0.0, "z": 0.0}, 0.0),
        ("Ext_Tz", {"x": 0.0, "y": 0.0, "z": 1.0}, {"x": 0.0, "y": 0.0, "z": 0.0}, 0.0),
        ("Ext_Rx", {"x": 0.0, "y": 0.0, "z": 0.0}, {"x": 1.0, "y": 0.0, "z": 0.0}, 0.0),
        ("Ext_Ry", {"x": 0.0, "y": 0.0, "z": 0.0}, {"x": 0.0, "y": 1.0, "z": 0.0}, 0.0),
        ("Ext_Rz", {"x": 0.0, "y": 0.0, "z": 0.0}, {"x": 0.0, "y": 0.0, "z": 1.0}, 0.0),
    ]
    final = [{"label": name, "frequency": 0.0,
              "vector": {"T": T, "R": R, "V": V}}
             for name, T, R, V in mode_specs]

    scorer = _StubScorer(th.v_weighting)
    scored_3 = {m["name"]: m for m in classify_all_modes(scorer, final, th, scheme="threeway")}
    scored_b = {m["name"]: m for m in classify_all_modes(scorer, final, th, scheme="binary")}

    assert scored_3["Ext_Tx"]["classification"] == "Tx"    # 0.1 <= tau_B=0.2 -> clean
    assert scored_b["Ext_Tx"]["classification"] == "Tx*"   # 0.1 >  tau_purity=0.05 -> mixed
    assert scored_b["Ext_Tx"]["annotation"] == f"vibration={vib_label(0.1, th, scheme='binary')}"

    # Every other slot's V=0 passes gate 2 under EITHER bar -- confirms the
    # divergence above is Ext_Tx-specific, not a wholesale scheme effect.
    for name in ("Ext_Ty", "Ext_Tz", "Ext_Rx", "Ext_Ry", "Ext_Rz"):
        axis = name.replace("Ext_", "")
        assert scored_3[name]["classification"] == axis
        assert scored_b[name]["classification"] == axis


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
