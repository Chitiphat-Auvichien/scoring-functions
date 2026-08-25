"""Golden-reference regression tests for src/classifier.py (Algorithm 1).

Run from ``Github/scoring-functions/``:
    py -m pytest tests/                (with pytest)
    py tests/test_classifier.py       (standalone; no pytest needed)

Pins the label-level validation targets from IMPLEMENTATION_PLAN.md: water's
6 constructed T/R references win their own Step-3 slot, water's 3 real
vibrations split BENDING/STRETCHING/STRETCHING, and benzene EMIT 34/35
(stretching-character externals) vs EMIT 36 (the documented out-of-plane-
bending blind spot -- see IMPLEMENTATION_PLAN.md Decision X) all resolve as
designed.

**2026-08-25 restructuring** (see src/classifier.py's module docstring): the
old two-gate purity test (gate2_bar/Thresholds.tau_purity/the starred "Tx*"
mixed-external label/rescheme_external_label) is RETIRED. Step 2 (vib_label)
now runs UNCONDITIONALLY on every mode; Step 3 (tr_label/tr_score) is
OPTIONAL (identify_tr=True by default) and applies NO purity gate at all --
every Hungarian-assignment winner gets the bare slot name and its own raw
signed score, full stop. The gate2-specific tests that used to live here
(test_gate2_bar_scheme_dependent, test_gate2_bar_threeway_matches_pre_tau_
purity_behavior, test_rescheme_external_label_gate2_scheme_divergence,
test_rescheme_external_label_internal_label_passthrough,
test_gate2_scheme_divergence_synthetic_case, and the _StubScorer class that
only they used) are DELETED, not repinned -- their entire subject matter no
longer exists. New tests below cover the restructured Step 2/Step 3 split.

Threshold policy (decided Phase 3, 2026-07, unchanged by the 2026-08-25
restructuring): this file explicitly pins Thresholds() -- the hardcoded
PROVISIONAL constants (tau_TR=0.95, tau_S=0.9, tau_B=0.2, tau_SB=0.42) --
rather than relying on classify_all_modes()'s default (which now auto-loads
the Phase-3 CALIBRATED values via Thresholds.calibrated() when
data/results/thresholds.json exists). This keeps these regression goldens
fixed and reproducible even if a future recalibration (e.g. an expanded
library) changes thresholds.json's numbers. tau_TR is DIAGNOSTIC-ONLY as of
2026-08-25 (Step 3 no longer gates on it) -- kept pinned here anyway for
provenance, since nothing in this file's assertions depends on its value
anymore.

**2026-07-07 basename update:** "water"/"benzene" -> the finalized roster
basenames `H2O`/`C6H6` (data/mol_list_method.csv).
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from main import load_inputs, build_scorer_and_final              # noqa: E402
from src.classifier import (                                       # noqa: E402
    classify_all_modes, classify_to_rows, Thresholds,
    vib_label, vib_label_binary, predicted_category, predicted_category_column,
    STRETCHING, BENDING, MIXED_STRETCH_BEND,
)

TOL = 5e-4  # 3 decimal places

# Explicit provisional thresholds -- see module docstring for why this file
# does not rely on classify_all_modes()'s calibrated-by-default resolution.
_PROVISIONAL = Thresholds()


def _classify(mol, mode_type, scheme="threeway", identify_tr=True):
    """`scheme` pinned to "threeway" by default: every golden target below
    (water B/S, benzene EMIT 34/35/36) was derived/pinned under the
    three-way scheme, so this keeps them fixed and reproducible even though
    classify_all_modes()'s own default flipped to "binary" (the 2026-08
    paper-standard switch) -- see
    test_classify_all_modes_default_scheme_is_now_binary below for the
    dedicated check that the new global default is genuinely binary."""
    raw, _ = load_inputs(mol, mode_type, os.path.join(ROOT, "data"))
    scorer, final = build_scorer_and_final(raw, mode_type)
    scored = classify_all_modes(scorer, final, _PROVISIONAL, scheme=scheme,
                                 identify_tr=identify_tr)
    return {m["name"]: m for m in scored}


def test_water_normal_externals_win_own_slot():
    """The 6 ideal T/R reference modes each win their OWN Step-3 slot (the
    bare slot name itself, e.g. "Tx") with |tr_score| very close to 1 --
    exact Eckart-Sayvetz completeness, no purity gate involved anymore."""
    t = _classify("H2O", "normal")
    for lbl in ("Tx", "Ty", "Tz", "Rx", "Ry", "Rz"):
        assert t[lbl]["tr_label"] == lbl, t[lbl]["tr_label"]
        assert abs(abs(t[lbl]["tr_score"]) - 1.0) < 1e-6


def test_water_normal_vibrations():
    """Vib1 = bend -> BENDING; Vib2/Vib3 = stretches -> STRETCHING (Step 2,
    vib_label), with |s_AB| summing to V (s_AB is signed; V sums
    magnitudes). None of water's 3 real vibrational modes win a Step-3 slot
    (the ideal T/R references always do instead)."""
    t = _classify("H2O", "normal")
    assert t["Vib 1"]["vib_label"] == BENDING
    assert t["Vib 1"]["tr_label"] is None
    assert len(t["Vib 1"]["bonds"]) == 2  # bonds are attached for every mode, not just S/SB

    for name in ("Vib 2", "Vib 3"):
        assert t[name]["vib_label"] == STRETCHING
        assert t[name]["tr_label"] is None
        bonds = t[name]["bonds"]
        assert len(bonds) == 2  # water has 2 O-H bonds
        assert abs(sum(abs(b["s_AB"]) for b in bonds) - t[name]["V"]) < 1e-6


def test_step2_runs_unconditionally_including_synthetic_tr_rows():
    """Step 2 (vib_label) is computed for EVERY mode in the pool, including
    the 6 synthetic ideal T/R reference rows -- which trivially get "B"
    since their V_Stretch is exactly 0. This is expected, not a bug: a
    downstream consumer must treat tr_label as authoritative for T/R
    identity, never vib_label (see src/classifier.py's module docstring)."""
    t = _classify("H2O", "normal")
    for lbl in ("Tx", "Ty", "Tz", "Rx", "Ry", "Rz"):
        assert abs(t[lbl]["V"]) < 1e-9
        assert t[lbl]["vib_label"] == BENDING  # trivial consequence of V=0, not meaningful T/R info


def test_identify_tr_false_leaves_tr_fields_blank():
    """identify_tr=False skips Step 3 entirely: every mode's tr_label/
    tr_score stay None, but vib_label is still fully populated (Step 2 is
    unconditional and independent of identify_tr)."""
    t = _classify("H2O", "normal", identify_tr=False)
    for m in t.values():
        assert m["tr_label"] is None
        assert m["tr_score"] is None
        assert m["vib_label"] in (STRETCHING, BENDING, MIXED_STRETCH_BEND)


def test_benzene_emit_34_35_36_spot_check():
    """The Decision-X benzene EMIT 34/35/36 spot check under the
    restructured (2026-08-25) algorithm, values confirmed against a live
    rerun (not assumed from the exploration pass that originally proposed
    them -- s[T]=1.0 exactly for all three, matching the "perfect axis
    alignment" framing in IMPLEMENTATION_PLAN.md's Decision X; V_Stretch
    0.667/0.577/0 is a SEPARATE quantity from tr_score, not to be confused
    with it).

    EMIT 36 (A2u, out-of-plane bending, the documented blind spot): clean
    tr_label="Tz", tr_score very close to 1.0, vib_label="B" (V_Stretch=0 --
    Step 2 cannot see the out-of-plane bending residual, same blind spot as
    before, just re-expressed without a purity flag to fail).

    EMIT 34/35 (E1u, in-plane, stretching-character): both win a T-axis
    slot (tr_label in {"Tx","Ty"}) with |tr_score|~=1.0 -- the SAME argmax
    winners as before 2026-08-25 (Step 3's assignment mechanism itself is
    unchanged, only the removed gate), no more asterisk. Under the default
    binary scheme with tau_SB=0.42 (both 0.667 and 0.577 are >= 0.42),
    vib_label="S" for both -- a real, checkable behavior change from the
    pre-2026-08-25 threeway "SB" framing (V_Stretch itself, tab:water's own
    ground truth, is completely unaffected -- this is a label-layer-only
    change).
    """
    t = _classify("C6H6", "emit", scheme="binary")

    e36 = t["EMIT 36"]
    assert e36["tr_label"] == "Tz"
    assert abs(abs(e36["tr_score"]) - 1.0) < 1e-6
    assert abs(e36["V"]) < TOL
    assert e36["vib_label"] == BENDING

    e34, e35 = t["EMIT 34"], t["EMIT 35"]
    assert e34["tr_label"] in ("Tx", "Ty")
    assert e35["tr_label"] in ("Tx", "Ty")
    assert e34["tr_label"] != e35["tr_label"]  # one-to-one assignment, distinct slots
    assert abs(abs(e34["tr_score"]) - 1.0) < 1e-6
    assert abs(abs(e35["tr_score"]) - 1.0) < 1e-6
    assert abs(e34["V"] - 0.6666662706312129) < TOL
    assert abs(e35["V"] - 0.5773504978406205) < TOL

    tau_SB = Thresholds().tau_SB  # provisional default 0.42
    assert e34["V"] >= tau_SB and e35["V"] >= tau_SB
    assert e34["vib_label"] == STRETCHING
    assert e35["vib_label"] == STRETCHING


def test_co2_linear_no_spurious_onaxis_mode():
    """Linear-molecule candidate pool must have exactly 3N modes, not 3N+1.

    Regression for a bug found by formula-auditor: construct_R() always builds
    3 ideal rotation references, but for a linear molecule MIT() places the
    molecular (smallest-moment) axis on the new X axis, so the ideal "Rx"
    reference is an all-zero vector (n_R=2 excludes it, spec). Left in the
    pool, it fell through Step 3 (no slot claims an all-zero mode) into
    Step 2's bare vib_label and was mislabeled BENDING (V=0 <= tau_B).
    build_scorer_and_final() must drop it before scoring, and the on-axis
    slot must simply be absent (not present-and-clean, not
    present-and-mislabeled).
    """
    t = _classify("CO2", "normal")
    assert len(t) == 9  # 3N for a 3-atom linear molecule, not 3N+1
    assert "Rx" not in t
    assert t["Ry"]["tr_label"] == "Ry"
    assert t["Rz"]["tr_label"] == "Rz"
    # the 2 bends and 2 stretches must be the REAL vibrational modes, not a
    # spurious all-zero placeholder
    vib_labels = {name: m["vib_label"] for name, m in t.items() if name.startswith("Vib")}
    assert sorted(vib_labels.values()) == [BENDING, BENDING, STRETCHING, STRETCHING]


def test_classify_to_rows_shape():
    """classify_to_rows() produces the documented CSV columns for both mol/mode_type combos,
    with one 's_AB[iLabel-jLabel]' column per bond (not a joined string) -- and every row
    (every mode) carries the same bond columns, since a molecule's bond list doesn't
    depend on which mode is being scored."""
    base_cols = {"Mode", "Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "V_Stretch",
                 "Mu", "K", "Irrep", "vib_label", "tr_label", "tr_score"}
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
    for EVERY mode in the pool (never 'SB') -- and specifically diverge on
    these two modes (Vib 13 -> B, Vib 14 -> S), confirming the two schemes
    are genuinely different, not just synonyms. Step 3's tr_label/tr_score
    are identical between the two calls (scheme never touches Step 3).
    """
    raw, _ = load_inputs("C6H6", "normal", os.path.join(ROOT, "data"))
    scorer, final = build_scorer_and_final(raw, "normal")
    # tau_SB pinned explicitly at 0.5 (not the bare Thresholds() default,
    # which is 0.42 -- see Thresholds' docstring). This test is about
    # scheme-divergence mechanics on a straddling pair, not about the
    # current canonical tau_SB value, so it stays independent of that
    # constant rather than being repinned every time tau_SB moves.
    thresholds = Thresholds(tau_SB=0.5)  # provisional tau_S=0.9, tau_B=0.2, tau_SB=0.5

    scored_3 = {m["name"]: m for m in classify_all_modes(scorer, final, thresholds, scheme="threeway")}
    scored_b = {m["name"]: m for m in classify_all_modes(scorer, final, thresholds, scheme="binary")}

    # Both real modes' V_Stretch sit strictly inside the threeway SB band.
    assert 0.2 < scored_3["Vib 13"]["V"] < 0.9
    assert 0.2 < scored_3["Vib 14"]["V"] < 0.9
    assert scored_3["Vib 13"]["vib_label"] == MIXED_STRETCH_BEND
    assert scored_3["Vib 14"]["vib_label"] == MIXED_STRETCH_BEND

    # Binary scheme diverges: same V_Stretch, different verdict either side
    # of tau_SB=0.5 (Vib 13's V~0.451 < 0.5 -> BENDING; Vib 14's V~0.509 >= 0.5 -> STRETCHING).
    assert scored_b["Vib 13"]["vib_label"] == BENDING
    assert scored_b["Vib 14"]["vib_label"] == STRETCHING

    # No mode in the WHOLE pool is ever labeled "SB" under the binary scheme.
    for m in scored_b.values():
        assert m["vib_label"] != MIXED_STRETCH_BEND

    # Step 3 (tr_label/tr_score) is untouched by scheme -- identical winners
    # and scores under either call.
    for name in scored_3:
        assert scored_3[name]["tr_label"] == scored_b[name]["tr_label"], name
        assert scored_3[name]["tr_score"] == scored_b[name]["tr_score"], name


def test_classify_all_modes_default_scheme_is_now_binary():
    """2026-08 paper-standard switch: classify_all_modes()'s own default
    (no `scheme` kwarg at all, unlike `_classify()` above which pins
    "threeway" explicitly) must now be BINARY, not threeway. Benzene EMIT
    34/35 (V_Stretch 0.667/0.577, both >= tau_SB=0.42) are STRETCHING under
    binary, never MIXED_STRETCH_BEND. tr_label/tr_score are identical to the
    threeway-scheme result (scheme never touches Step 3)."""
    raw, _ = load_inputs("C6H6", "emit", os.path.join(ROOT, "data"))
    scorer, final = build_scorer_and_final(raw, "emit")
    scored_default = {m["name"]: m for m in classify_all_modes(scorer, final, _PROVISIONAL)}  # no scheme=
    scored_threeway = {m["name"]: m for m in classify_all_modes(scorer, final, _PROVISIONAL, scheme="threeway")}
    assert scored_default["EMIT 34"]["vib_label"] == STRETCHING
    assert scored_default["EMIT 35"]["vib_label"] == STRETCHING
    assert scored_default["EMIT 34"]["tr_label"] == scored_threeway["EMIT 34"]["tr_label"]
    assert scored_default["EMIT 34"]["tr_score"] == scored_threeway["EMIT 34"]["tr_score"]
    # Confirms this genuinely differs from the threeway-pinned result above,
    # not just a coincidentally-identical label.
    assert vib_label(0.667, _PROVISIONAL) != vib_label(0.667, _PROVISIONAL, scheme="threeway")


def test_predicted_category_helper():
    """predicted_category()/predicted_category_column() implement the
    canonical rule the whole codebase now shares: tr_label if non-empty,
    else vib_label. Empty string/None/NaN-like tr_label all fall through to
    vib_label."""
    assert predicted_category("Tx", "S") == "Tx"
    assert predicted_category("", "S") == "S"
    assert predicted_category(None, "B") == "B"
    assert predicted_category(float("nan"), "B") == "B"  # not a str -> falls through

    import pandas as pd
    tr_col = pd.Series(["Tx", "", None, float("nan")])
    vib_col = pd.Series(["S", "B", "S", "B"])
    out = predicted_category_column(tr_col, vib_col)
    assert list(out) == ["Tx", "B", "S", "B"]


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
