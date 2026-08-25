"""Regression tests for src/flag_validation.py (Phase-6 systematic flag
precision/recall check -- IMPLEMENTATION_PLAN.md Phase 6, item 1).

Run from ``Github/scoring-functions/``:
    py -m pytest tests/                (with pytest)
    py tests/test_flag_validation.py  (standalone; no pytest needed)

**2026-08-25 RESTRUCTURING NOTE** (see src/flag_validation.py's own module
docstring for the full mechanism): Step 3 (T/R identification) no longer has
ANY purity gate, so there is no more "mixed-external flag" for this module
to validate -- `is_mixed_external()` on a freshly-produced `tr_label` is
always False by construction. The pre-2026-08-25 headline numbers (TP=5,
FP=0, FN=12, TN=19, precision 1.0, recall 5/17~0.294) are RETIRED, not
repinned to a stale historical value -- the correct, current numbers are
degenerate (TP=0 everywhere, precision undefined/NaN, recall 0.0), and that
degeneracy is now the thing being pinned/tested below, not a bug. The
geometry-backed library external (T/R) rows still correctly show FP=0
(trivial, ground truth always CLEAN by Eckart-Sayvetz completeness) --
UNCHANGED from before, since that count was already 0 either way.

**2026-07-07 basename update:** benzene is now resolved through
library_ingest.resolve_log_basename("C6H6", ...) inside
src/flag_validation.py itself (no test-side change needed for the benzene
tests below -- content-preserving rename, numbers unchanged). The library
external-reference population grew from 25 to 72 geometry-backed molecules
(427 external rows) now that the full mol_list_method.csv roster is
Gaussian-direct, including OBr4 once its data/gjf .com connectivity gap
was fixed (see IMPLEMENTATION_PLAN.md's 2026-07-07 RESUME HERE) --
test_library_external_references_never_false_positive's counts are updated
accordingly.

**2026-07-09 OH4/OF4 exclusion:** OH4 and OF4 removed from the roster (not
genuine stationary points at this project's MP2/3-21G level -- imaginary/
negative frequencies -- so their normal modes cannot be validly compared to
the TeH4 ideal see-saw template). The library external-reference population
shrank from 72 to 70 geometry-backed molecules (427 to 415 external rows,
-12 = 2 molecules x (n_T=3 + n_R=3) for non-linear AB4) --
test_library_external_references_never_false_positive's counts are updated
accordingly.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from src.flag_validation import (                                    # noqa: E402
    benzene_emit_flag_confusion, library_external_flag_confusion,
    ground_truth_label, GT_EXT_LO, GT_EXT_HI,
)


def test_benzene_emit_confusion_counts_pinned():
    """2026-08-25: the flag no longer exists (Step 3 has no purity gate), so
    predicted_positive is False for all 36 modes by construction -- TP=0,
    FP=0 (never fires, trivially), FN=17 (every genuinely-mixed mode is
    missed, since nothing is ever flagged), TN=19, precision undefined
    (0/0 -> NaN), recall 0.0. This is the current, correct, DEGENERATE
    result, not a repin of the old headline number -- see this module's
    RESTRUCTURING NOTE above and src/flag_validation.py's own docstring."""
    detail, stats = benzene_emit_flag_confusion()
    assert stats["n"] == 36
    assert stats["TP"] == 0
    assert stats["FP"] == 0
    assert stats["FN"] == 17
    assert stats["TN"] == 19
    assert stats["precision"] != stats["precision"]  # NaN
    assert stats["recall"] == 0.0
    assert len(detail) == 36


def test_benzene_emit_36_is_a_false_negative_not_a_true_negative():
    """The documented Decision-X blind spot must show up as FN here (the
    classifier calls it CLEAN_TRANSLATION, but projection shows ~76.7%
    translational + genuine out-of-plane-bending residual, i.e. ground truth
    MIXED) -- this is the expected, intentionally-reported limitation, not a
    bug to fix."""
    detail, _ = benzene_emit_flag_confusion()
    row = detail[detail["Mode"] == "EMIT 36"].iloc[0]
    assert row["classifier_label"] == "Tz"
    assert row["ground_truth"] == "MIXED"
    assert row["cell"] == "FN"


def test_benzene_emit_34_35_are_false_negatives():
    """2026-08-25: EMIT 34/35 (Tx/Ty=1, s[V_S]=0.667/0.577) win their T-axis
    slot cleanly (bare "Tx"/"Ty", no more starred label -- Step 3 has no
    purity gate to fail) and are STILL genuinely externally-mixed by
    projection (~76.7% translational, rest vibrational) -- but since there
    is no more flag mechanism at all, they are now FALSE NEGATIVES (missed),
    not true positives. This is the direct, mechanical consequence of
    retiring the purity gate, not a regression to fix."""
    detail, _ = benzene_emit_flag_confusion()
    expected = {"EMIT 34": "Tx", "EMIT 35": "Ty"}
    for name, label in expected.items():
        row = detail[detail["Mode"] == name].iloc[0]
        assert row["classifier_label"] == label
        assert row["cell"] == "FN"


def test_emit_2_vs_9_ry_inversion_still_visible_in_projection_and_assignment():
    """EMIT 2 has MORE genuine Ry projection character (38.7%) than EMIT 9
    (14.1%) -- the documented score/projection ranking inversion -- yet
    EMIT 9 wins the Ry slot (its |s[Ry]| SCORE is larger), while EMIT 2 does
    not win any slot at all. This underlying ranking inversion (a fact about
    src/projection.py's projection vs. Step 3's own score-based assignment)
    is UNCHANGED by the 2026-08-25 restructuring -- what changed is that
    BOTH modes are now FN (there is no more flag for EMIT 9 to be a TP of),
    since Step 3 no longer has a purity gate to pass or fail. See
    test_benzene_emit_confusion_counts_pinned for the mechanism."""
    detail, _ = benzene_emit_flag_confusion()
    row2 = detail[detail["Mode"] == "EMIT 2"].iloc[0]
    row9 = detail[detail["Mode"] == "EMIT 9"].iloc[0]
    assert row2["M_ext"] > row9["M_ext"]
    assert row2["cell"] == "FN"
    assert row9["cell"] == "FN"
    assert row2["classifier_label"] == "B"    # never won a slot
    assert row9["classifier_label"] == "Ry"   # won the Ry slot, but cleanly (no flag exists to fire)


def test_ground_truth_label_thresholds():
    assert GT_EXT_LO == 0.05
    assert GT_EXT_HI == 0.95
    assert ground_truth_label(0.0) == "CLEAN"
    assert ground_truth_label(0.03) == "CLEAN"
    assert ground_truth_label(0.5) == "MIXED"
    assert ground_truth_label(0.9674) == "CLEAN"
    assert ground_truth_label(1.0) == "CLEAN"


def test_library_external_references_never_false_positive():
    """The 85 geometry-backed library molecules' REAL normal-mode T/R
    references (506 rows total) are exact-by-construction Eckart-Sayvetz
    references (ground truth always CLEAN); this checks -- rather than
    assumes -- that the classifier never flags a single one of them
    MIXED_EXTERNAL_WITH_VIBRATION (FP=0), the much easier degenerate case
    named in the task's parenthetical. `library_external_flag_confusion()`
    filters on `kind == 'external'` only, not on the roster's `mol_type`/
    `ideal` column, so this count spans every roster molecule regardless of
    category (ideal/non-ideal/multi-centre/test) -- it's about T/R
    reference correctness, not the ideal/non-ideal calibration scope.

    Roster history: 68/404 (2026-07-09 OH4/OF4 removal, 2026-07-23 SnO2/FH3
    removal + CO2 re-tag) held until 2026-08-11, when 9 'test'-category
    transferability molecules were added (none linear, so +9*6=54 external
    rows) and library_scores.csv was regenerated to include them: 68 -> 77
    molecules, 404 -> 458 external rows. **2026-08-14:** 8 more 'test'-category
    molecules added (none linear, so +8*6=48 external rows; H2O also moved
    'non-ideal' -> 'test' but keeps contributing its own 6 external rows
    either way -- this test spans every roster molecule regardless of
    category): 77 -> 85 molecules, 458 -> 506 external rows."""
    ext, stats = library_external_flag_confusion()
    assert stats["n"] == 506
    assert stats["n_molecules"] == 85
    assert stats["FP"] == 0
    assert stats["TP"] == 0
    assert stats["FN"] == 0
    assert stats["TN"] == 506


if __name__ == "__main__":
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
