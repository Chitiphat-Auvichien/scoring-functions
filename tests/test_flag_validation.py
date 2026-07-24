"""Regression tests for src/flag_validation.py (Phase-6 systematic flag
precision/recall check -- IMPLEMENTATION_PLAN.md Phase 6, item 1).

Run from ``Github/scoring-functions/``:
    py -m pytest tests/                (with pytest)
    py tests/test_flag_validation.py  (standalone; no pytest needed)

Pins the exact confusion counts computed this session (see
src/flag_validation.py's module docstring for the ground-truth criterion and
its rationale): benzene's 36 EMIT modes give TP=5, FP=0, FN=12, TN=19
(precision 1.0, recall 5/17 ~ 0.294), and the geometry-backed library
external (T/R) rows give FP=0 (trivial, ground truth always CLEAN by
Eckart-Sayvetz completeness).

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
    """The headline systematic result: precision 1.0 (flag never fires on a
    genuinely clean mode), recall 5/17~0.294 (it misses most genuinely mixed
    modes) -- across ALL 36 modes, not just the prior 3-mode spot-check."""
    detail, stats = benzene_emit_flag_confusion()
    assert stats["n"] == 36
    assert stats["TP"] == 5
    assert stats["FP"] == 0
    assert stats["FN"] == 12
    assert stats["TN"] == 19
    assert stats["precision"] == 1.0
    assert abs(stats["recall"] - 5 / 17) < 1e-9
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


def test_benzene_emit_34_35_are_true_positives():
    """EMIT 34/35 (Tx/Ty=1, s[V_S]=0.667/0.577) are correctly flagged with the
    axis-specific mixed-external label ("Tx*"/"Ty*") by the classifier AND are
    genuinely externally-mixed by projection (~76.7% translational, rest
    vibrational) -- a true positive, matching the prior anecdotal spot-check."""
    detail, _ = benzene_emit_flag_confusion()
    expected = {"EMIT 34": "Tx*", "EMIT 35": "Ty*"}
    for name, label in expected.items():
        row = detail[detail["Mode"] == name].iloc[0]
        assert row["classifier_label"] == label
        assert row["cell"] == "TP"


def test_emit_2_vs_9_ry_inversion_produces_opposite_confusion_cells():
    """EMIT 2 has MORE genuine Ry projection character (38.7%) than EMIT 9
    (14.1%) -- the documented score/projection ranking inversion -- yet
    EMIT 9 wins the Ry slot (its |s[Ry]| SCORE is larger) and is flagged
    (TP), while EMIT 2, despite being the more genuinely mixed mode, is
    missed (FN). This is direct evidence the low recall is a systematic
    consequence of Step 2's one-to-one assignment, not an isolated case."""
    detail, _ = benzene_emit_flag_confusion()
    row2 = detail[detail["Mode"] == "EMIT 2"].iloc[0]
    row9 = detail[detail["Mode"] == "EMIT 9"].iloc[0]
    assert row2["M_ext"] > row9["M_ext"]
    assert row2["cell"] == "FN"
    assert row9["cell"] == "TP"


def test_ground_truth_label_thresholds():
    assert GT_EXT_LO == 0.05
    assert GT_EXT_HI == 0.95
    assert ground_truth_label(0.0) == "CLEAN"
    assert ground_truth_label(0.03) == "CLEAN"
    assert ground_truth_label(0.5) == "MIXED"
    assert ground_truth_label(0.9674) == "CLEAN"
    assert ground_truth_label(1.0) == "CLEAN"


def test_library_external_references_never_false_positive():
    """The 68 geometry-backed library molecules' REAL normal-mode T/R
    references (404 rows total) are exact-by-construction Eckart-Sayvetz
    references (ground truth always CLEAN); this checks -- rather than
    assumes -- that the classifier never flags a single one of them
    MIXED_EXTERNAL_WITH_VIBRATION (FP=0), the much easier degenerate case
    named in the task's parenthetical.

    Roster history behind 68/404 (2026-07-09: OH4/OF4 removed -- neither is
    a genuine stationary point at this project's MP2/3-21G level
    (imaginary/negative frequencies), so their normal modes cannot be
    validly compared to the TeH4 ideal see-saw template; 72 -> 70 molecules,
    427 -> 415 external rows, -12 for OH4+OF4's 2x(n_T=3+n_R=3) non-linear-
    AB4 references). 2026-07-23 (commit `8a12c32`) removed `SnO2` (linear,
    n_T=3+n_R=2=5 external rows) and `FH3` (T-shape, non-linear, 6 external
    rows) from the roster entirely (both had negative/imaginary
    frequencies), and separately re-tagged `CO2` `ideal`: 70 -> 68
    molecules, 415 -> 404 external rows (-11, not -12, since SnO2 is linear
    and only loses 5 rows, not 6). The CO2 ideal/non-ideal tag has no effect
    on this count -- `library_external_flag_confusion()` filters on
    `kind == 'external'` only, not on the `ideal` column -- so 2026-07-24's
    CO2 -> non-ideal correction (this session, undoing 8a12c32's re-tag to
    match linear siblings CS2/CSe2/CTe2) leaves 68/404 unchanged; verified
    directly against the regenerated `library_scores.csv` rather than
    assumed."""
    ext, stats = library_external_flag_confusion()
    assert stats["n"] == 404
    assert stats["n_molecules"] == 68
    assert stats["FP"] == 0
    assert stats["TP"] == 0
    assert stats["FN"] == 0
    assert stats["TN"] == 404


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
