"""Golden-reference regression tests for src/projection.py (eq:emitproj).

Run from ``Github/scoring-functions/``:
    py -m pytest tests/                (with pytest)
    py tests/test_projection.py       (standalone; no pytest needed)

Pins the mass-weighting convention (IMPLEMENTATION_PLAN.md Phase-0 "Pin the
projection / mass-weighting convention" decision): both the normal-mode
reference basis Q (ideal T/R + real vibrational normal modes) and the raw
EMIT eigenvectors Theta are re-expressed in mass-weighted Cartesian
coordinates (each atom's 3-vector scaled by sqrt(mass_A)) and renormalized to
unit length before projecting (Theta_tilde = Q^T Theta). Validated during
development to <=3.2e-4 absolute agreement, on EVERY one of benzene's 36 EMIT
modes and all 9 grouped columns, against the pre-existing hand-derived
data/results/benzene_EMIT_contributions.csv ground-truth file. This test
pins the specific spot-check values IMPLEMENTATION_PLAN.md Phase 2 names:
EMIT 34-36 approx 77% translational, EMIT 2 (39% Ry) vs EMIT 9 (14% Ry).

**2026-07-07 basename update:** "benzene" -> the finalized roster basename
`C6H6` (data/mol_list_method.csv); the rename was
content-preserving (re-verified against a live run), so none of the pinned
numbers below changed.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from main import run_projection_pipeline, load_inputs, build_scorer_and_final  # noqa: E402
from src.projection import build_reference_basis_cartesian, project_emit_cartesian  # noqa: E402

TOL = 2e-3  # generous vs. the observed <=3.2e-4 max deviation from ground truth


def _rows_by_mode(df):
    return {row["Mode"]: row for _, row in df.iterrows()}


def test_benzene_emit_34_35_36_translational_fraction():
    """EMIT 34/35/36 (E1u/A2u externals) project to ~77% onto their
    respective ideal translation axis -- the manuscript's headline number."""
    df, _, _ = run_projection_pipeline("C6H6", write=False)
    rows = _rows_by_mode(df)
    assert abs(rows["EMIT 34"]["C2_Tx"] - 0.767349) < TOL
    assert abs(rows["EMIT 35"]["C2_Ty"] - 0.767349) < TOL
    assert abs(rows["EMIT 36"]["C2_Tz"] - 0.767349) < TOL


def test_benzene_emit_2_vs_9_ry_inversion():
    """EMIT 2 has more Ry projection character than EMIT 9 (39% vs 14%),
    the inverse of the |s[Ry]| score ordering -- the manuscript's
    non-monotonicity example (IMPLEMENTATION_PLAN.md Phase 2)."""
    df, _, _ = run_projection_pipeline("C6H6", write=False)
    rows = _rows_by_mode(df)
    ry2 = rows["EMIT 2"]["C2_Ry"]
    ry9 = rows["EMIT 9"]["C2_Ry"]
    assert abs(ry2 - 0.386632) < TOL
    assert abs(ry9 - 0.141039) < TOL
    assert ry2 > ry9  # projection ordering is opposite the |s[Ry]| score ordering


def test_projected_fractions_sum_to_one():
    """Every EMIT mode's grouped fractions (Parseval, since Q is
    near-orthonormal) sum to ~1 -- the sanity check project_emit() itself
    asserts internally; re-checked here at the DataFrame level."""
    df, _, _ = run_projection_pipeline("C6H6", write=False)
    cols = ["C2_Tx", "C2_Ty", "C2_Tz", "C2_Rx", "C2_Ry", "C2_Rz",
            "C2_VS", "C2_VB", "C2_VMix"]
    totals = df[cols].sum(axis=1)
    assert (totals - 1.0).abs().max() < 0.01


def test_full_projection_file_has_per_mode_detail():
    """The 'full' output carries per-individual-reference-mode detail (not
    just the grouped external/internal fractions) -- the richer
    projection-coefficients data file Phase 2 calls for."""
    df, df_full, _ = run_projection_pipeline("C6H6", write=False)
    # 6 ideal T/R slots + 30 real vibrational normal modes for benzene (3N=36).
    expected_ref_cols = {"Mode", "Eigenvalue", "Tx", "Ty", "Tz", "Rx", "Ry", "Rz"}
    assert expected_ref_cols <= set(df_full.columns)
    vib_cols = [c for c in df_full.columns if c.startswith("Vib ")]
    assert len(vib_cols) == 30
    assert len(df_full) == len(df) == 36


def test_cartesian_overlap_pathway_is_raw_and_not_orthonormal():
    """Ocart_Tx.."Ocart_Rz" are deliberately NOT squared (unlike the
    mass-weighted C2_* columns, whose squaring is only meaningful because
    that basis is near-orthonormal -- Parseval). EMIT 34's Ocart_Tx should
    equal the raw signed overlap, i.e. +sqrt of the old squared value (both
    Q and Theta are positively aligned on this pure-translation mode).

    Ocart_VS/Ocart_VB, by contrast, ARE normalized (hard-classified S/B
    |overlap| totals, L1-normalized by their S+B sum -- same magnitude-
    weighted mechanism as Vscore, see project_emit_cartesian's docstring),
    so they sum to exactly 1 and stay in [0,1] like V_Stretch even though
    this basis is non-orthonormal. Ocart_Sum (raw Tx..Rz + the normalized
    VS/VB) is still not expected to be anywhere near 1, in contrast to the
    mass-weighted C2_* fractions' ~1 (Parseval) for that same mode -- pinned
    regression values from a live run, demonstrating why the mass-weighted
    pathway is the one to trust for T/R character specifically."""
    df, _, _ = run_projection_pipeline("C6H6", write=False)
    rows = _rows_by_mode(df)
    assert abs(rows["EMIT 34"]["Ocart_Tx"] - 0.875991) < TOL
    assert abs(rows["EMIT 34"]["Ocart_Sum"] - 1.875991) < TOL
    assert abs(rows["EMIT 34"]["Ocart_VS"] + rows["EMIT 34"]["Ocart_VB"] - 1.0) < TOL
    mass_weighted_cols = ["C2_Tx", "C2_Ty", "C2_Tz", "C2_Rx", "C2_Ry", "C2_Rz",
                           "C2_VS", "C2_VB", "C2_VMix"]
    assert abs(sum(rows["EMIT 34"][c] for c in mass_weighted_cols) - 1.0) < 0.01
    # sanity: unsquared overlap re-squared recovers what C2_Tx would have been
    # under this (non-orthonormal) basis had it been squared like C2_ is.
    assert abs(rows["EMIT 34"]["Ocart_Tx"] ** 2 - 0.767360) < TOL


def test_ocart_vs_vb_sum_to_one():
    """Ocart_VS + Ocart_VB == 1 exactly for every EMIT mode -- L1
    (magnitude) normalization by the S+B total |overlap|, with no reference
    mode excluded from the sum, so the two totals always partition the full
    internal |overlap| sum with nothing dropped or double-counted."""
    df, _, _ = run_projection_pipeline("C6H6", write=False)
    totals = df["Ocart_VS"] + df["Ocart_VB"]
    assert (totals - 1.0).abs().max() < TOL


def test_ocart_vs_vb_l1_no_exclusion():
    """No reference normal mode is dropped from the Ocart_VS/VB sum --
    every internal mode's |overlap| counts toward whichever bucket
    classifier.vib_label(s[V_S]) assigns it, including the 10/30 benzene
    real normal modes with "medium" s[V_S] in [0.20, 0.80] (two other
    schemes were tried and reverted: excluding those medium-V modes
    entirely, and a continuous v_m-proportional split instead of a hard
    S/B bucket -- both changed the numbers only slightly versus this plain
    L1 version, so the simpler version is what's kept). Pinned regression
    value for EMIT 34."""
    df, _, _ = run_projection_pipeline("C6H6", write=False)
    rows = _rows_by_mode(df)
    assert abs(rows["EMIT 34"]["Ocart_VS"] - 0.596629) < TOL
    assert abs(rows["EMIT 34"]["Ocart_VB"] - 0.403371) < TOL


def test_full_cartesian_projection_has_per_mode_detail():
    """The Cartesian-overlap 'full' output has the same per-individual-
    reference-mode shape as the mass-weighted one (build_reference_basis_cartesian
    / project_emit_cartesian, tested directly since run_projection_pipeline
    only writes the cartesian file rather than returning it)."""
    raw_n, _ = load_inputs("C6H6", "normal")
    scorer_n, final_n = build_scorer_and_final(raw_n, "normal")
    raw_e, _ = load_inputs("C6H6", "emit")
    _, final_e = build_scorer_and_final(raw_e, "emit")

    ref = build_reference_basis_cartesian(scorer_n, final_n)
    rows, full_rows = project_emit_cartesian(ref, final_e)

    assert len(rows) == len(full_rows) == 36
    expected_row_cols = {"Mode", "Eigenvalue", "Ocart_Tx", "Ocart_Ty", "Ocart_Tz",
                          "Ocart_Rx", "Ocart_Ry", "Ocart_Rz", "Ocart_VS",
                          "Ocart_VB", "Ocart_Sum"}
    assert expected_row_cols <= set(rows[0].keys())
    expected_ref_cols = {"Mode", "Eigenvalue", "Tx", "Ty", "Tz", "Rx", "Ry", "Rz"}
    assert expected_ref_cols <= set(full_rows[0].keys())
    vib_keys = [k for k in full_rows[0] if k.startswith("Vib ")]
    assert len(vib_keys) == 30


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
