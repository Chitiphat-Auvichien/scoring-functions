"""Regression tests for src/benzene_validation.py -- the manuscript's PRIMARY
classification-vs-reference validation (benzene's real NORMAL modes vs.
literature/group-theory ref_label; see JCC/Scoring_Manuscript_Plan_2026-07-01.pdf
for the mandated Results & Discussion positioning, and the module's own
docstring for why this is non-circular, unlike the excluded EMIT confusion
matrix in src/flag_validation.py).

Run from ``Github/scoring-functions/``:
    py -m pytest tests/                          (with pytest)
    py tests/test_benzene_validation.py         (standalone; no pytest needed)

Design note (mirrors tests/test_calibrate.py / tests/test_excel_ingest.py):
reads the already-committed ``data/results/library_scores.csv`` golden
directly (fast, plain pandas) rather than re-running the ~1-minute Excel
ingest; both `benzene_validation` functions are pure re-derivations from that
one CSV, so this is a faithful regression of the module's actual logic, not
just of the CSV.
"""
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import pandas as pd                                                # noqa: E402

from src.benzene_validation import (                               # noqa: E402
    benzene_normal_reference_detail, benzene_normal_reference_summary,
    benzene_mixed_bond_diagnostic, MOLECULE,
)

LIB_CSV = os.path.join(ROOT, "data", "results", "library_scores.csv")


def _lib():
    return pd.read_csv(LIB_CSV)


def test_benzene_normal_detail_covers_all_36_modes():
    detail = benzene_normal_reference_detail(_lib())
    assert len(detail) == 36  # 6 external + 30 internal, all literature-labeled
    assert set(detail["ref_label"].unique()) == {"translation", "rotation", "stretch", "bend"}
    assert detail["predicted_label"].notna().all()


def test_benzene_normal_summary_matches_ad_hoc_session_numbers():
    """Pins the exact headline numbers reported this session (Task A):
    6/6 external correct; 7/7 literature stretch modes recalled; 18/23
    literature bend modes recalled; ZERO crossings into the opposite clean
    category in either direction."""
    detail = benzene_normal_reference_detail(_lib())
    summary = benzene_normal_reference_summary(detail)
    by_ref = summary.set_index("ref_label")

    assert int(by_ref.loc["translation", "n"]) == 3
    assert int(by_ref.loc["translation", "n_correct"]) == 3
    assert by_ref.loc["translation", "recall"] == 1.0

    assert int(by_ref.loc["rotation", "n"]) == 3
    assert int(by_ref.loc["rotation", "n_correct"]) == 3
    assert by_ref.loc["rotation", "recall"] == 1.0

    assert int(by_ref.loc["stretch", "n"]) == 7
    assert int(by_ref.loc["stretch", "n_correct"]) == 7
    assert by_ref.loc["stretch", "recall"] == 1.0

    assert int(by_ref.loc["bend", "n"]) == 23
    assert int(by_ref.loc["bend", "n_correct"]) == 18
    assert abs(by_ref.loc["bend", "recall"] - 18 / 23) < 1e-9
    assert int(by_ref.loc["bend", "n_migrated_to_mixed"]) == 5

    # The zero-crossings claim, computed, not eyeballed.
    assert int(by_ref.loc["stretch", "n_crossed_opposite"]) == 0
    assert int(by_ref.loc["bend", "n_crossed_opposite"]) == 0


def test_benzene_normal_detail_raises_on_incomplete_merge():
    """Fail-loud check: if a ref_label row lacks a predicted_label (merge
    never happened), the function must raise, not silently validate against
    a partial dataset."""
    lib_df = _lib().copy()
    mask = (lib_df["molecule"] == MOLECULE) & (lib_df["ref_label"] == "bend")
    idx = lib_df[mask].index[0]
    lib_df.loc[idx, "predicted_label"] = None
    try:
        benzene_normal_reference_detail(lib_df)
        assert False, "expected ValueError for an incomplete geometry merge"
    except ValueError:
        pass


def test_benzene_normal_detail_raises_if_molecule_absent():
    lib_df = _lib()
    lib_df = lib_df[lib_df["molecule"] != MOLECULE]
    try:
        benzene_normal_reference_detail(lib_df)
        assert False, "expected ValueError when C6H6 has no ref_label rows"
    except ValueError:
        pass


def test_benzene_mixed_bond_diagnostic_finds_exactly_the_five_known_modes():
    bond_detail, pairs = benzene_mixed_bond_diagnostic(_lib())
    assert sorted(bond_detail["mode_index"].tolist()) == [13, 14, 19, 23, 24]


def test_benzene_mixed_bond_diagnostic_ch_contributions_are_near_zero():
    """The 5 mixed modes' C-H bond contributions are ~0 (0.0005-0.0025 per
    bond in the ad hoc session finding) -- essentially all of V_Stretch sits
    on the 6 C-C ring bonds, not noise scattered across all 12 bonds."""
    bond_detail, _ = benzene_mixed_bond_diagnostic(_lib())
    assert (bond_detail["ch_total"] < 0.02).all(), bond_detail[["mode_index", "ch_total"]]
    assert (bond_detail["cc_fraction_of_V"] > 0.95).all(), bond_detail[["mode_index", "cc_fraction_of_V"]]


def test_benzene_mixed_bond_diagnostic_mode_19_uniform_across_ring():
    """Mode 19 (1319.2678 cm-1) shows a perfectly uniform 6-fold-symmetric
    C-C contribution -- computed as a low coefficient of variation across
    the 6 ring bonds, not eyeballed."""
    bond_detail, _ = benzene_mixed_bond_diagnostic(_lib())
    row = bond_detail[bond_detail["mode_index"] == 19].iloc[0]
    cc_vals = [row[f"s_AB[{b}]"] for b in
               ("C1-C2", "C2-C3", "C3-C4", "C4-C5", "C5-C6", "C1-C6")]
    cc_vals = pd.Series(cc_vals)
    cv = cc_vals.std() / cc_vals.mean()
    assert cv < 0.02, cv  # essentially flat across all 6 ring bonds


def test_benzene_mixed_bond_diagnostic_detects_near_degenerate_complementary_pairs():
    """The (13,14) and (23,24) near-degenerate pairs (freq splitting
    <0.03 cm-1) must be detected within a 1 cm-1 tolerance and show a
    strong NEGATIVE (complementary) C-C bond-pattern correlation -- turning
    the "complementary alternating pattern" claim from an eyeballed
    observation into a computed fact."""
    bond_detail, pairs = benzene_mixed_bond_diagnostic(_lib(), freq_tol=1.0)
    pairs_set = {(int(r["mode_i"]), int(r["mode_j"])) for _, r in pairs.iterrows()}
    assert pairs_set == {(13, 14), (23, 24)}  # mode 19 pairs with nothing

    for i, j in ((13, 14), (23, 24)):
        row = pairs[(pairs["mode_i"] == i) & (pairs["mode_j"] == j)].iloc[0]
        assert row["complementary"]
        assert row["cc_pattern_correlation"] < -0.9, (i, j, row["cc_pattern_correlation"])
        assert row["delta_freq"] < 0.05


def test_benzene_mixed_bond_diagnostic_raises_if_no_mixed_modes():
    lib_df = _lib().copy()
    mask = lib_df["molecule"] == MOLECULE
    lib_df.loc[mask, "predicted_label"] = lib_df.loc[mask, "predicted_label"].replace(
        "MIXED_STRETCH_BEND", "BENDING")
    try:
        benzene_mixed_bond_diagnostic(lib_df)
        assert False, "expected ValueError when no MIXED_STRETCH_BEND modes exist"
    except ValueError:
        pass


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
