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
    benzene_mixed_bond_diagnostic, benzene_worked_examples, MOLECULE,
    benzene_internal_confusion_matrix, benzene_sb_vs_stretch_bond_diagnostic,
)
from src.classifier import MIXED_STRETCH_BEND, BENDING               # noqa: E402

LIB_CSV = os.path.join(ROOT, "data", "results", "library_scores.csv")


def _lib():
    return pd.read_csv(LIB_CSV)


def test_benzene_normal_detail_covers_all_36_modes():
    detail = benzene_normal_reference_detail(_lib())
    assert len(detail) == 36  # 6 external + 30 internal, all literature-labeled
    # "SB" (2026-07-05 literature relabeling of modes 21/22 -- a genuine
    # literature-sourced 3rd class, see src/csv_label_ingest.py) joins the
    # original 4 categories.
    assert set(detail["ref_label"].unique()) == {"translation", "rotation", "stretch", "bend", "SB"}
    assert detail["predicted_label"].notna().all()


def test_benzene_normal_summary_matches_ad_hoc_session_numbers():
    """Pins the exact headline numbers (Task A), **updated 2026-07-05** for
    the literature relabeling (mode_index 19/23/24 bend->stretch, 21/22
    bend->SB -- see src/csv_label_ingest.py): 6/6 external correct; 7/10
    literature stretch modes recalled (the 3 newly-stretch modes 19/23/24
    still migrate to the predicted MIXED bucket, unchanged from before);
    16/18 literature bend modes recalled (the 5 original bend "misses"
    minus the 3 that moved to stretch = 2 remaining migrated-to-mixed,
    modes 13/14); 0/2 literature SB modes recalled (modes 21/22 are BOTH
    predicted a clean BENDING, not MIXED -- the tau_B purity-gate "bending
    blind spot"; see benzene_internal_confusion_matrix for the full 3-class
    table). ZERO crossings into the opposite clean category in either
    direction (bend<->stretch), unaffected by the relabeling."""
    detail = benzene_normal_reference_detail(_lib())
    summary = benzene_normal_reference_summary(detail)
    by_ref = summary.set_index("ref_label")

    assert int(by_ref.loc["translation", "n"]) == 3
    assert int(by_ref.loc["translation", "n_correct"]) == 3
    assert by_ref.loc["translation", "recall"] == 1.0

    assert int(by_ref.loc["rotation", "n"]) == 3
    assert int(by_ref.loc["rotation", "n_correct"]) == 3
    assert by_ref.loc["rotation", "recall"] == 1.0

    assert int(by_ref.loc["stretch", "n"]) == 10
    assert int(by_ref.loc["stretch", "n_correct"]) == 7
    assert abs(by_ref.loc["stretch", "recall"] - 0.7) < 1e-9
    assert int(by_ref.loc["stretch", "n_migrated_to_mixed"]) == 3

    assert int(by_ref.loc["bend", "n"]) == 18
    assert int(by_ref.loc["bend", "n_correct"]) == 16
    assert abs(by_ref.loc["bend", "recall"] - 16 / 18) < 1e-9
    assert int(by_ref.loc["bend", "n_migrated_to_mixed"]) == 2

    assert int(by_ref.loc["SB", "n"]) == 2
    assert int(by_ref.loc["SB", "n_correct"]) == 0
    assert by_ref.loc["SB", "recall"] == 0.0
    # Both 21/22 predicted BENDING, not MIXED -- not a "migration to mixed"
    # by this stat's definition (that flag means "incorrect AND predicted
    # mixed"; here it's "incorrect AND predicted bend").
    assert int(by_ref.loc["SB", "n_migrated_to_mixed"]) == 0

    # The zero-crossings claim, computed, not eyeballed.
    assert int(by_ref.loc["stretch", "n_crossed_opposite"]) == 0
    assert int(by_ref.loc["bend", "n_crossed_opposite"]) == 0
    assert int(by_ref.loc["SB", "n_crossed_opposite"]) == 0


def test_benzene_internal_confusion_matrix_3x3_matches_hand_derived_table():
    """Pins the exact 3x3 (bend/stretch/SB x bend/stretch/mixed) confusion
    table and per-category recall for benzene's 30 internal normal modes,
    author-hand-derived and independently re-verified this session:
      ref bend    (n=18): 16 correct (bend), 2 -> mixed (13,14)      recall 0.889
      ref stretch (n=10):  7 correct (stretch), 3 -> mixed (19,23,24) recall 0.700
      ref SB      (n=2):   0 correct, BOTH -> bend (21,22)            recall 0.000
    Zero bend<->stretch crossings (n_crossed_opposite==0 for both)."""
    confusion_table, per_category = benzene_internal_confusion_matrix(_lib())

    assert list(confusion_table.index) == ["bend", "stretch", "SB"]
    assert list(confusion_table.columns) == ["bend", "stretch", "mixed"]
    assert confusion_table.loc["bend"].tolist() == [16, 0, 2]
    assert confusion_table.loc["stretch"].tolist() == [0, 7, 3]
    assert confusion_table.loc["SB"].tolist() == [2, 0, 0]

    by_ref = per_category.set_index("ref_label")
    assert int(by_ref.loc["bend", "n"]) == 18
    assert int(by_ref.loc["bend", "n_correct"]) == 16
    assert abs(by_ref.loc["bend", "recall"] - 16 / 18) < 1e-9
    assert int(by_ref.loc["bend", "n_crossed_opposite"]) == 0

    assert int(by_ref.loc["stretch", "n"]) == 10
    assert int(by_ref.loc["stretch", "n_correct"]) == 7
    assert abs(by_ref.loc["stretch", "recall"] - 0.7) < 1e-9
    assert int(by_ref.loc["stretch", "n_crossed_opposite"]) == 0

    assert int(by_ref.loc["SB", "n"]) == 2
    assert int(by_ref.loc["SB", "n_correct"]) == 0
    assert by_ref.loc["SB", "recall"] == 0.0
    assert int(by_ref.loc["SB", "n_crossed_opposite"]) == 0


def test_benzene_sb_vs_stretch_bond_diagnostic_finds_exactly_21_22_23_24():
    """The SB-vs-stretch contrast diagnostic (Task B') is keyed off
    ref_label, not predicted_bucket -- modes 21/22 (ref SB, predicted clean
    bend) must be found even though they are NOT in the predicted-mixed
    bucket. The contrast pair is specifically 23/24 (the near-degenerate
    literature-stretch/predicted-mixed pair), not mode 19 (a lone
    literature-stretch/predicted-mixed mode, already the dedicated SB
    worked example elsewhere)."""
    result = benzene_sb_vs_stretch_bond_diagnostic(_lib())
    assert sorted(result["mode_index"].tolist()) == [21, 22, 23, 24]
    by_mode = result.set_index("mode_index")
    assert by_mode.loc[21, "case"] == "blind_spot_bend"
    assert by_mode.loc[22, "case"] == "blind_spot_bend"
    assert by_mode.loc[23, "case"] == "overflagged_mixed"
    assert by_mode.loc[24, "case"] == "overflagged_mixed"
    assert by_mode.loc[21, "ref_label"] == "SB"
    assert by_mode.loc[21, "predicted_bucket"] == "bend"
    assert by_mode.loc[23, "ref_label"] == "stretch"
    assert by_mode.loc[23, "predicted_bucket"] == "mixed"


def test_benzene_sb_vs_stretch_bond_diagnostic_raises_if_no_sb_modes():
    lib_df = _lib().copy()
    mask = (lib_df["molecule"] == MOLECULE) & (lib_df["ref_label"] == "SB")
    lib_df.loc[mask, "ref_label"] = "bend"
    try:
        benzene_sb_vs_stretch_bond_diagnostic(lib_df)
        assert False, "expected ValueError when no ref_label=='SB' modes exist"
    except ValueError:
        pass


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


def test_benzene_worked_examples_identifies_ring_breathing_and_ch_stretch():
    """Manuscript claim (Scoring_Manuscript_Plan_2026-07-02.pdf step 3a): pins
    the two newly-identified named worked-example mode indices, protected the
    same way mode 13/14/19/23/24 are pinned above. Ring-breathing = mode 12
    (992.58 cm-1, literature ~992 cm-1, V_Stretch==1.0000, essentially all of
    it on the 6 C-C ring bonds, uniform to <0.3% CV). Representative C-H
    stretch = mode 30 (3223.17 cm-1, the highest-frequency STRETCHING mode,
    V_Stretch==1.0000, essentially all of it on the 6 C-H bonds, non-
    degenerate/isolated in frequency)."""
    df = benzene_worked_examples(_lib())
    by_role = df.set_index("role")

    assert int(by_role.loc["ring_breathing", "mode_index"]) == 12
    assert abs(by_role.loc["ring_breathing", "freq"] - 992.5825) < 0.01
    assert abs(by_role.loc["ring_breathing", "V_Stretch"] - 1.0) < 1e-6
    assert by_role.loc["ring_breathing", "cc_fraction_of_V"] > 0.95
    assert by_role.loc["ring_breathing", "cc_cv"] < 0.02  # uniform across all 6 C-C bonds
    assert by_role.loc["ring_breathing", "ch_total"] < 0.02  # negligible C-H character

    assert int(by_role.loc["ch_stretch", "mode_index"]) == 30
    assert abs(by_role.loc["ch_stretch", "freq"] - 3223.172) < 0.01
    assert abs(by_role.loc["ch_stretch", "V_Stretch"] - 1.0) < 1e-6
    assert by_role.loc["ch_stretch", "ch_total"] > 0.95  # dominant C-H character
    assert by_role.loc["ch_stretch", "cc_total"] < 0.02  # negligible C-C character
    # Well separated in frequency from the ring-breathing pick, and (unlike
    # the 26/27, 28/29 near-degenerate pairs among the other S-labeled modes)
    # not itself part of a near-degenerate pair.
    assert by_role.loc["ch_stretch", "freq"] - by_role.loc["ring_breathing", "freq"] > 2000
    assert by_role.loc["ch_stretch", "near_degenerate_partner"] is None


def test_benzene_worked_examples_raises_if_fewer_than_two_stretch_modes():
    lib_df = _lib().copy()
    mask = (lib_df["molecule"] == MOLECULE) & (lib_df["predicted_label"] == "S")
    idxs = lib_df[mask].index
    # Collapse all but one STRETCHING mode into BENDING, leaving only 1.
    lib_df.loc[idxs[1:], "predicted_label"] = BENDING
    try:
        benzene_worked_examples(lib_df)
        assert False, "expected ValueError with fewer than 2 STRETCHING modes"
    except ValueError:
        pass


def test_benzene_mixed_bond_diagnostic_raises_if_no_mixed_modes():
    lib_df = _lib().copy()
    mask = lib_df["molecule"] == MOLECULE
    lib_df.loc[mask, "predicted_label"] = lib_df.loc[mask, "predicted_label"].replace(
        MIXED_STRETCH_BEND, BENDING)
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
