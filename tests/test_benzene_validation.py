"""Regression tests for src/benzene_validation.py -- the manuscript's PRIMARY
classification-vs-reference validation (benzene's real NORMAL modes vs.
literature/group-theory ref_label; see JCC/Scoring_Manuscript_Plan_2026-07-01.pdf
for the mandated Results & Discussion positioning, and the module's own
docstring for why this is non-circular, unlike the excluded EMIT confusion
matrix in src/flag_validation.py).

Run from ``Github/scoring-functions/``:
    py -m pytest tests/                          (with pytest)
    py tests/test_benzene_validation.py         (standalone; no pytest needed)

Design note (mirrors tests/test_calibrate.py / tests/test_library_ingest.py):
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
    """Pins the exact headline numbers (Task A). History: 6/6 external
    correct throughout. Under the three-way scheme (see
    test_benzene_internal_confusion_matrix_3x3_matches_hand_derived_table
    for that pinned table): 6/10 stretch, 13/18 bend, 2/2 SB recalled, with
    the remainder migrating to the predicted MIXED bucket, zero
    opposite-category (bend<->stretch) crossings.

    **Re-derived 2026-08-19 (binary-classification-scheme default switch)**:
    `benzene_normal_reference_detail()`'s own default (`scheme=None`) now
    tracks `library_scores.csv`'s own scheme, binary as of this switch --
    forcing every internal mode to bend or stretch, never "mixed". Since
    literature ref_label=='SB' (modes 21/22) can only ever be "correct"
    against a predicted "mixed" bucket (`_expected_pred_bucket`), which
    binary never produces, SB recall is now 0/2 by construction (NOT a
    regression -- the binary scheme has no answer for "is this mode
    genuinely mixed?", only "is it more stretch-like or bend-like?"). The
    previously-migrated-to-mixed modes are now forced into a definite S/B
    answer instead: stretch recall rises to 10/10 (perfect -- every
    literature-stretch mode's V_Stretch clears tau_SB=0.5) and bend recall
    is 17/18, with exactly ONE opposite-category crossing this time: mode 14
    (1056.39 cm-1) crosses bend->stretch (its V_Stretch clears tau_SB even
    though its literature label is bend) -- the three-way scheme's genuine
    3-class table (via benzene_internal_confusion_matrix, scheme="threeway"
    default) is unaffected and still shows zero crossings, since this
    binary-forced call is a different question, not a threeway "miss".
    n_migrated_to_mixed is 0 across every category under binary by
    construction (binary never produces a "mixed" predicted bucket at
    all).

    **Re-derived again 2026-08-20 (canonical tau_SB default changed
    0.50->0.42)**: `data/results/thresholds.json`'s top-level `tau_SB` key
    (read by `Thresholds.calibrated()`, which `library_scores.csv`'s own
    `predicted_label` column was regenerated against) moved down from 0.50
    to 0.42, so more borderline internal modes now clear the stretch
    cutoff. Benzene's near-degenerate bend pair at 1056.3901 cm-1 (modes 13
    and 14, V_Stretch 0.451177 and 0.508777) is the exact mechanism: mode
    14 (0.508777) was already above the old 0.50 cutoff and crossed
    bend->stretch before this change; mode 13 (0.451177) sat between the
    two cutoffs (0.42 < 0.451177 < 0.50) and now crosses too. So bend
    n_crossed_opposite rises from 1 to 2 and bend n_correct falls from 17
    to 16 (recall 17/18 -> 16/18 = 0.888889). Nothing else in this table
    moves: stretch stays 10/10 (every literature-stretch mode's V_Stretch
    already cleared even the old, higher 0.50 cutoff, so lowering it
    further changes nothing there), and SB stays 0/2 by the same
    binary-scheme construction argument as the 2026-08-19 paragraph above
    (predicted "mixed" is structurally unreachable under binary,
    independent of where the S/B cutoff itself sits)."""
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
    assert int(by_ref.loc["stretch", "n_correct"]) == 10
    assert by_ref.loc["stretch", "recall"] == 1.0
    assert int(by_ref.loc["stretch", "n_migrated_to_mixed"]) == 0

    assert int(by_ref.loc["bend", "n"]) == 18
    assert int(by_ref.loc["bend", "n_correct"]) == 16
    assert abs(by_ref.loc["bend", "recall"] - 16 / 18) < 1e-9
    assert int(by_ref.loc["bend", "n_migrated_to_mixed"]) == 0
    assert int(by_ref.loc["bend", "n_crossed_opposite"]) == 2

    assert int(by_ref.loc["SB", "n"]) == 2
    assert int(by_ref.loc["SB", "n_correct"]) == 0
    assert by_ref.loc["SB", "recall"] == 0.0
    assert int(by_ref.loc["SB", "n_migrated_to_mixed"]) == 0
    assert int(by_ref.loc["SB", "n_crossed_opposite"]) == 0

    assert int(by_ref.loc["stretch", "n_crossed_opposite"]) == 0


def test_benzene_normal_summary_threeway_matches_prior_pinned_numbers():
    """The three-way scheme's own numbers (explicit scheme="threeway",
    bypassing the now-binary-default library_scores.csv column via the
    cheap V_Stretch rescheme -- see benzene_normal_reference_detail's
    docstring), pinned 2026-08-14 under reduced-mass V-score weighting and
    UNCHANGED by the 2026-08-19 binary-default switch (V_Stretch itself
    never depends on scheme, only which cutoffs are applied to it): 6/10
    stretch, 13/18 bend, 2/2 SB recalled -- see
    test_benzene_internal_confusion_matrix_3x3_matches_hand_derived_table
    for the full 3x3 breakdown these numbers are consistent with. Zero
    opposite-category (bend<->stretch) crossings."""
    detail = benzene_normal_reference_detail(_lib(), scheme="threeway")
    summary = benzene_normal_reference_summary(detail)
    by_ref = summary.set_index("ref_label")

    assert int(by_ref.loc["stretch", "n"]) == 10
    assert int(by_ref.loc["stretch", "n_correct"]) == 6
    assert abs(by_ref.loc["stretch", "recall"] - 0.6) < 1e-9
    assert int(by_ref.loc["stretch", "n_migrated_to_mixed"]) == 4

    assert int(by_ref.loc["bend", "n"]) == 18
    assert int(by_ref.loc["bend", "n_correct"]) == 13
    assert abs(by_ref.loc["bend", "recall"] - 13 / 18) < 1e-9
    assert int(by_ref.loc["bend", "n_migrated_to_mixed"]) == 5

    assert int(by_ref.loc["SB", "n"]) == 2
    assert int(by_ref.loc["SB", "n_correct"]) == 2
    assert by_ref.loc["SB", "recall"] == 1.0
    assert int(by_ref.loc["SB", "n_migrated_to_mixed"]) == 0

    assert int(by_ref.loc["stretch", "n_crossed_opposite"]) == 0
    assert int(by_ref.loc["bend", "n_crossed_opposite"]) == 0
    assert int(by_ref.loc["SB", "n_crossed_opposite"]) == 0


def test_benzene_internal_confusion_matrix_3x3_matches_hand_derived_table():
    """Pins the exact 3x3 (bend/stretch/SB x bend/stretch/mixed) confusion
    table and per-category recall for benzene's 30 internal normal modes,
    author-hand-derived and independently re-verified this session:
      ref bend    (n=18): 13 correct (bend), 5 -> mixed              recall 0.722
      ref stretch (n=10):  6 correct (stretch), 4 -> mixed           recall 0.600
      ref SB      (n=2):   2 correct (mixed)                         recall 1.000
    Zero bend<->stretch crossings (n_crossed_opposite==0 for both).

    **Re-derived 2026-08-14 (reduced-mass V-score weighting)** from
    16/0/2, 0/7/3, 2/0/0. Benzene mixes C-C (mu 6.005) with C-H (mu 0.930),
    so it genuinely re-scores; every one of the six changes is a move INTO
    the mixed column. The SB row is the headline: it goes from 2/0/0 (both
    reference-mixed modes called clean BENDING -- what this file previously
    called the tau_B "bending blind spot") to 0/0/2, i.e. the mixed bucket
    now recovers both. The pure buckets pay for it, 16->13 and 7->6."""
    confusion_table, per_category = benzene_internal_confusion_matrix(_lib())

    assert list(confusion_table.index) == ["bend", "stretch", "SB"]
    assert list(confusion_table.columns) == ["bend", "stretch", "mixed"]
    assert confusion_table.loc["bend"].tolist() == [13, 0, 5]
    assert confusion_table.loc["stretch"].tolist() == [0, 6, 4]
    assert confusion_table.loc["SB"].tolist() == [0, 0, 2]

    by_ref = per_category.set_index("ref_label")
    assert int(by_ref.loc["bend", "n"]) == 18
    assert int(by_ref.loc["bend", "n_correct"]) == 13
    assert abs(by_ref.loc["bend", "recall"] - 13 / 18) < 1e-9
    assert int(by_ref.loc["bend", "n_crossed_opposite"]) == 0

    assert int(by_ref.loc["stretch", "n"]) == 10
    assert int(by_ref.loc["stretch", "n_correct"]) == 6
    assert abs(by_ref.loc["stretch", "recall"] - 0.6) < 1e-9
    assert int(by_ref.loc["stretch", "n_crossed_opposite"]) == 0

    assert int(by_ref.loc["SB", "n"]) == 2
    assert int(by_ref.loc["SB", "n_correct"]) == 2
    assert by_ref.loc["SB", "recall"] == 1.0
    assert int(by_ref.loc["SB", "n_crossed_opposite"]) == 0


def test_benzene_sb_vs_stretch_bond_diagnostic_finds_exactly_21_22_23_24():
    """The SB-vs-stretch contrast diagnostic (Task B') is keyed off
    ref_label, not predicted_bucket, so modes 21/22 are found whatever they
    were predicted as. The contrast pair is specifically 23/24 (the
    near-degenerate literature-stretch/predicted-mixed pair), not mode 19 (a
    lone literature-stretch/predicted-mixed mode, already the dedicated SB
    worked example elsewhere).

    **2026-08-14 (reduced-mass V-score weighting)**: 21/22 are now predicted
    MIXED rather than clean BENDING, so their `case` is 'recovered_mixed'
    instead of 'blind_spot_bend'. That field is now derived from the actual
    prediction rather than hardcoded from ref_label -- otherwise the table
    would keep reporting a blind spot that has stopped happening."""
    result = benzene_sb_vs_stretch_bond_diagnostic(_lib())
    assert sorted(result["mode_index"].tolist()) == [21, 22, 23, 24]
    by_mode = result.set_index("mode_index")
    assert by_mode.loc[21, "case"] == "recovered_mixed"
    assert by_mode.loc[22, "case"] == "recovered_mixed"
    assert by_mode.loc[23, "case"] == "overflagged_mixed"
    assert by_mode.loc[24, "case"] == "overflagged_mixed"
    assert by_mode.loc[21, "ref_label"] == "SB"
    assert by_mode.loc[21, "predicted_bucket"] == "mixed"
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


# Modes the classifier calls MIXED. Under reduced-mass weighting (2026-08-14)
# this grew from 5 to 11: C-C bonds carry mu 6.005 against C-H's 0.930, so ring
# motion that the unweighted score treated as marginal now dominates.
_MIXED_MODES = [13, 14, 16, 17, 18, 19, 21, 22, 23, 24, 25]
# Mode 25 (3182.29 cm-1) is the odd one out: a C-H stretch, not a ring mode.
# It entered the mixed bucket from the other side -- its score FELL (0.982 ->
# 0.894, below tau_S) because the C-C bonds it barely moves now count for more.
_CH_STRETCH_MODE = 25


def test_benzene_mixed_bond_diagnostic_finds_exactly_the_known_modes():
    bond_detail, pairs = benzene_mixed_bond_diagnostic(_lib())
    assert sorted(bond_detail["mode_index"].tolist()) == _MIXED_MODES


def test_benzene_mixed_bond_diagnostic_ch_contributions_are_near_zero():
    """For the RING mixed modes, C-H bond contributions are ~0 -- essentially
    all of V_Stretch sits on the 6 C-C ring bonds, not noise scattered across
    all 12 bonds.

    Mode 25 is deliberately excluded and asserted separately: it is the C-H
    stretch, so C-H carrying everything is the correct answer for it, not a
    violation. Lumping it in would have this test assert that benzene's C-H
    stretch has no C-H character."""
    bond_detail, _ = benzene_mixed_bond_diagnostic(_lib())
    ring = bond_detail[bond_detail["mode_index"] != _CH_STRETCH_MODE]
    assert (ring["ch_total"] < 0.02).all(), ring[["mode_index", "ch_total"]]
    assert (ring["cc_fraction_of_V"] > 0.95).all(), ring[["mode_index", "cc_fraction_of_V"]]

    ch = bond_detail[bond_detail["mode_index"] == _CH_STRETCH_MODE].iloc[0]
    assert ch["cc_total"] == 0.0
    assert ch["ch_total"] > 0.85
    assert ch["cc_fraction_of_V"] == 0.0


def test_benzene_mixed_bond_diagnostic_mode_19_uniform_across_ring():
    """Mode 19 (1319.2678 cm-1) shows a perfectly uniform 6-fold-symmetric
    |C-C| contribution -- computed as a low coefficient of variation across
    the 6 ring bonds, not eyeballed. s_AB is signed (stretch vs. compress);
    mode 19 alternates sign around the ring with equal magnitude, so
    uniformity is a magnitude property, not a raw-signed-value one."""
    bond_detail, _ = benzene_mixed_bond_diagnostic(_lib())
    row = bond_detail[bond_detail["mode_index"] == 19].iloc[0]
    cc_vals = [abs(row[f"s_AB[{b}]"]) for b in
               ("C1-C2", "C2-C3", "C3-C4", "C4-C5", "C5-C6", "C1-C6")]
    cc_vals = pd.Series(cc_vals)
    cv = cc_vals.std() / cc_vals.mean()
    assert cv < 0.02, cv  # essentially flat across all 6 ring bonds


def test_benzene_mixed_bond_diagnostic_detects_near_degenerate_complementary_pairs():
    """The near-degenerate pairs (freq splitting <0.03 cm-1) must be detected
    within a 1 cm-1 tolerance and show a strong NEGATIVE (complementary) C-C
    bond-pattern correlation -- turning the "complementary alternating
    pattern" claim from an eyeballed observation into a computed fact.

    **2026-08-14 (reduced-mass weighting)**: two more degenerate pairs,
    (17,18) and (21,22), joined the mixed bucket, and both show the same
    complementary structure (r = -0.9996 and -0.9970). Modes 16, 19 and 25
    are non-degenerate and pair with nothing."""
    bond_detail, pairs = benzene_mixed_bond_diagnostic(_lib(), freq_tol=1.0)
    pairs_set = {(int(r["mode_i"]), int(r["mode_j"])) for _, r in pairs.iterrows()}
    assert pairs_set == {(13, 14), (17, 18), (21, 22), (23, 24)}

    for i, j in ((13, 14), (17, 18), (21, 22), (23, 24)):
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
    # 2026-08-14 reduced-mass weighting: 0.02 -> 0.05. The mode's actual
    # motion is unchanged; what changed is that its slight C-C participation
    # is now weighted by mu(C-C)=6.005 against mu(C-H)=0.930, so the same
    # motion books ~6.5x the score. cc_total 0.0438, still clearly negligible
    # beside ch_total 0.956. This is the same mechanism that pushes the OTHER
    # C-H stretch, mode 25, below tau_S into the mixed bucket.
    assert by_role.loc["ch_stretch", "cc_total"] < 0.05  # negligible C-C character
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
    """scheme="binary" explicitly: benzene_mixed_bond_diagnostic's own
    default (scheme="threeway") cheaply RE-DERIVES predicted_label from
    V_Stretch (rescheme_internal_label), which would ignore a hand-mutated
    predicted_label column and still find the real threeway "mixed" modes --
    so the no-mixed-modes edge case has to be forced via the scheme itself
    (binary never produces "mixed" at all, by construction) rather than by
    mutating the (now-irrelevant, since it gets rescheme'd away)
    predicted_label column the way this test used to."""
    lib_df = _lib().copy()
    try:
        benzene_mixed_bond_diagnostic(lib_df, scheme="binary")
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
