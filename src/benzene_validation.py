"""Benzene NORMAL-modes-vs-reference validation -- the manuscript's PRIMARY
classification-vs-reference result.

Per the mandated Results & Discussion structure in
``JCC/Scoring_Manuscript_Plan_2026-07-01.pdf`` (one directory above this
repo), "Benzene normal modes -> low-frequency stretching" is the paper's
PRIMARY classification-vs-reference validation -- not benzene EMIT, which is
positioned only as an extreme/rare edge-case stress test (see
``src/flag_validation.py``'s systematic 36-mode EMIT confusion matrix, built
and tested but EXCLUDED from the manuscript by author decision: its "ground
truth" is a threshold cut on continuous, genuinely-mixed EMIT projection
fractions -- circular reasoning for an accuracy claim. That decision is
recorded, not re-litigated, here or anywhere in this module.)

Why benzene's real NORMAL modes are a non-circular ground truth (unlike
EMIT): the literature/group-theory vibrational assignment for each of
benzene's 30 real normal modes is external to this framework entirely (it
predates this code, e.g. Wilson's classic benzene mode numbering) --
ingested as ``ref_label`` in ``data/results/library_scores.csv`` for molecule
``C6H6`` by ``src/excel_ingest.py``. Comparing the classifier's own
``predicted_label`` against that independent label is therefore a genuine
accuracy check, not a threshold circularity.

Module functions, each reusing already-computed ``library_scores.csv``
columns -- no new scores are computed here, no classifier/thresholds are
re-run:

1. ``benzene_normal_reference_detail`` / ``_summary`` / ``run_benzene_normal_
   validation`` (Task A): per-mode + per-category (translation/rotation/
   stretch/bend/SB) recall against ``ref_label``, plus an explicit
   "crossed_opposite" flag (predicted the OPPOSITE clean category, e.g. a
   literature bend predicted STRETCHING) so the "zero crossings" claim is a
   computed count, not an eyeballed one. "SB" (2026-07-05: a genuine
   literature-sourced mixed reference class for benzene modes 21/22, see
   ``src/csv_label_ingest.py``) counts as correctly classified iff the
   predicted bucket is "mixed" (`_expected_pred_bucket`), since the engine
   has no predicted bucket literally named "SB".

1'. ``benzene_internal_confusion_matrix`` / ``run_benzene_internal_confusion``
    (Task A', 2026-07-05): the genuine 3x3 (ref bend/stretch/SB x predicted
    bend/stretch/mixed) confusion table for benzene's 30 internal modes,
    enabled by the literature "SB" class -- mirrors
    ``src/calibrate.py::confusion_matrix_stats``'s ``pd.crosstab`` pattern,
    scoped to benzene only (that function's external-row/ideal-tier
    assumptions don't apply here).

2. ``benzene_mixed_bond_diagnostic`` / ``run_benzene_bond_diagnostic``
   (Task B): for whichever benzene normal modes land in the MIXED_STRETCH_
   BEND bucket (found from (1)'s own output, not hardcoded mode numbers), a
   per-bond C-C vs. C-H breakdown (parsed from the semicolon-joined ``s_AB``
   column) plus a computed near-degenerate-pair check: any two mixed modes
   within ``freq_tol`` cm-1 of each other get their 6 ring-bond ``s_AB``
   vectors correlated (Pearson r) -- a strong NEGATIVE correlation is the
   computed signature of the "complementary alternating pattern" a D6h
   doubly-degenerate (E-type) mode pair is expected to show, replacing an
   eyeballed claim with a number.

2'. ``benzene_sb_vs_stretch_bond_diagnostic`` / ``run_benzene_sb_vs_stretch_
    bond_diagnostic`` (Task B', 2026-07-05): a sibling to Task B, NOT an
    extension of it, keyed off ``ref_label`` rather than
    ``predicted_bucket`` -- contrasts modes 21/22 (literature SB, predicted
    a clean BENDING: the tau_B two-gate purity test's "bending blind spot")
    against the near-degenerate pair 23/24 (literature stretch, predicted
    MIXED_STRETCH_BEND) with the same per-bond evidence, so the manuscript
    can explain both miss mechanisms side by side.

This module is deliberately benzene-specific (single ``MOLECULE`` constant,
hardcoded C-C ring-bond ordering) -- generalizing the bond diagnostic to
arbitrary molecules is explicitly out of scope (task instruction); benzene is
the manuscript's specific worked example.
"""
import os
from itertools import combinations

import numpy as np
import pandas as pd

from src.classifier import classification_bucket

MOLECULE = "C6H6"

# Literature ref_label -> the predicted BUCKET that counts as a correct
# match. Needed because the literal literature "SB" class (2026-07-05
# relabeling of modes 21/22 -- a genuine literature-sourced mixed reference
# label, not just the classifier's own predicted bucket) has no identically-
# named predicted bucket: `classification_bucket()` collapses the engine's
# own MIXED_STRETCH_BEND label to bucket "mixed", not "SB". translation/
# rotation map to themselves (classification_bucket() already returns those
# exact strings for external slots), so this dict only needs an explicit
# override for "SB"; everything else falls back to ref_label unchanged via
# `_expected_pred_bucket()`.
_REF_TO_EXPECTED_PRED = {"SB": "mixed"}


def _expected_pred_bucket(ref_label):
    return _REF_TO_EXPECTED_PRED.get(ref_label, ref_label)

# Ring bonds in cyclic order, matching data/gjf/benzene.com's connectivity
# (atoms 1-6 are the ring carbons; 1-2,2-3,3-4,4-5,5-6,6-1 are the C-C bonds,
# 1-7..6-12 the C-H bonds). Excel's own atom labels (e.g. "C1", "H7") are
# already embedded in library_scores.csv's s_AB string (src/excel_ingest.py's
# _bond_string()), so no separate atom-type lookup is needed here.
_CC_RING_BONDS = ("C1-C2", "C2-C3", "C3-C4", "C4-C5", "C5-C6", "C1-C6")
# C-H bonds, one per ring carbon (Ci-H(i+6), matching the same connectivity).
_CH_BONDS = ("C1-H7", "C2-H8", "C3-H9", "C4-H10", "C5-H11", "C6-H12")


def _load_lib(lib_df, data_dir):
    if lib_df is not None:
        return lib_df
    return pd.read_csv(os.path.join(data_dir, "results", "library_scores.csv"))


# --------------------------------------------------------------------------
# Task A: primary classification-vs-reference validation
# --------------------------------------------------------------------------

def benzene_normal_reference_detail(lib_df=None, data_dir="data"):
    """One row per benzene (C6H6) normal mode carrying a literature
    ``ref_label`` (the 6 external T/R rows + the 30 internal
    stretch/bend/SB rows -- ``"SB"`` is a genuine literature-sourced mixed
    class for 2 modes, 2026-07-05 relabeling; see
    ``src/csv_label_ingest.py``), comparing the classifier's
    ``predicted_label`` against ``ref_label``.

    ``correct`` uses ``_expected_pred_bucket(ref_label)`` rather than
    ``ref_label`` directly, so a literature ``"SB"`` row counts as correct
    iff the predicted bucket is ``"mixed"`` (the framework's own structural
    equivalent -- there is no predicted bucket literally named ``"SB"``).

    Columns: mode_index, kind, freq, ref_label, predicted_label,
    predicted_bucket, correct, crossed_opposite, migrated_to_mixed.

    Raises ValueError if benzene has no ref_label rows at all (ingest not
    run) or if any ref_label row lacks a predicted_label (the geometry
    merge did not happen -- e.g. a frequency mismatch gate failure) --
    fail loud rather than silently validate against an incomplete merge.
    """
    df = _load_lib(lib_df, data_dir)
    b = df[(df["molecule"] == MOLECULE) & df["ref_label"].notna()].copy()
    if b.empty:
        raise ValueError(
            f"No {MOLECULE} rows with a ref_label in library_scores.csv -- "
            "has src.excel_ingest.build_library_scores() been run?")
    missing_pred = b["predicted_label"].isna()
    if missing_pred.any():
        bad = b.loc[missing_pred, "mode_index"].tolist()
        raise ValueError(
            f"{len(bad)} {MOLECULE} row(s) have a ref_label but no "
            f"predicted_label (mode_index={bad}) -- "
            "attach_geometry_classification() did not merge these (frequency "
            "mismatch?). Cannot validate against an incomplete merge.")

    b["predicted_bucket"] = b["predicted_label"].map(classification_bucket)
    expected_pred_bucket = b["ref_label"].map(_expected_pred_bucket)
    b["correct"] = b["predicted_bucket"] == expected_pred_bucket
    b["crossed_opposite"] = (
        ((b["ref_label"] == "bend") & (b["predicted_bucket"] == "stretch")) |
        ((b["ref_label"] == "stretch") & (b["predicted_bucket"] == "bend"))
    )
    b["migrated_to_mixed"] = (~b["correct"]) & (b["predicted_bucket"] == "mixed")

    cols = ["mode_index", "kind", "freq", "ref_label", "predicted_label",
            "predicted_bucket", "correct", "crossed_opposite", "migrated_to_mixed"]
    return b[cols].reset_index(drop=True)


def benzene_normal_reference_summary(detail_df):
    """Per-category (translation/rotation/stretch/bend/SB) n / n_correct /
    recall + crossing counts, derived from `benzene_normal_reference_detail`'s
    output (never recomputed independently of it, so the two are always
    self-consistent). "SB" (2026-07-05 literature relabeling, modes 21/22)
    is included alongside the original 4 categories -- its "correct"/
    "n_migrated_to_mixed" already use the SB-aware criterion computed in
    `benzene_normal_reference_detail` (predicted bucket=="mixed" counts as
    correct for ref_label=="SB", so `n_migrated_to_mixed` is always 0 for
    "SB" by construction: being predicted "mixed" IS correct for this
    category, not a migration away from a clean match).
    """
    rows = []
    for ref in ("translation", "rotation", "stretch", "bend", "SB"):
        sub = detail_df[detail_df["ref_label"] == ref]
        n = len(sub)
        n_correct = int(sub["correct"].sum())
        n_crossed = int(sub["crossed_opposite"].sum())
        n_mixed = int(sub["migrated_to_mixed"].sum())
        rows.append({
            "ref_label": ref,
            "n": n,
            "n_correct": n_correct,
            "recall": n_correct / n if n else float("nan"),
            "n_crossed_opposite": n_crossed,
            "n_migrated_to_mixed": n_mixed,
        })
    return pd.DataFrame(rows)


def run_benzene_normal_validation(lib_df=None, data_dir="data", write=True):
    """Headless entry point (Task A). Returns
    (detail_df, summary_df, (path_detail, path_summary) or None).
    """
    detail = benzene_normal_reference_detail(lib_df, data_dir)
    summary = benzene_normal_reference_summary(detail)
    paths = None
    if write:
        p_detail = os.path.join(data_dir, "results", "benzene_normal_reference_detail.csv")
        p_summary = os.path.join(data_dir, "results", "benzene_normal_reference_summary.csv")
        detail.to_csv(p_detail, index=False)
        summary.to_csv(p_summary, index=False)
        paths = (p_detail, p_summary)
    return detail, summary, paths


# --------------------------------------------------------------------------
# Task A' (2026-07-05): genuine 3-class internal confusion matrix, enabled by
# the literature "SB" relabeling of modes 21/22 (see src/csv_label_ingest.py).
# Reuses src.calibrate.confusion_matrix_stats's pd.crosstab pattern, scoped
# to benzene's 30 internal (stretch/bend/SB) modes -- NOT a copy of that
# function, since its translation/rotation-specific assumptions (external
# rows, ideal/non-ideal tiers) don't apply here at all.
# --------------------------------------------------------------------------

INTERNAL_REF_CATEGORIES = ("bend", "stretch", "SB")
INTERNAL_PRED_CATEGORIES = ("bend", "stretch", "mixed")


def benzene_internal_confusion_matrix(lib_df=None, data_dir="data"):
    """3x3 reference (bend/stretch/SB) x predicted-bucket (bend/stretch/
    mixed) confusion table + per-category recall for benzene's 30 internal
    normal modes -- the first genuine 3-class ground truth in this pipeline
    (literature "SB" is an independent literature-sourced mixed label, not
    the classifier's own predicted bucket).

    Built entirely from `benzene_normal_reference_detail`'s own ref_label/
    predicted_bucket/correct columns (never recomputes score/classification
    logic itself), so it is always self-consistent with Task A's per-
    category summary.

    Returns (confusion_table, per_category):
      confusion_table -- 3x3 DataFrame, index=['bend','stretch','SB'],
        columns=['bend','stretch','mixed'], reindexed with fill_value=0 so
        the shape is always exactly 3x3.
      per_category -- one row per ref in ('bend','stretch','SB'): n,
        n_correct (using `_expected_pred_bucket` -- predicted=='mixed'
        counts as correct for ref=='SB'), recall, n_crossed_opposite (bend
        <-> stretch crossings only; always 0 for 'SB', which has no
        "opposite" clean category to cross into).
    """
    detail = benzene_normal_reference_detail(lib_df, data_dir)
    internal = detail[detail["kind"] == "internal"].copy()

    confusion_table = pd.crosstab(internal["ref_label"], internal["predicted_bucket"])
    confusion_table = confusion_table.reindex(
        index=list(INTERNAL_REF_CATEGORIES), columns=list(INTERNAL_PRED_CATEGORIES), fill_value=0)

    rows = []
    for ref in INTERNAL_REF_CATEGORIES:
        sub = internal[internal["ref_label"] == ref]
        n = len(sub)
        n_correct = int(sub["correct"].sum())
        n_crossed = int(sub["crossed_opposite"].sum())
        rows.append({
            "ref_label": ref,
            "n": n,
            "n_correct": n_correct,
            "recall": n_correct / n if n else float("nan"),
            "n_crossed_opposite": n_crossed,
        })
    per_category = pd.DataFrame(rows)
    return confusion_table, per_category


def run_benzene_internal_confusion(lib_df=None, data_dir="data", write=True):
    """Headless entry point (Task A'). Returns
    (confusion_table, per_category, (path_table, path_summary) or None)."""
    confusion_table, per_category = benzene_internal_confusion_matrix(lib_df, data_dir)
    paths = None
    if write:
        p_table = os.path.join(data_dir, "results", "benzene_internal_confusion_matrix.csv")
        p_summary = os.path.join(data_dir, "results", "benzene_internal_confusion_summary.csv")
        confusion_table.to_csv(p_table)
        per_category.to_csv(p_summary, index=False)
        paths = (p_table, p_summary)
    return confusion_table, per_category, paths


# --------------------------------------------------------------------------
# Task B: bond-contribution / degenerate-pair diagnostic for the mixed modes
# --------------------------------------------------------------------------

def _parse_bond_string(s):
    """'"C1-C2:0.0342;C1-C6:0.0539;..."' -> {'C1-C2': 0.0342, ...}."""
    out = {}
    if not isinstance(s, str) or not s:
        return out
    for part in s.split(";"):
        atoms, val = part.split(":")
        out[atoms] = float(val)
    return out


def _is_cc_bond(label):
    a, b = label.split("-")
    return a[0] == "C" and b[0] == "C"


def _bond_row_stats(r):
    """Shared per-mode C-C/C-H bond-total breakdown (mode_index, freq,
    V_Stretch, cc_total, ch_total, cc_fraction_of_V, one 's_AB[<bond>]' per
    C-C ring bond) from one `library_scores.csv` row -- used by both
    `benzene_mixed_bond_diagnostic` (Task B, keyed off predicted_bucket==
    'mixed') and `benzene_sb_vs_stretch_bond_diagnostic` (Task B',
    2026-07-05, keyed off ref_label=='SB'/'stretch'), so the two diagnostics
    can never silently drift apart in how they compute cc_total/ch_total.
    """
    bonds = _parse_bond_string(r["s_AB"])
    cc_total = sum(v for k, v in bonds.items() if _is_cc_bond(k))
    ch_total = sum(v for k, v in bonds.items() if not _is_cc_bond(k))
    row = {
        "mode_index": int(r["mode_index"]),
        "freq": float(r["freq"]),
        "V_Stretch": float(r["V_Stretch"]),
        "cc_total": cc_total,
        "ch_total": ch_total,
        "cc_fraction_of_V": cc_total / r["V_Stretch"] if r["V_Stretch"] else float("nan"),
    }
    for bond in _CC_RING_BONDS:
        row[f"s_AB[{bond}]"] = bonds.get(bond, 0.0)
    return row


def benzene_mixed_bond_diagnostic(lib_df=None, data_dir="data", freq_tol=1.0):
    """Per-bond C-C vs. C-H breakdown for benzene's MIXED_STRETCH_BEND-
    predicted normal modes, plus a computed near-degenerate-pair check.

    The set of "mixed" modes is FOUND from
    `benzene_normal_reference_detail` (predicted_bucket == "mixed"), not
    hardcoded, so this stays correct across recalibration.

    Returns (bond_detail_df, pairs_df):
      bond_detail_df -- one row per mixed mode: mode_index, freq, V_Stretch,
        cc_total, ch_total, cc_fraction_of_V, plus one 's_AB[<bond>]' column
        per C-C ring bond (cyclic order).
      pairs_df -- one row per pair of mixed modes with |freq_i - freq_j| <=
        freq_tol: mode_i, mode_j, freq_i, freq_j, delta_freq,
        cc_pattern_correlation (Pearson r of the two modes' 6-bond C-C
        vectors), complementary (bool, r < -0.5). Empty (0 rows, correct
        columns) if no pair falls within freq_tol of each other.

    Raises ValueError if no MIXED_STRETCH_BEND modes are found at all
    (nothing to diagnose).
    """
    detail = benzene_normal_reference_detail(lib_df, data_dir)
    mixed_modes = detail.loc[detail["predicted_bucket"] == "mixed", "mode_index"].tolist()
    if not mixed_modes:
        raise ValueError(
            f"No MIXED_STRETCH_BEND {MOLECULE} normal modes found -- nothing "
            "to diagnose (has calibration or the reference labels changed?).")

    df = _load_lib(lib_df, data_dir)
    mixed_modes_str = {str(m) for m in mixed_modes}
    lib_rows = df[(df["molecule"] == MOLECULE) &
                  df["mode_index"].astype(str).isin(mixed_modes_str)]

    rows = []
    bond_vectors = {}
    for _, r in lib_rows.iterrows():
        row = _bond_row_stats(r)
        bonds = _parse_bond_string(r["s_AB"])
        bond_vectors[row["mode_index"]] = np.array([bonds.get(b, 0.0) for b in _CC_RING_BONDS])
        rows.append(row)
    bond_detail_df = pd.DataFrame(rows).sort_values("mode_index").reset_index(drop=True)

    freqs = dict(zip(bond_detail_df["mode_index"], bond_detail_df["freq"]))
    pair_rows = []
    for i, j in combinations(sorted(bond_vectors.keys()), 2):
        delta_freq = abs(freqs[i] - freqs[j])
        if delta_freq <= freq_tol:
            corr = float(np.corrcoef(bond_vectors[i], bond_vectors[j])[0, 1])
            pair_rows.append({
                "mode_i": i, "mode_j": j,
                "freq_i": freqs[i], "freq_j": freqs[j],
                "delta_freq": delta_freq,
                "cc_pattern_correlation": corr,
                "complementary": corr < -0.5,
            })
    pairs_df = pd.DataFrame(
        pair_rows,
        columns=["mode_i", "mode_j", "freq_i", "freq_j", "delta_freq",
                 "cc_pattern_correlation", "complementary"],
    )
    return bond_detail_df, pairs_df


def run_benzene_bond_diagnostic(lib_df=None, data_dir="data", freq_tol=1.0, write=True):
    """Headless entry point (Task B). Returns
    (bond_detail_df, pairs_df, (path_bonds, path_pairs) or None).
    """
    bond_detail_df, pairs_df = benzene_mixed_bond_diagnostic(lib_df, data_dir, freq_tol)
    paths = None
    if write:
        p_bonds = os.path.join(data_dir, "results", "benzene_mixed_bond_diagnostic.csv")
        p_pairs = os.path.join(data_dir, "results", "benzene_mixed_degenerate_pairs.csv")
        bond_detail_df.to_csv(p_bonds, index=False)
        pairs_df.to_csv(p_pairs, index=False)
        paths = (p_bonds, p_pairs)
    return bond_detail_df, pairs_df, paths


# --------------------------------------------------------------------------
# Task B' (2026-07-05): per-bond contrast between benzene's two "miss"
# mechanisms exposed by the literature SB relabeling -- modes 21/22
# (ref_label=='SB', but predicted a CLEAN "bend": the two-gate purity test's
# tau_B bending-blind-spot, V_Stretch 0.09072/0.08119 both <= tau_B=0.17327)
# vs. modes 23/24 (ref_label=='stretch', but predicted MIXED_STRETCH_BEND --
# the classifier correctly flags mixed character where the literature calls
# these pure stretches). Keyed off ref_label, NOT off predicted_bucket=='
# mixed' (unlike Task B/benzene_mixed_bond_diagnostic): 21/22 are
# specifically NOT predicted-mixed, so finding them via the predicted bucket
# would find nothing -- that IS the point of this diagnostic.
# --------------------------------------------------------------------------

def benzene_sb_vs_stretch_bond_diagnostic(lib_df=None, data_dir="data"):
    """Per-bond C-C vs. C-H breakdown contrasting benzene's literature-'SB'
    modes (21/22 -- predicted a clean BENDING, the purity-gate blind spot)
    against the literature-'stretch' modes the classifier itself calls
    mixed (23/24 -- predicted MIXED_STRETCH_BEND), so the manuscript can
    explain both miss mechanisms side by side with the same per-bond
    evidence `benzene_mixed_bond_diagnostic` already uses for the other 3
    predicted-mixed modes (13/14/19).

    Returns one row per mode (21, 22, 23, 24): mode_index, freq, V_Stretch,
    ref_label, predicted_label, predicted_bucket, case
    ('blind_spot_bend' for 21/22, 'overflagged_mixed' for 23/24), cc_total,
    ch_total, cc_fraction_of_V, plus one 's_AB[<bond>]' column per C-C ring
    bond (reusing `_bond_row_stats`, shared with Task B).

    Raises ValueError if benzene has no ref_label=='SB' modes at all (the
    literature relabeling has not been ingested -- see
    src/csv_label_ingest.py).
    """
    detail = benzene_normal_reference_detail(lib_df, data_dir)
    sb_modes = detail.loc[detail["ref_label"] == "SB", "mode_index"].astype(int).tolist()
    if not sb_modes:
        raise ValueError(
            f"No {MOLECULE} normal modes with ref_label=='SB' found -- "
            "nothing to contrast (has the literature relabeling been "
            "ingested? see src/csv_label_ingest.py's module docstring).")
    # The literature-stretch/predicted-mixed set is {19, 23, 24} (Task B);
    # the manuscript specifically wants the NEAR-DEGENERATE PAIR (23, 24) as
    # the direct contrast to 21/22 (also a near-degenerate pair), not the
    # lone mode 19 (already the dedicated worked SB example elsewhere, per
    # Task E's docstring) -- derived by keeping only candidates sharing an
    # (exact-to-4dp) frequency with another candidate, not hardcoded.
    candidates = detail[(detail["kind"] == "internal") &
                         (detail["ref_label"] == "stretch") &
                         (detail["predicted_bucket"] == "mixed")].copy()
    freq_group_sizes = candidates.groupby("freq")["mode_index"].transform("count")
    contrast_modes = candidates.loc[freq_group_sizes > 1, "mode_index"].astype(int).tolist()

    case_by_mode = {m: "blind_spot_bend" for m in sb_modes}
    case_by_mode.update({m: "overflagged_mixed" for m in contrast_modes})
    internal_detail = detail[detail["kind"] == "internal"].copy()
    internal_detail["mode_index"] = internal_detail["mode_index"].astype(int)
    detail_by_mode = internal_detail.set_index("mode_index")

    df = _load_lib(lib_df, data_dir)
    target_modes_str = {str(m) for m in sb_modes + contrast_modes}
    lib_rows = df[(df["molecule"] == MOLECULE) &
                  df["mode_index"].astype(str).isin(target_modes_str)]

    rows = []
    for _, r in lib_rows.iterrows():
        row = _bond_row_stats(r)
        mode_index = row["mode_index"]
        d = detail_by_mode.loc[mode_index]
        row["ref_label"] = d["ref_label"]
        row["predicted_label"] = d["predicted_label"]
        row["predicted_bucket"] = d["predicted_bucket"]
        row["case"] = case_by_mode.get(mode_index)
        rows.append(row)

    cols = ["mode_index", "freq", "V_Stretch", "ref_label", "predicted_label",
            "predicted_bucket", "case", "cc_total", "ch_total", "cc_fraction_of_V"] + \
           [f"s_AB[{b}]" for b in _CC_RING_BONDS]
    return pd.DataFrame(rows).sort_values("mode_index").reset_index(drop=True)[cols]


def run_benzene_sb_vs_stretch_bond_diagnostic(lib_df=None, data_dir="data", write=True):
    """Headless entry point (Task B'). Returns (df, path or None)."""
    result = benzene_sb_vs_stretch_bond_diagnostic(lib_df, data_dir)
    path = None
    if write:
        path = os.path.join(data_dir, "results", "benzene_sb_vs_stretch_bond_diagnostic.csv")
        result.to_csv(path, index=False)
    return result, path


# --------------------------------------------------------------------------
# Task E (2026-07-02): worked-example gallery mode identification
# --------------------------------------------------------------------------
#
# Per Scoring_Manuscript_Plan_2026-07-02.pdf step 3a, the "Benzene normal
# modes" section is being reworked into a descriptive worked-example gallery
# needing three NAMED modes: a ring-breathing mode, a representative C-H
# stretch, and an SB (mixed stretch/bend) example. The SB example is already
# in hand (mode 19, 1319.27 cm-1, described in tab:benzenemixed) and is not
# re-derived here. This identifies the other two by computation.
#
# Primary criterion (author-confirmed 2026-07-02, overriding an earlier
# bond-uniformity-first heuristic): among benzene's 7 `S` (STRETCHING)
# -labeled normal modes there is a clear ~2200 cm-1 frequency gap -- one mode
# sits far below 3000 cm-1, the other six sit at/above ~3180 cm-1. The
# low-frequency one, with V_Stretch essentially exactly 1.000, is the
# ring-breathing mode (matches the literature ~992 cm-1 assignment); the
# highest-frequency one is the representative C-H stretch. Per-bond `s_AB`
# C-C vs. C-H totals are reported as SUPPORTING evidence for both picks
# (reusing the same parsing helpers as Task B), not as the primary test.

def benzene_worked_examples(lib_df=None, data_dir="data", freq_tol=1.0):
    """Identify, by computation, benzene's ring-breathing and representative
    C-H-stretch normal modes for the manuscript's worked-example gallery.

    Returns a 2-row DataFrame, one row each for role "ring_breathing" and
    "ch_stretch": mode_index, freq, V_Stretch, role, cc_total, ch_total,
    cc_fraction_of_V, cc_min, cc_max, cc_cv, ch_min, ch_max, ch_cv,
    near_degenerate_partner (mode_index of another `S`-labeled mode within
    `freq_tol` cm-1, or None), partner_freq_diff (NaN if none).

    Method: among C6H6's internal normal modes with predicted_bucket ==
    "stretch", sort by frequency; the LOWEST-frequency one is
    ring_breathing, the HIGHEST-frequency one is ch_stretch. Raises
    ValueError if fewer than 2 STRETCHING-labeled modes exist (need two
    distinct picks) -- fail loud rather than silently returning a
    degenerate/duplicate identification.
    """
    df = _load_lib(lib_df, data_dir)
    internal = df[(df["molecule"] == MOLECULE) & (df["kind"] == "internal")].copy()
    internal["mode_index"] = internal["mode_index"].astype(int)
    internal["predicted_bucket"] = internal["predicted_label"].map(classification_bucket)

    stretch = internal[internal["predicted_bucket"] == "stretch"].sort_values("freq")
    if len(stretch) < 2:
        raise ValueError(
            f"Only {len(stretch)} STRETCHING-labeled {MOLECULE} normal mode(s) "
            "found -- cannot identify distinct ring-breathing and C-H-stretch "
            "worked examples (has calibration or the reference data changed?).")

    ring_row = stretch.iloc[0]
    ch_row = stretch.iloc[-1]

    def _bond_stats(row):
        bonds = _parse_bond_string(row["s_AB"])
        cc_total = sum(v for k, v in bonds.items() if _is_cc_bond(k))
        ch_total = sum(v for k, v in bonds.items() if not _is_cc_bond(k))
        cc_vec = np.array([bonds.get(b, 0.0) for b in _CC_RING_BONDS])
        ch_vec = np.array([bonds.get(b, 0.0) for b in _CH_BONDS])
        cc_cv = float(cc_vec.std() / cc_vec.mean()) if cc_vec.mean() else float("nan")
        ch_cv = float(ch_vec.std() / ch_vec.mean()) if ch_vec.mean() else float("nan")
        return {
            "cc_total": cc_total,
            "ch_total": ch_total,
            "cc_fraction_of_V": cc_total / row["V_Stretch"] if row["V_Stretch"] else float("nan"),
            "cc_min": float(cc_vec.min()),
            "cc_max": float(cc_vec.max()),
            "cc_cv": cc_cv,
            "ch_min": float(ch_vec.min()),
            "ch_max": float(ch_vec.max()),
            "ch_cv": ch_cv,
        }

    # Near-degenerate partner check for the C-H stretch pick (an E1u/E2g-style
    # doubly-degenerate partner would sit within freq_tol cm-1 of ch_row).
    others = stretch[stretch["mode_index"] != ch_row["mode_index"]]
    diffs = (others["freq"] - ch_row["freq"]).abs()
    if len(diffs) and diffs.min() <= freq_tol:
        partner_idx = int(others.loc[diffs.idxmin(), "mode_index"])
        partner_diff = float(diffs.min())
    else:
        partner_idx = None
        partner_diff = float("nan")

    rows = []
    for role, row in (("ring_breathing", ring_row), ("ch_stretch", ch_row)):
        out = {
            "mode_index": int(row["mode_index"]),
            "freq": float(row["freq"]),
            "V_Stretch": float(row["V_Stretch"]),
            "role": role,
        }
        out.update(_bond_stats(row))
        if role == "ch_stretch":
            out["near_degenerate_partner"] = partner_idx
            out["partner_freq_diff"] = partner_diff
        else:
            out["near_degenerate_partner"] = None
            out["partner_freq_diff"] = float("nan")
        rows.append(out)

    cols = ["mode_index", "freq", "V_Stretch", "role", "cc_total", "ch_total",
            "cc_fraction_of_V", "cc_min", "cc_max", "cc_cv", "ch_min", "ch_max",
            "ch_cv", "near_degenerate_partner", "partner_freq_diff"]
    return pd.DataFrame(rows)[cols]


def run_benzene_worked_examples(lib_df=None, data_dir="data", freq_tol=1.0, write=True):
    """Headless entry point (Task E). Returns (df, path or None)."""
    result = benzene_worked_examples(lib_df, data_dir, freq_tol)
    path = None
    if write:
        path = os.path.join(data_dir, "results", "benzene_worked_examples.csv")
        result.to_csv(path, index=False)
    return result, path


if __name__ == "__main__":
    d, s, p = run_benzene_normal_validation()
    print(s.to_string(index=False))
    ct, pc, p_ct = run_benzene_internal_confusion()
    print(ct.to_string())
    print(pc.to_string(index=False))
    bd, pr, p2 = run_benzene_bond_diagnostic()
    print(bd.to_string(index=False))
    print(pr.to_string(index=False))
    sb, p4 = run_benzene_sb_vs_stretch_bond_diagnostic()
    print(sb.to_string(index=False))
    we, p3 = run_benzene_worked_examples()
    print(we.to_string(index=False))
