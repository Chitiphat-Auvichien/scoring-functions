"""Benzene real-normal-modes-vs-literature-reference validation -- the
manuscript's PRIMARY classification-vs-reference result (benzene EMIT is a
separate, deliberately-excluded stress test; see src/flag_validation.py's
docstring for why its "ground truth" would be circular).

Ground truth: the literature/group-theory vibrational assignment for each of
benzene's 30 real normal modes predates this code (e.g. Wilson's classic
numbering) and is ingested as ``ref_label`` for molecule C6H6 in
``library_scores.csv``. Comparing the classifier's own ``predicted_label``
against that independent label is a genuine, non-circular accuracy check.

Every function here reuses already-computed ``library_scores.csv`` columns
-- no re-parsing/re-running of Step 1-3 is ever done here. The three-way-only
SB diagnostics below additionally cheaply RE-DERIVE ``predicted_label`` for
every internal row from its own ``V_Stretch`` (``classifier.
rescheme_internal_label``, an explicit ``scheme="threeway"`` default),
since ``library_scores.csv``'s own ``predicted_label`` column is built under
the global default scheme (binary, as of the 2026-08 switch) and
MIXED_STRETCH_BEND never occurs there -- see
``benzene_normal_reference_detail``'s ``scheme``/``thresholds`` parameters:
- ``benzene_normal_reference_detail``/``_summary``: per-mode and per-category
  recall against ``ref_label``, including an explicit "crossed_opposite"
  flag and literature "SB" (mixed) handling for modes 21/22. ``scheme=None``
  by default -- tracks the cached CSV's own scheme (binary).
- ``benzene_internal_confusion_matrix``: 3x3 (bend/stretch/SB x
  bend/stretch/mixed) confusion table for the 30 internal modes.
  ``scheme="threeway"`` by default (V_Stretch rescheme).
- ``benzene_mixed_bond_diagnostic``: per-bond C-C vs. C-H breakdown for
  MIXED_STRETCH_BEND-predicted modes, plus a near-degenerate-pair
  correlation check. ``scheme="threeway"`` by default (V_Stretch rescheme).
- ``benzene_sb_vs_stretch_bond_diagnostic``: contrasts the SB "bending blind
  spot" modes (21/22) against the stretch modes the classifier calls mixed
  (23/24), with the same per-bond evidence. ``scheme="threeway"`` by default
  (V_Stretch rescheme).
- ``benzene_worked_examples``: identifies ring-breathing and representative
  C-H-stretch modes for the manuscript's worked-example gallery. Tracks the
  cached CSV's own scheme (binary) like ``benzene_normal_reference_detail``.

Deliberately benzene-specific (single ``MOLECULE`` constant, hardcoded C-C
ring-bond ordering) -- not generalized to arbitrary molecules, by design.
"""
import os
from itertools import combinations

import numpy as np
import pandas as pd

from src.classifier import classification_bucket, Thresholds, rescheme_internal_label
from src.library_ingest import load_library_scores
from src.scoring import parse_bond_string

MOLECULE = "C6H6"

# Literature ref_label -> predicted bucket that counts as a match. "SB"
# (literature mixed) has no identically-named predicted bucket -- the
# engine's MIXED_STRETCH_BEND collapses to bucket "mixed" -- so it needs an
# explicit override; everything else falls back to ref_label unchanged.
_REF_TO_EXPECTED_PRED = {"SB": "mixed"}


def _expected_pred_bucket(ref_label):
    return _REF_TO_EXPECTED_PRED.get(ref_label, ref_label)

# Ring bonds in cyclic order, matching data/gjf/benzene.com's connectivity.
_CC_RING_BONDS = ("C1-C2", "C2-C3", "C3-C4", "C4-C5", "C5-C6", "C1-C6")
# C-H bonds, one per ring carbon (Ci-H(i+6), same connectivity).
_CH_BONDS = ("C1-H7", "C2-H8", "C3-H9", "C4-H10", "C5-H11", "C6-H12")


# --------------------------------------------------------------------------
# Task A: primary classification-vs-reference validation
# --------------------------------------------------------------------------

def benzene_normal_reference_detail(lib_df=None, data_dir="data", scheme=None, thresholds=None):
    """One row per benzene normal mode with a literature ``ref_label`` (6
    external T/R + 30 internal stretch/bend/SB), comparing
    ``predicted_label`` against ``ref_label`` (via
    ``_expected_pred_bucket`` so literature "SB" counts as correct iff
    predicted bucket is "mixed").

    Columns: mode_index, kind, freq, ref_label, predicted_label,
    predicted_bucket, correct, crossed_opposite, migrated_to_mixed.

    `scheme` (None by default): None reads ``predicted_label`` exactly as
    already computed in ``library_scores.csv`` -- tracks whichever scheme
    ``library_ingest`` last classified the roster under (binary by default,
    the paper-standard scheme as of the 2026-08 switch). An explicit
    "threeway"/"binary" cheaply RE-DERIVES ``predicted_label`` for every
    internal row from its own ``V_Stretch`` (Step 1's score -- always
    scheme-independent) via ``classifier.rescheme_internal_label``, with NO
    re-parse/re-run of Step 1-3 -- needed by the threeway-only SB
    diagnostics below (``benzene_internal_confusion_matrix``,
    ``benzene_mixed_bond_diagnostic``, ``benzene_sb_vs_stretch_bond_diagnostic``),
    since MIXED_STRETCH_BEND never occurs under the binary scheme and those
    functions would otherwise degenerate to an empty "mixed" bucket once
    library_scores.csv's canonical column is binary-only. A Step-2
    external-slot row (e.g. an internal reference mode that won a T/R slot)
    is left untouched, since scheme never touches Step 2/3 -- see
    ``rescheme_internal_label``'s own docstring.

    Raises ValueError if benzene has no ref_label rows (ingest not run) or
    if any ref_label row lacks a predicted_label (incomplete geometry merge).
    """
    df = load_library_scores(data_dir, lib_df)
    b = df[(df["molecule"] == MOLECULE) & df["ref_label"].notna()].copy()
    if b.empty:
        raise ValueError(
            f"No {MOLECULE} rows with a ref_label in library_scores.csv -- "
            "has src.library_ingest.build_library_scores() been run?")

    if scheme is not None:
        th = thresholds or Thresholds.calibrated()
        b["predicted_label"] = b.apply(
            lambda row: rescheme_internal_label(row["predicted_label"], row["V_Stretch"], th, scheme),
            axis=1)

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
    recall + crossing counts, derived from `benzene_normal_reference_detail`
    (never recomputed independently). `n_migrated_to_mixed` is always 0 for
    "SB" by construction: predicted "mixed" already counts as correct there.
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
# Genuine 3-class internal confusion matrix, enabled by the literature "SB"
# label. Mirrors src.calibrate.confusion_matrix_stats's crosstab pattern but
# is not a copy -- that function's external-row/ideal-tier logic doesn't apply.
# --------------------------------------------------------------------------

INTERNAL_REF_CATEGORIES = ("bend", "stretch", "SB")
INTERNAL_PRED_CATEGORIES = ("bend", "stretch", "mixed")


def benzene_internal_confusion_matrix(lib_df=None, data_dir="data", scheme="threeway",
                                       thresholds=None):
    """3x3 reference (bend/stretch/SB) x predicted-bucket (bend/stretch/
    mixed) confusion table + per-category recall for benzene's 30 internal
    normal modes. Built entirely from `benzene_normal_reference_detail`'s
    own columns.

    `scheme` defaults to "threeway" (unlike `benzene_normal_reference_detail`
    itself, whose own default None just tracks the cached CSV): this table's
    entire point is a genuine 3-class (bend/stretch/SB) breakdown, which is
    only meaningful under the three-way scheme -- under "binary" the "mixed"
    predicted column would always be empty. Re-derives benzene's predicted
    labels under `scheme` via `benzene_normal_reference_detail`'s cheap
    V_Stretch-based rescheme path (`thresholds` defaulting to
    `Thresholds.calibrated()`), independent of whatever scheme the cached
    library_scores.csv (`lib_df`) happens to be built under.

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
    detail = benzene_normal_reference_detail(lib_df, data_dir, scheme=scheme, thresholds=thresholds)
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


def run_benzene_internal_confusion(lib_df=None, data_dir="data", write=True,
                                    scheme="threeway", thresholds=None):
    """Headless entry point (Task A'). Returns
    (confusion_table, per_category, (path_table, path_summary) or None)."""
    confusion_table, per_category = benzene_internal_confusion_matrix(
        lib_df, data_dir, scheme=scheme, thresholds=thresholds)
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

def _is_cc_bond(label):
    a, b = label.split("-")
    return a[0] == "C" and b[0] == "C"


def _bond_row_stats(r):
    """Shared per-mode C-C/C-H bond-total breakdown (mode_index, freq,
    V_Stretch, cc_total, ch_total, cc_fraction_of_V, one 's_AB[<bond>]' per
    C-C ring bond) from one `library_scores.csv` row -- shared by both bond
    diagnostics below so they can't silently drift apart."""
    bonds = parse_bond_string(r["s_AB"])
    # s_AB is signed (positive = stretching, negative = compressing); these
    # roll-up totals want magnitude sums (no cancellation between bonds), so
    # cc_fraction_of_V keeps its previously-validated meaning.
    cc_total = sum(abs(v) for k, v in bonds.items() if _is_cc_bond(k))
    ch_total = sum(abs(v) for k, v in bonds.items() if not _is_cc_bond(k))
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


def benzene_mixed_bond_diagnostic(lib_df=None, data_dir="data", freq_tol=1.0,
                                   scheme="threeway", thresholds=None):
    """Per-bond C-C vs. C-H breakdown for benzene's MIXED_STRETCH_BEND-
    predicted normal modes, plus a computed near-degenerate-pair check.

    The set of "mixed" modes is FOUND from
    `benzene_normal_reference_detail` (predicted_bucket == "mixed"), not
    hardcoded, so this stays correct across recalibration.

    `scheme` defaults to "threeway" (see `benzene_internal_confusion_matrix`
    for why): MIXED_STRETCH_BEND never occurs under "binary", so this
    diagnostic would always raise ValueError there -- it is inherently a
    three-way-scheme tool. Re-derives benzene's predicted labels under
    `scheme` (`thresholds` defaulting to `Thresholds.calibrated()`).

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
    detail = benzene_normal_reference_detail(lib_df, data_dir, scheme=scheme, thresholds=thresholds)
    mixed_modes = detail.loc[detail["predicted_bucket"] == "mixed", "mode_index"].tolist()
    if not mixed_modes:
        raise ValueError(
            f"No MIXED_STRETCH_BEND {MOLECULE} normal modes found -- nothing "
            "to diagnose (has calibration or the reference labels changed?).")

    df = load_library_scores(data_dir, lib_df)
    mixed_modes_str = {str(m) for m in mixed_modes}
    lib_rows = df[(df["molecule"] == MOLECULE) &
                  df["mode_index"].astype(str).isin(mixed_modes_str)]

    rows = []
    bond_vectors = {}
    for _, r in lib_rows.iterrows():
        row = _bond_row_stats(r)
        bonds = parse_bond_string(r["s_AB"])
        # s_AB is signed; this correlation measures magnitude-pattern
        # complementarity between near-degenerate partners (which C-C bonds
        # are strongly vs. weakly perturbed), not phase/sign agreement, so
        # use magnitudes here too -- consistent with cc_total/ch_total above
        # and unchanged from before s_AB became signed.
        bond_vectors[row["mode_index"]] = np.array([abs(bonds.get(b, 0.0)) for b in _CC_RING_BONDS])
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


def run_benzene_bond_diagnostic(lib_df=None, data_dir="data", freq_tol=1.0, write=True,
                                 scheme="threeway", thresholds=None):
    """Headless entry point (Task B). Returns
    (bond_detail_df, pairs_df, (path_bonds, path_pairs) or None).
    """
    bond_detail_df, pairs_df = benzene_mixed_bond_diagnostic(
        lib_df, data_dir, freq_tol, scheme=scheme, thresholds=thresholds)
    paths = None
    if write:
        p_bonds = os.path.join(data_dir, "results", "benzene_mixed_bond_diagnostic.csv")
        p_pairs = os.path.join(data_dir, "results", "benzene_mixed_degenerate_pairs.csv")
        bond_detail_df.to_csv(p_bonds, index=False)
        pairs_df.to_csv(p_pairs, index=False)
        paths = (p_bonds, p_pairs)
    return bond_detail_df, pairs_df, paths


# --------------------------------------------------------------------------
# Per-bond contrast between benzene's two "miss" mechanisms exposed by the
# literature SB label: modes 21/22 (ref 'SB', predicted a clean bend --
# V_Stretch 0.09072/0.08119, both <= tau_B=0.17327, the purity gate's
# bending blind spot) vs. modes 23/24 (ref 'stretch', predicted mixed).
# Keyed off ref_label, not predicted_bucket: 21/22 are NOT predicted-mixed,
# so finding them via the predicted bucket would find nothing.
# --------------------------------------------------------------------------

def benzene_sb_vs_stretch_bond_diagnostic(lib_df=None, data_dir="data",
                                           scheme="threeway", thresholds=None):
    """Per-bond C-C vs. C-H breakdown contrasting benzene's literature-'SB'
    modes (21/22, bending blind spot) against the literature-'stretch' modes
    the classifier calls mixed (23/24), with the same per-bond evidence
    `benzene_mixed_bond_diagnostic` uses for the other predicted-mixed modes.

    `scheme` defaults to "threeway" (see `benzene_internal_confusion_matrix`
    for why): the "recovered_mixed"/"overflagged_mixed" cases this
    diagnostic derives are three-way-scheme concepts (the "mixed" bucket
    never occurs under "binary"). Re-derives benzene's predicted labels
    under `scheme` (`thresholds` defaulting to `Thresholds.calibrated()`).

    Returns one row per mode (21, 22, 23, 24): mode_index, freq, V_Stretch,
    ref_label, predicted_label, predicted_bucket, case
    ('overflagged_mixed' for 23/24; for the SB modes the case is derived from
    what was actually predicted -- 'blind_spot_bend' if the classifier called
    them clean bending, 'recovered_mixed' if it called them mixed), cc_total,
    ch_total, cc_fraction_of_V, plus one 's_AB[<bond>]' column per C-C bond.

    Raises ValueError if benzene has no ref_label=='SB' modes (literature
    relabeling not ingested).
    """
    detail = benzene_normal_reference_detail(lib_df, data_dir, scheme=scheme, thresholds=thresholds)
    sb_modes = detail.loc[detail["ref_label"] == "SB", "mode_index"].astype(int).tolist()
    if not sb_modes:
        raise ValueError(
            f"No {MOLECULE} normal modes with ref_label=='SB' found -- "
            "nothing to contrast (has the literature relabeling been "
            "ingested? see src/csv_label_ingest.py's module docstring).")
    # Contrast set: the near-degenerate PAIR among literature-stretch/
    # predicted-mixed modes (excludes lone mode 19, the dedicated worked
    # SB example elsewhere) -- derived by frequency-sharing, not hardcoded.
    candidates = detail[(detail["kind"] == "internal") &
                         (detail["ref_label"] == "stretch") &
                         (detail["predicted_bucket"] == "mixed")].copy()
    freq_group_sizes = candidates.groupby("freq")["mode_index"].transform("count")
    contrast_modes = candidates.loc[freq_group_sizes > 1, "mode_index"].astype(int).tolist()

    # The SB modes' case is DERIVED from what the classifier actually did, not
    # asserted from ref_label. Under the original unweighted V-score both were
    # called clean BENDING -- the tau_B "bending blind spot" this diagnostic
    # was written to expose. Under reduced-mass weighting they come out MIXED,
    # matching the literature, so hardcoding 'blind_spot_bend' would have this
    # table reporting a failure that is no longer happening.
    # Internal rows only -- external rows carry slot names ("Tx") in
    # mode_index, not integers. A dict rather than a Series, because integer
    # .get() on a Series is a positional lookup, not a label one.
    _internal = detail[detail["kind"] == "internal"]
    sb_bucket = dict(zip(_internal["mode_index"].astype(int),
                         _internal["predicted_bucket"]))
    case_by_mode = {
        m: ("blind_spot_bend" if sb_bucket.get(m) == "bend" else
            "recovered_mixed" if sb_bucket.get(m) == "mixed" else
            f"sb_predicted_{sb_bucket.get(m)}")
        for m in sb_modes
    }
    case_by_mode.update({m: "overflagged_mixed" for m in contrast_modes})
    internal_detail = detail[detail["kind"] == "internal"].copy()
    internal_detail["mode_index"] = internal_detail["mode_index"].astype(int)
    detail_by_mode = internal_detail.set_index("mode_index")

    df = load_library_scores(data_dir, lib_df)
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


def run_benzene_sb_vs_stretch_bond_diagnostic(lib_df=None, data_dir="data", write=True,
                                               scheme="threeway", thresholds=None):
    """Headless entry point (Task B'). Returns (df, path or None)."""
    result = benzene_sb_vs_stretch_bond_diagnostic(lib_df, data_dir, scheme=scheme, thresholds=thresholds)
    path = None
    if write:
        path = os.path.join(data_dir, "results", "benzene_sb_vs_stretch_bond_diagnostic.csv")
        result.to_csv(path, index=False)
    return result, path


# --------------------------------------------------------------------------
# Worked-example gallery mode identification: ring-breathing + representative
# C-H stretch (the SB example, mode 19, is already fixed elsewhere and not
# re-derived here). Criterion: among benzene's 7 STRETCHING-labeled modes
# there is a clear ~2200 cm-1 frequency gap (one mode far below 3000 cm-1,
# the rest at/above ~3180 cm-1) -- the low one (V_Stretch ~1.000) is
# ring-breathing (matches the literature ~992 cm-1 assignment), the highest
# is the C-H stretch. Per-bond C-C/C-H totals are reported as supporting
# evidence only, not the primary test.

def benzene_worked_examples(lib_df=None, data_dir="data", freq_tol=1.0):
    """Identify benzene's ring-breathing and representative C-H-stretch
    normal modes for the manuscript's worked-example gallery (see section
    comment above for the criterion).

    Returns a 2-row DataFrame, one row each for role "ring_breathing" and
    "ch_stretch": mode_index, freq, V_Stretch, role, cc_total, ch_total,
    cc_fraction_of_V, cc_min, cc_max, cc_cv, ch_min, ch_max, ch_cv,
    near_degenerate_partner (mode_index of another `S`-labeled mode within
    `freq_tol` cm-1, or None), partner_freq_diff (NaN if none).

    Raises ValueError if fewer than 2 STRETCHING-labeled modes exist.
    """
    df = load_library_scores(data_dir, lib_df)
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
        bonds = parse_bond_string(row["s_AB"])
        # s_AB is signed; these totals want magnitude sums (see
        # _bond_row_stats above for the same rationale).
        cc_total = sum(abs(v) for k, v in bonds.items() if _is_cc_bond(k))
        ch_total = sum(abs(v) for k, v in bonds.items() if not _is_cc_bond(k))
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
