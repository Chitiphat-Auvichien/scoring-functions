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

Two deliverables (module functions), each reusing already-computed
``library_scores.csv`` columns -- no new scores are computed here, no
classifier/thresholds are re-run:

1. ``benzene_normal_reference_detail`` / ``_summary`` / ``run_benzene_normal_
   validation`` (Task A): per-mode + per-category (translation/rotation/
   stretch/bend) recall against ``ref_label``, plus an explicit
   "crossed_opposite" flag (predicted the OPPOSITE clean category, e.g. a
   literature bend predicted STRETCHING) so the "zero crossings" claim is a
   computed count, not an eyeballed one.

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

# Ring bonds in cyclic order, matching data/gjf/benzene.com's connectivity
# (atoms 1-6 are the ring carbons; 1-2,2-3,3-4,4-5,5-6,6-1 are the C-C bonds,
# 1-7..6-12 the C-H bonds). Excel's own atom labels (e.g. "C1", "H7") are
# already embedded in library_scores.csv's s_AB string (src/excel_ingest.py's
# _bond_string()), so no separate atom-type lookup is needed here.
_CC_RING_BONDS = ("C1-C2", "C2-C3", "C3-C4", "C4-C5", "C5-C6", "C1-C6")


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
    stretch/bend rows), comparing the classifier's ``predicted_label``
    against ``ref_label``.

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
    b["correct"] = b["predicted_bucket"] == b["ref_label"]
    b["crossed_opposite"] = (
        ((b["ref_label"] == "bend") & (b["predicted_bucket"] == "stretch")) |
        ((b["ref_label"] == "stretch") & (b["predicted_bucket"] == "bend"))
    )
    b["migrated_to_mixed"] = (~b["correct"]) & (b["predicted_bucket"] == "mixed")

    cols = ["mode_index", "kind", "freq", "ref_label", "predicted_label",
            "predicted_bucket", "correct", "crossed_opposite", "migrated_to_mixed"]
    return b[cols].reset_index(drop=True)


def benzene_normal_reference_summary(detail_df):
    """Per-category (translation/rotation/stretch/bend) n / n_correct /
    recall + crossing counts, derived from `benzene_normal_reference_detail`'s
    output (never recomputed independently of it, so the two are always
    self-consistent).
    """
    rows = []
    for ref in ("translation", "rotation", "stretch", "bend"):
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
        mode_index = int(r["mode_index"])
        bonds = _parse_bond_string(r["s_AB"])
        cc_total = sum(v for k, v in bonds.items() if _is_cc_bond(k))
        ch_total = sum(v for k, v in bonds.items() if not _is_cc_bond(k))
        cc_vec = np.array([bonds.get(b, 0.0) for b in _CC_RING_BONDS])
        bond_vectors[mode_index] = cc_vec
        row = {
            "mode_index": mode_index,
            "freq": float(r["freq"]),
            "V_Stretch": float(r["V_Stretch"]),
            "cc_total": cc_total,
            "ch_total": ch_total,
            "cc_fraction_of_V": cc_total / r["V_Stretch"] if r["V_Stretch"] else float("nan"),
        }
        for bond in _CC_RING_BONDS:
            row[f"s_AB[{bond}]"] = bonds.get(bond, 0.0)
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


if __name__ == "__main__":
    d, s, p = run_benzene_normal_validation()
    print(s.to_string(index=False))
    bd, pr, p2 = run_benzene_bond_diagnostic()
    print(bd.to_string(index=False))
    print(pr.to_string(index=False))
