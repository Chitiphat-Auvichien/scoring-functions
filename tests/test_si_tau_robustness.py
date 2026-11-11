"""Tests for src/si_tau_robustness.py (SI: robustness of tau_SB and the S/B ambiguity range,
JCC revision, referee comments R1.9/R1.10).

Run from ``Github/scoring-functions/``:  py -m pytest tests/test_si_tau_robustness.py

Fast checks with small replicate counts on the tracked data/results/library_scores.csv:
  * the vectorised staircase and tied-run search agree with a direct count / loop;
  * SC set: minimum-error gap (0.407, 0.444], 1/409 SC errors at the adopted 0.43;
  * D1 (molecule bootstrap) keeps 0.43 optimal in every resample;
  * L1: only the SBr4 135 relabelling removes 0.43 from the optimum;
  * the ambiguity range (C1 + score noise 0.125, 95 %) is roughly [0.31, 0.53].
The flip-unit checks need the degenerate-set CSV and are skipped without it.
"""
import os
import shutil
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.chdir(ROOT)

from src import si_tau_robustness as TR  # noqa: E402
from src.figures import AMBIGUITY_RANGE  # noqa: E402

LIB = os.path.join("data", "results", "library_scores.csv")
TAU = 0.43
GI = int(np.argmin(np.abs(TR.GRID - TAU)))
# Degenerate sets: the repo file, or (until it is generated) the R1.4 side-report copy, same columns.
_UNIT_CANDIDATES = (TR.DEFAULT_UNITS_CSV,
                    os.path.join("..", "..", "JCC", "2nd submission", "side_reports",
                                 "R1.4_degenerate_modes", "degenerate_sets.csv"))
UNITS_CSV = next((p for p in _UNIT_CANDIDATES if os.path.exists(p)), None)
needs_units = pytest.mark.skipif(UNITS_CSV is None, reason="degenerate-set CSV not available")


@pytest.fixture(scope="module")
def sb():
    sc, te = TR.load_sb_modes(LIB)
    return sc, te


def test_staircase_matches_direct_count(sb):
    sc, _ = sb
    s, y = sc["V_Stretch"].to_numpy(float), sc["y"].to_numpy(float)
    E = TR.error_counts(s[None], y[None], np.ones((1, len(s))))[0]
    direct = np.array([np.sum(((y == 1) & (s < g)) | ((y == 0) & (s >= g))) for g in TR.GRID])
    assert np.array_equal(E, direct)


def test_optimum_rows_matches_loop():
    rng = np.random.default_rng(0)
    E = rng.integers(0, 3, size=(50, TR.NG)).astype(float)
    lo, hi, emin = TR.optimum_rows(E)
    for r in range(len(E)):
        ok = E[r] <= E[r].min() + 1e-9
        idx = np.flatnonzero(ok)
        b = idx[np.argmin(np.abs(TR.GRID[idx] - 0.5) + 1e-9 * TR.GRID[idx])]
        l_, h_ = b, b
        while l_ > 0 and ok[l_ - 1]:
            l_ -= 1
        while h_ < TR.NG - 1 and ok[h_ + 1]:
            h_ += 1
        assert (lo[r], hi[r]) == (l_, h_)


def test_sc_plateau_and_errors(sb):
    sc, te = sb
    assert len(sc) == 409 and sc["molecule"].nunique() == 66
    s, y = sc["V_Stretch"].to_numpy(float), sc["y"].to_numpy(float)
    res = TR.recalibrate(s, y, np.ones_like(s), GI)
    lo_g, hi_g = res["lo"][0], res["hi"][0]
    gap = (s[s < lo_g].max(), s[s >= hi_g].min())
    assert 0.407 < gap[0] < 0.408 and 0.444 < gap[1] < 0.445       # (0.40706, 0.44422]
    assert round(0.5 * (gap[0] + gap[1]), 2) == TAU                  # tau_SB = rounded midpoint
    assert res["emin"][0] == 1 and res["e_tau"][0] == 1             # 1/409 SC errors at 0.43
    st, yt = te["V_Stretch"].to_numpy(float), te["y"].to_numpy(float)
    assert int(np.sum((st >= TAU) != (yt == 1))) == 19              # per-member scores: 19/379


def test_d1_keeps_tau_optimal(sb):
    sc, _ = sb
    s, y = sc["V_Stretch"].to_numpy(float), sc["y"].to_numpy(float)
    mols = sc["molecule"].unique()
    codes = pd.Categorical(sc["molecule"], categories=mols).codes
    w = TR.molecule_weights(codes, len(mols), 2000, np.random.default_rng(TR.SEED_DLC))
    res = TR.recalibrate(s, y, w, GI)
    assert np.all(res["e_tau"] - res["emin"] < 1e-9)


@needs_units
def test_single_flips_and_ambiguity_range(tmp_path):
    for f in ("library_scores.csv", "thresholds.json", "C6H6_normal.csv"):
        shutil.copy(os.path.join("data", "results", f), tmp_path / f)
    out = TR.run(data_dir="data", out_dir=str(tmp_path), results_dir=str(tmp_path), n_rep=1000,
                 units_csv=UNITS_CSV, make_figures=False)
    assert out["tau_SB"] == TAU and out["sc_errors_at_tau"] == 1
    assert out["n_units"] == 332 and not out["mixed_label_sets"]
    assert out["l1_excluding_tau"] == ["SBr4 135"]
    summ = out["summary"].set_index("id")
    assert summ.loc["D1", "pct_tau_optimal"] == 100
    assert summ.loc["L1", "n_rep"] == 332
    lo, hi = out["ci95"]
    assert abs(lo - AMBIGUITY_RANGE[0]) < 0.02 and abs(hi - AMBIGUITY_RANGE[1]) < 0.02
    assert lo < TAU < hi
    for k in ("summary_csv", "ci_csv", "modes_csv", "tests_tex"):
        assert os.path.getsize(out["paths"][k]) > 0
    tex = open(out["paths"]["tests_tex"], encoding="utf-8").read()
    assert r"\begin{tabular}" in tex and "332 (all)" in tex
