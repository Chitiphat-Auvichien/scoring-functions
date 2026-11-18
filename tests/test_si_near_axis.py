"""Checks for src/si_near_axis.py (R1.6 SI figure: N of NH3 moved off the C3 axis, mode 3)."""

import os

import numpy as np
import pytest

from src import si_near_axis as sna
from src.scoring import ON_AXIS_TOL

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, "data")


@pytest.fixture(scope="module")
def maps():
    w = sna.nh3(DATA)
    df = sna.maps(w, per_dec_tol=2, per_dec_notol=1, n_theta=12)
    return w, df


def test_frame_and_matches_modescorer(maps):
    w, df = maps
    assert w["k"] == 2                                          # C3 axis = z in the program frame
    assert np.hypot(*w["X0"][w["n"], [w["ka"], w["kb"]]]) < 1e-15   # N put exactly on the axis
    assert sna.validate(w, df, n_sample=10) < 1e-10              # own eq:rscore = ModeScorer.calculate_scores


def test_inside_tolerance_score_is_invariant(maps):
    _, df = maps
    tol = df[df.tolerance == "1e-5 A"]
    inside = tol[tol.r_Q_N_A <= ON_AXIS_TOL]
    s_h = inside.s_R_z.iloc[0]
    assert (inside.contribution_N == 0).all() and (inside.N_minus_NQ == 3).all()
    np.testing.assert_allclose(inside.s_R_z, s_h, atol=1e-4)    # H atoms alone, wherever N sits inside
    assert abs(s_h - 0.171) < 1e-3


def test_outside_and_without_tolerance(maps):
    _, df = maps
    tol = df[df.tolerance == "1e-5 A"]
    s_h = tol[tol.r_Q_N_A <= ON_AXIS_TOL].s_R_z.iloc[0] * 3       # sum of the H summands
    out = tol[(tol.r_Q_N_A > ON_AXIS_TOL) & (tol.r_Q_N_A < 1e-4)]
    assert (out.N_minus_NQ == 4).all()
    np.testing.assert_allclose(out.s_R_z, (s_h + out.contribution_N) / 4, atol=1e-3)
    no = df[(df.tolerance == "none") & (df.rho_A > 0)]
    assert no.contribution_N.abs().max() <= 1 + 1e-12            # bounded at every scale
    assert no.contribution_N.max() > 0.9 and no.contribution_N.min() < -0.9   # full size, sign set by direction
