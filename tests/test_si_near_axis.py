"""Checks for src/si_near_axis.py (R1.6 SI figure: s[R_C2] of water as O is moved off the C2 axis)."""

import os

import numpy as np
import pytest

from src import si_near_axis as sna
from src.scoring import ON_AXIS_TOL

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, "data")


@pytest.fixture(scope="module")
def scan():
    w = sna.water(DATA)
    df, dev = sna.scan(w, n_r=15, n_psi=9, validate_every=1)
    return w, df, dev


def test_frame_and_matches_modescorer(scan):
    w, df, dev = scan
    assert w["r_O_axis"] < 1e-12                      # O exactly on the C2 axis before it is moved
    assert w["out_of_plane_max"] < 1e-12              # planar molecule, out-of-plane axis found
    assert dev < 1e-10                                # own eq:rscore = ModeScorer.calculate_scores
    assert np.allclose(df.r_Q_O, df.r_target, rtol=1e-6)


def test_step_at_the_cutoff(scan):
    _, df, _ = scan
    below = df[df.r_Q_O <= ON_AXIS_TOL]
    near = df[(df.r_Q_O > ON_AXIS_TOL) & (df.r_Q_O <= 1e-4)]
    assert below.s_nu_as.abs().max() < 1e-3
    sinp = np.abs(np.sin(np.radians(near.psi_deg)))
    np.testing.assert_allclose(3 * near.s_nu_as.abs(), sinp, atol=2e-3)   # O's full term: sin(psi)/3
    assert df[["s_delta", "s_nu_s"]].abs().to_numpy().max() < 5e-3      # O moves along the axis: term 0


def test_ideal_rotation_and_no_nq(scan):
    w, df, _ = scan
    np.testing.assert_allclose(df.s_R_C2, 1.0, atol=1e-12)
    assert abs(sna.no_nq_rotation(w) - 2.0 / 3.0) < 1e-12


def test_recomputed_axis_follows_o_out_of_plane(scan):
    w, _, _ = scan
    rec = sna.recomputed_axis_distance(w)
    assert rec["out_of_plane"] < 1e-12 and rec["in_plane"] > 1e-4
