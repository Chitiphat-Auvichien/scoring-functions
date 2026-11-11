"""Fast checks of src/si_frame.py (SI analyses of the frame, referee comments R1.5 and R1.6).

Run from ``Github/scoring-functions/``:  py -m pytest tests/test_si_frame.py

Uses a small set of roster molecules covering every case (asymmetric, linear, symmetric,
spherical; a planar D3h test molecule whose allowed frames move an on-axis H off its axis).
"""
import os
import sys

import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
os.chdir(ROOT)

import main as M  # noqa: E402
from src import si_frame as SF  # noqa: E402
from src.scoring import EPS_DENOM, ON_AXIS_TOL, ModeScorer  # noqa: E402

MOLS = ["H2O", "CO2", "NH3", "B3N3H6", "CH4", "SnH4"]


@pytest.fixture(scope="module")
def tables():
    return SF.compute(molecules=MOLS, n_sphere=2000)


@pytest.fixture(scope="module")
def inputs():
    return {m: SF.molecule_inputs(m) for m in ("NH3", "CH4", "CO2")}


def test_cases_detected(tables):
    S, _, _ = tables
    case = S.drop_duplicates("molecule").set_index("molecule").case.to_dict()
    assert case == {"H2O": "asymmetric", "CO2": "linear", "NH3": "symmetric", "B3N3H6": "symmetric",
                    "CH4": "spherical", "SnH4": "spherical"}


def test_axis_scores_match_modescorer(inputs):
    """The vectorised per-direction scores equal ModeScorer in a frame (eq:tscore, eq:rscore)."""
    for d in inputs.values():
        F = d["F"]
        T, R = SF.axis_scores(d["X"], d["vecs"], F)
        sc = ModeScorer(d["sym"], d["X"] @ F.T, d["bonds"])
        for k, v in enumerate(d["vecs"]):
            r = sc.calculate_scores(v @ F.T)
            assert np.allclose(T[k], [r["T"][q] for q in "xyz"], atol=1e-12)
            assert np.allclose(R[k], [r["R"][q] for q in "xyz"], atol=1e-12)


def test_frame_is_the_programs(inputs):
    """The frame analysed is the one main.build_scorer_and_final applies."""
    for mol, d in inputs.items():
        raw, _ = M.load_inputs(mol, "normal")
        scorer, _ = M.build_scorer_and_final(raw, "normal")
        assert np.allclose(scorer.frame, d["F"], atol=1e-10)


def test_asymmetric_spread_is_zero(tables):
    S, _, _ = tables
    asym = S[S.case == "asymmetric"]
    assert len(asym) and (asym.spread == 0).all()


def test_spread_bounds(tables):
    S, _, _ = tables
    assert (S.spread <= S.max_abs + 1e-12).all() and (S.max_abs < 1).all()
    tops = S[S.case.isin(["symmetric", "spherical"])]
    assert (tops.spread_allowed <= tops.spread + 1e-9).all()       # allowed frames are a subset
    # 1-D irreps keep |s| over the allowed frames (B3N3H6: the R1.6 on-axis exception)
    one = tops[(tops.set_size == 1) & (tops.molecule != "B3N3H6")]
    assert one.spread_allowed.max() < 1e-4                         # modes are printed to 5 decimals


def test_linear_spread(tables):
    S, _, _ = tables
    v = S[(S.molecule == "CO2") & (S.max_abs >= SF.MIN_ABS)].spread
    assert len(v) == 2 and np.allclose(v, 1 / 3, atol=1e-3)


def test_selected_frame_is_allowed(tables):
    _, P, _ = tables
    n_sel = P[P.selected].groupby("molecule").frame_index.nunique()
    assert (n_sel == 1).all() and len(n_sel) == len(MOLS)


def test_onaxis_classes_separated_by_cutoff(tables):
    """On-axis (by symmetry) pairs lie within ON_AXIS_TOL, off-axis pairs outside, in the selected
    frame and in every allowed frame; the class itself never uses the cut-off."""
    _, P, _ = tables
    for Q in (P[P.selected], P):
        assert (Q[Q.on_by_symmetry].r_Q <= ON_AXIS_TOL).all()
        assert (Q[~Q.on_by_symmetry].r_Q > ON_AXIS_TOL).all()


def test_onaxis_known_pairs(tables):
    _, P, _ = tables
    sel = P[P.selected & P.on_by_symmetry]
    # NH3: N on the C3 axis; CO2: every atom on the molecular axis, C at the CM
    assert set(sel[sel.molecule == "NH3"].atom) == {"N1"}
    co2 = sel[sel.molecule == "CO2"]
    assert len(co2[co2.axis == "x"]) == 3                             # the molecular axis
    assert set(co2[co2.reason == "at the CM"].atom) == {"C1"} and len(co2[co2.atom == "C1"]) == 3
    # B3N3H6: an allowed (not selected) frame puts an on-axis H ~1.1e-6 A off its axis
    b = P[(P.molecule == "B3N3H6") & P.on_by_symmetry]
    assert 1e-6 < b.r_Q.max() < ON_AXIS_TOL


def test_displacement_regions(tables):
    _, P, D = tables
    n_expected = sum(len(P[(P.molecule == m) & P.selected]) * D[D.molecule == m]["mode"].nunique() for m in MOLS)
    assert len(D) == n_expected
    assert ((D.region == "not counted") == (D.r_Q <= ON_AXIS_TOL)).all()
    assert ((D.region == "not counted") == D.on_by_symmetry).all()
    t0 = D[D.region == "term 0"]
    assert (t0.d <= EPS_DENOM).all() and (D[D.region == "full term"].d > EPS_DENOM).all()


def test_roster_filter():
    r = SF.roster()
    full = pd.read_csv("data/mol_list_method.csv")
    assert set(r.mol_type) <= set(SF.MOL_TYPES)
    assert len(r) == int(full.mol_type.isin(SF.MOL_TYPES).sum())
