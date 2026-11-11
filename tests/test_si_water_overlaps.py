"""Checks for src/si_water_overlaps.py (R1.8 SI figure: water MS scores vs o_Cart and o_MW)."""

import os

import numpy as np
import pytest

from main import load_inputs, build_scorer_and_final
from src import si_water_overlaps as swo

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(REPO, "data")


@pytest.fixture(scope="module")
def water():
    raw, _ = load_inputs("H2O", "normal", data_dir=DATA)
    scorer, final = build_scorer_and_final(raw, "normal")
    df = swo.build_table(scorer, final, "H2O").set_index("mode")
    return scorer, final, df


def _block(df, q, rows):
    return df.loc[rows, [f"{q}_{s}" for s in swo.SLOTS]].to_numpy(float)


def test_rows_and_columns(water):
    _, _, df = water
    assert list(df.index) == swo.SLOTS + ["delta", "nu_s", "nu_as"]
    for q in swo.QUANTITIES:
        for s in swo.SLOTS:
            assert f"{q}_{s}" in df.columns


def test_omw_ideal_block_is_identity_and_vibrations_vanish(water):
    _, _, df = water
    np.testing.assert_allclose(_block(df, "oMW", swo.SLOTS), np.eye(6), atol=1e-10)
    # Eckart conditions: normal vibrations carry no mass-weighted T/R component.
    assert np.abs(_block(df, "oMW", ["delta", "nu_s", "nu_as"])).max() < 1e-3


def test_ms_matches_table2_up_to_sign(water):
    _, _, df = water
    assert abs(df.loc["Tz", "MS_Rx"]) == pytest.approx(1 / 3, abs=5e-3)
    assert abs(df.loc["Rx", "MS_Tz"]) == pytest.approx(1 / 3, abs=5e-3)
    assert abs(df.loc["Tx", "MS_Rz"]) == pytest.approx(0.04, abs=5e-3)
    assert abs(df.loc["delta", "MS_Ty"]) == pytest.approx(0.23, abs=5e-3)
    assert abs(df.loc["nu_as", "MS_Tx"]) == pytest.approx(0.20, abs=5e-3)
    assert abs(df.loc["nu_as", "MS_Rz"]) == pytest.approx(0.30, abs=5e-3)


def test_ocart_cross_overlaps_match_slide(water):
    _, _, df = water
    assert abs(df.loc["Tz", "oCart_Rx"]) == pytest.approx(0.76, abs=5e-3)
    assert abs(df.loc["Tx", "oCart_Rz"]) == pytest.approx(0.43, abs=5e-3)
    assert abs(df.loc["nu_as", "oCart_Tx"]) == pytest.approx(0.61, abs=5e-3)
    np.testing.assert_allclose(np.diag(_block(df, "oCart", swo.SLOTS)), 1.0, atol=1e-10)


def _general_score_m2(scorer, d):
    """Side-report eq. general with w^A = m^A, moment m = 2, written out independently:
    T = P.Q / (M |d|_rms),  R = L_Q / (I_Q (sum m |d|^2 / I_Q)^(1/2))."""
    m = np.asarray(scorer.masses, float)
    r = np.asarray(scorer.coords, float)
    r = r - (m[:, None] * r).sum(0) / m.sum()
    M = m.sum()
    amp = np.sqrt(np.sum(m * np.sum(d ** 2, axis=1)) / M)
    T = (m @ d) / (M * amp)
    R = []
    for q in range(3):
        rp = r.copy()
        rp[:, q] = 0.0
        IQ = np.sum(m * np.sum(rp ** 2, axis=1))
        LQ = np.sum(m * np.cross(r, d)[:, q])
        R.append(LQ / (IQ * np.sqrt(np.sum(m * np.sum(d ** 2, axis=1)) / IQ)))
    return np.concatenate([T, R])


def test_omw_equals_m2_general_score(water):
    scorer, final, _ = water
    rng = np.random.default_rng(0)
    vecs = [np.asarray(md["vector"], float) for md in final]
    vecs += [rng.normal(size=vecs[0].shape) for _ in range(5)]
    refs = swo.ideal_references(scorer)
    m = np.asarray(scorer.masses, float)
    for d in vecs:
        o = [swo.weighted_overlap(d, e, m) for e in refs]
        np.testing.assert_allclose(o, _general_score_m2(scorer, d), atol=1e-10)


def test_weighted_overlap_zero_vector():
    assert swo.weighted_overlap(np.zeros((3, 3)), np.ones((3, 3)), np.ones(3)) == 0.0
