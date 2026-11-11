"""Tests for src/si_degenerate.py (SI, referee comment R1.4: s[V_S] inside degenerate sets).

Run from ``Github/scoring-functions/``:  py -m pytest tests/test_si_degenerate.py

Reads the tracked intermediates in data/intermediate/. Reference values are the R1.4 side report's
degenerate_sets.csv (side_reports/R1.4_degenerate_modes, full-resolution scan); the reduced scan here
(200 rotations, short optimisations) must reproduce them within the scan tolerance.
"""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
DATA = os.path.join(ROOT, "data")

from src import si_degenerate as SD  # noqa: E402

TAU = 0.43
# side report degenerate_sets.csv: Vmin, Vmax, V_mean, member_avg_min, member_avg_max, w_max
SIDE_REPORT = {
    ("C6H6", 1056.39): (0.441343, 0.509063, 0.486575, 0.475219, 0.492299, 0.067707),
    ("C4H4", 1123.58): (0.380425, 0.515003, 0.481263, 0.380428, 0.506214, 0.081579),
}


@pytest.fixture(scope="module")
def all_sets():
    diag = {}
    sets = SD.find_degenerate_sets(DATA, diagnostics=diag)
    return sets, diag


@pytest.fixture(scope="module")
def scanned():
    sets = SD.find_degenerate_sets(DATA, molecules=["C6H6", "C4H4"], keep_modes=True)
    pick = [s for s in sets if any(s["molecule"] == m and abs(s["freq"] - f) < 0.5 for m, f in SIDE_REPORT)]
    D, samples = SD.scan_sets(pick, TAU, n_rot=200, n_start=20, n_refine=3)
    return D, samples


def test_152_sets_by_group_and_class(all_sets):
    sets, _ = all_sets
    assert len(sets) == 152
    groups = [s["group"] for s in sets]
    assert (groups.count("test"), groups.count("SC"), groups.count("benzene")) == (71, 71, 10)
    classes = [s["sym_class"] for s in sets]
    assert {c: classes.count(c) for c in SD.SYM_CLASSES} == {"E3": 113, "T": 30, "E4": 5, "Pi": 4}
    assert all(s["deg"] == len(s["member_index"]) == len(s["irreps"].split("/")) for s in sets)
    assert all(s["freq_spread"] < SD.TOL_CM * (s["deg"] - 1) for s in sets)


def test_alcl3_unlabelled_pairs_resolved_by_c3(all_sets):
    sets, diag = all_sets
    assert len(diag["c3_character"]) == 2
    assert all(abs(x["c3_character"] + 1) < 1e-3 for x in diag["c3_character"])
    unl = [s for s in sets if s["status"] == "symmetry-unlabelled"]
    assert [(s["molecule"], round(s["freq"])) for s in unl] == [("AlCl3", 140), ("AlCl3", 586)]
    # dropped close clusters never carry a degenerate irrep
    for x in diag["close_pairs_not_degenerate"]:
        assert all((SD.irrep_dim(lab) or 1) == 1 for lab in x["irreps"].split("/"))


def test_benzene_1056_members(scanned):
    D, _ = scanned
    r = D[(D.molecule == "C6H6") & (abs(D.freq - 1056.39) < 0.5)].iloc[0]
    assert sorted(np.round(SD._members(r), 2)) == [0.45, 0.51]       # manuscript Table 3
    assert r.irreps == "E1U/E1U" and r.sym_class == "E3"


def test_member_mean_and_invariants(scanned):
    D, samples = scanned
    for (_, r), (w, dm) in zip(D.iterrows(), samples):
        m = SD._members(r)
        assert r.member_avg_orig == pytest.approx(m.mean(), abs=1e-6)     # member mean = mean of members
        assert r.w_gauss == pytest.approx(np.ptp(m), abs=1e-6)
        assert r.vscan_dev < 1e-10                                         # VScan == ModeScorer.Vscore
        assert r.Vmin <= m.min() + 1e-9 and r.Vmax >= m.max() - 1e-9
        assert r.Vmin - 1e-9 <= r.member_avg_min <= r.member_avg_max <= r.Vmax + 1e-9
        assert r.w_max <= r.spread + 1e-6 and r.w_min <= r.w_gauss + 1e-6
        assert len(w) == len(dm) == r.n_rot


def test_reduced_scan_reproduces_side_report(scanned):
    D, _ = scanned
    for (mol, f0), (vmin, vmax, vmean, mmlo, mmhi, wmax) in SIDE_REPORT.items():
        r = D[(D.molecule == mol) & (abs(D.freq - f0) < 0.5)].iloc[0]
        assert r.Vmin == pytest.approx(vmin, abs=2e-5)
        assert r.Vmax == pytest.approx(vmax, abs=2e-5)
        assert r.V_mean == pytest.approx(vmean, abs=2e-4)
        assert r.member_avg_min == pytest.approx(mmlo, abs=2e-3)
        assert r.member_avg_max == pytest.approx(mmhi, abs=2e-3)
        assert r.w_max == pytest.approx(wmax, abs=2e-3)
        assert r.meets_ambiguity_range
    # at tau_SB = 0.43 the T2 triple straddles the threshold; benzene 1056 (0.441-0.509) reads S in every basis
    by = {r.molecule: r for _, r in D.iterrows()}
    assert by["C4H4"].straddles_tau and by["C4H4"].labels_orig == "S/B/S"
    assert not by["C6H6"].straddles_tau and by["C6H6"].labels_orig == "S/S"
