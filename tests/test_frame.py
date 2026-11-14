"""Tests for src/frame.py, the principal-axis frame rule (JCC revision, referee comment R1.5).

Run from ``Github/scoring-functions/``:  py -m pytest tests/test_frame.py

Reads the tracked intermediates in data/intermediate/. The checks follow the
rule's reference prototype (checks I1-I12 and S1-S5 of the R1.5 analysis):
  * top detection by symmetry agrees with the roster point groups (85/85);
  * every Gaussian standard orientation is an allowed frame, kept with keep_input_frame;
  * 3N input (vibrations + T/R vectors built in an allowed frame) returns that frame;
  * benzene: the 36 EMIT modes give the published T/R picks (axes x<->y vs the
    submitted MIT frame), and the normal-mode and EMIT pools end in the same frame;
  * re-oriented / renumbered input gives the same assignment, labels and |scores|;
  * re-oriented input with the same numbering gives exactly the same signed frame;
  * asymmetric frames are right-handed; linear molecules: the vector turn is orientation-independent,
    within one scan step of the scan, and follows the densest direction of the modes' sum_A d^A.
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
from src import frame as FR  # noqa: E402
from src.classifier import Thresholds, classify_all_modes  # noqa: E402
from src.parser import IntermediateIO  # noqa: E402
from src.scoring import ModeScorer  # noqa: E402

PG_TOP = {"C3v": "symmetric", "C4v": "symmetric", "D3h": "symmetric", "D4h": "symmetric",
          "D6h": "symmetric", "Td": "spherical", "Oh": "spherical", "Dinfh": "linear", "Cinfv": "linear"}
RNG = np.random.default_rng(20261111)


def _roster():
    r = pd.read_csv(os.path.join(ROOT, "data", "mol_list_method.csv"))
    return r[r["mol_type"].isin(["ideal", "non-ideal", "test", "multi-centre"])]


def _load(mol, mode="normal"):
    return IntermediateIO.load(M.intermediate_path(M.resolve_dirs("data"), mol, mode))


def _com(raw):
    sc = ModeScorer(raw["atoms"], raw["coords"], raw["bonds"])
    return sc.coords.copy(), sc.symbols, sc.masses, raw["bonds"]


def _vecs(raw):
    return [np.asarray(m["vector"], dtype=float) for m in raw["modes"]]


def _random_rotation():
    Q, R = np.linalg.qr(RNG.normal(size=(3, 3)))
    Q = Q @ np.diag(np.sign(np.diag(R)))
    return Q if np.linalg.det(Q) > 0 else -Q


def _classify(X, sym, bonds, vecs, f, mode_type="normal"):
    sc = ModeScorer(sym, X @ f.T, bonds)
    vib = [dict(frequency=0.0, vector=v @ f.T, label=f"Vib {k + 1}") for k, v in enumerate(vecs)]
    final = sc.construct_T() + sc.construct_R() + vib if mode_type == "normal" else vib
    return classify_all_modes(sc, final, Thresholds())


TOPS = [r.molecule for r in _roster().itertuples() if PG_TOP.get(r.point_group) in ("symmetric", "spherical")]


def test_detect_top_matches_roster_point_groups():
    bad = []
    for r in _roster().itertuples():
        X, sym, m, _ = _com(_load(r.molecule))
        kind, _, _ = FR.detect_top(X, sym, m)
        if kind != PG_TOP.get(r.point_group, "asymmetric"):
            bad.append((r.molecule, r.point_group, kind))
    assert not bad
    assert len(TOPS) == 40


def test_gaussian_standard_orientation_is_allowed_and_kept():
    for mol in TOPS:
        raw = _load(mol)
        X, sym, m, bonds = _com(raw)
        kind, z, _ = FR.detect_top(X, sym, m)
        cands = FR.candidate_frames(X, sym, m, kind, z)
        assert any(np.allclose(c, np.eye(3), atol=1e-4) for c in cands), mol
        f, info = FR.select_frame(X, sym, m, bonds, _vecs(raw), keep_input_frame=True)
        assert np.allclose(f, np.eye(3), atol=1e-4), mol
        assert info["step"].endswith("input frame kept")


@pytest.mark.parametrize("tie_break", FR.TIE_BREAKS)
@pytest.mark.parametrize("mol", ["NH3", "BF3", "C6H6", "CH3CN"])
def test_3n_input_returns_the_frame_its_references_were_built_in(mol, tie_break):
    raw = _load(mol)
    X, sym, m, bonds = _com(raw)
    kind, z, _ = FR.detect_top(X, sym, m)
    for f0 in FR.candidate_frames(X, sym, m, kind, z)[::4]:   # one sign pattern per atom line
        s0 = ModeScorer(sym, X @ f0.T, bonds)
        tr = [md["vector"] @ f0 for md in s0.construct_T() + s0.construct_R()]   # back to input coords
        f, _ = FR.select_frame(X, sym, m, bonds, _vecs(raw) + tr, tie_break=tie_break)
        assert FR.line_set(f) == FR.line_set(f0), mol


def test_benzene_emit_picks_and_common_frame():
    raw_n, _ = M.load_inputs("C6H6", "normal")
    raw_e, _ = M.load_inputs("C6H6", "emit")
    sc_n, final_n = M.build_scorer_and_final(raw_n, "normal")
    sc_e, final_e = M.build_scorer_and_final(raw_e, "emit")
    assert np.allclose(sc_n.coords, sc_e.coords, atol=1e-6)       # one frame for the projection
    assert np.allclose(sc_e.frame, np.eye(3), atol=1e-4)          # Gaussian's standard orientation
    picks = {r["tr_label"]: r["name"] for r in classify_all_modes(sc_e, final_e, Thresholds()) if r["tr_label"]}
    assert picks == {"Tx": "EMIT 35", "Ty": "EMIT 34", "Tz": "EMIT 36",
                     "Rx": "EMIT 9", "Ry": "EMIT 6", "Rz": "EMIT 19"}
    assert sc_e.frame_info["tie_break"] == "tmax"
    assert sc_e.frame_info["values"][0] == pytest.approx(2.0, abs=1e-3)          # EMIT 34 + 35 on Tx, Ty
    X, sym, m, bonds = sc_e.coords @ sc_e.frame, sc_e.symbols, sc_e.masses, sc_e.bList   # back to the input frame
    vecs = [md["vector"] @ sc_e.frame for md in final_e]
    f_gap, info_gap = FR.select_frame(X, sym, m, bonds, vecs, tie_break="gap")
    assert np.allclose(f_gap, sc_e.frame, atol=1e-6)                             # both rules: the same frame
    assert info_gap["values"][0] == pytest.approx(5.7875, abs=1e-3)


def test_tie_break_default_and_option():
    assert FR.get_tie_break() == "tmax"
    raw = _load("PCl3")                                   # 3N-6 input: the two rules choose different lines
    X, sym, m, bonds = _com(raw)
    f_t, i_t = FR.select_frame(X, sym, m, bonds, _vecs(raw))
    f_g, i_g = FR.select_frame(X, sym, m, bonds, _vecs(raw), tie_break="gap")
    assert i_t["tie_break"] == "tmax" and i_g["tie_break"] == "gap"
    assert FR.line_set(f_t) != FR.line_set(f_g)
    assert i_g["values"][0] == pytest.approx(1.5438, abs=1e-3)
    FR.set_tie_break("gap")
    try:
        f_g2, _ = FR.select_frame(X, sym, m, bonds, _vecs(raw))
    finally:
        FR.set_tie_break("tmax")
    assert np.allclose(f_g2, f_g)
    with pytest.raises(ValueError):
        FR.set_tie_break("sum")


def test_t_max_is_the_best_assignment_of_distinct_modes():
    raw = _load("NH3")
    X, sym, m, bonds = _com(raw)
    vecs = _vecs(raw)
    sc = ModeScorer(sym, X, bonds)
    A = np.array([[abs(sc.calculate_scores(v)["T"][q]) for q in "xy"] for v in vecs])
    brute = max(A[a, 0] + A[b, 1] for a in range(len(A)) for b in range(len(A)) if a != b)
    assert FR.t_max(X, sym, bonds, vecs, np.eye(3)) == pytest.approx(brute, abs=1e-12)


@pytest.mark.parametrize("mol", ["NH3", "BF3", "C6H6", "CH4", "PCl5", "XeOH4"])
def test_reoriented_and_renumbered_input_same_assignment_labels_scores(mol):
    raw = _load(mol)
    X, sym, m, bonds = _com(raw)
    vecs = _vecs(raw)

    def summary(X, sym, bonds, vecs):
        mm = ModeScorer(sym, X, bonds).masses
        f, _ = FR.select_frame(X, sym, mm, bonds, vecs)
        res = _classify(X, sym, bonds, vecs, f)
        assign = sorted((r["tr_label"], r["name"]) for r in res if r["tr_label"])
        labels = [(r["name"], r["vib_label"]) for r in res]
        mags = np.array([sorted([abs(r["T"][q]) for q in "xyz"]) + sorted([abs(r["R"][q]) for q in "xyz"]) + [r["V"]]
                         for r in res])
        return assign, labels, mags

    a0, l0, s0 = summary(X, sym, bonds, vecs)
    Q = _random_rotation()
    p = RNG.permutation(len(X))
    inv = np.argsort(p)
    a1, l1, s1 = summary((X @ Q.T)[p], [sym[i] for i in p], [(int(inv[a]), int(inv[b])) for a, b in bonds],
                         [(v @ Q.T)[p] for v in vecs])
    assert a0 == a1 and l0 == l1
    assert np.abs(s1 - s0).max() < 1e-6


def test_reoriented_input_same_signed_frame_for_every_nonlinear_molecule():
    bad = []
    for r in _roster().itertuples():
        if PG_TOP.get(r.point_group) == "linear":
            continue
        raw = _load(r.molecule)
        X, sym, m, bonds = _com(raw)
        f0, _ = FR.select_frame(X, sym, m, bonds, _vecs(raw))
        Q = _random_rotation()
        f1, _ = FR.select_frame(X @ Q.T, sym, m, bonds, [v @ Q.T for v in _vecs(raw)])
        if not np.allclose(f1 @ Q, f0, atol=1e-4):
            bad.append(r.molecule)
    assert not bad


def test_asymmetric_frames_are_right_handed_on_the_eigenvector_lines():
    for r in _roster().itertuples():
        if PG_TOP.get(r.point_group, "asymmetric") != "asymmetric":
            continue
        raw = _load(r.molecule)
        X, sym, m, bonds = _com(raw)
        f, info = FR.select_frame(X, sym, m, bonds, _vecs(raw))
        assert info["kind"] == "asymmetric"
        assert np.linalg.det(f) > 0
        _, W = np.linalg.eigh(FR.inertia(X, m))
        assert np.allclose(np.abs(f @ W), np.eye(3), atol=1e-6), r.molecule   # x = a, y = b, z = c


def test_linear_molecule_scan_reproducible_and_3n_returns_its_frame():
    raw = _load("CO2")
    X, sym, m, bonds = _com(raw)
    vecs = _vecs(raw)
    f0, info = FR.select_frame(X, sym, m, bonds, vecs)
    assert info["kind"] == "linear" and np.linalg.det(f0) > 0
    assert np.allclose(np.abs(f0[0]), np.abs(FR._linear_start(X, m)[0]), atol=1e-12)   # x = molecular axis
    Q = _random_rotation()
    f1, _ = FR.select_frame(X @ Q.T, sym, m, bonds, [v @ Q.T for v in vecs])
    assert np.allclose(f1 @ Q, f0, atol=1e-6)                                        # same names and signs
    s0 = ModeScorer(sym, X @ f0.T, bonds)
    tr = [md["vector"] @ f0 for md in s0.construct_T() + s0.construct_R()]
    f3, _ = FR.select_frame(X, sym, m, bonds, vecs + tr)
    assert np.allclose(np.abs(f3 @ f0.T).max(axis=1), 1, atol=1e-6)


@pytest.mark.parametrize("mol", ["CO2", "CS2", "CSe2", "CTe2"])
def test_linear_vector_turn_follows_the_densest_direction_and_matches_the_scan(mol):
    raw = _load(mol)
    X, sym, m, bonds = _com(raw)
    vecs = _vecs(raw)
    fv, iv = FR.frame_linear(X, sym, m, bonds, vecs, linear_rule="vector")
    fs, _ = FR.frame_linear(X, sym, m, bonds, vecs, linear_rule="scan")
    assert iv["step"].startswith("linear: vector")
    assert np.all(np.sum(fv * fs, axis=1) > np.cos(np.radians(FR.LINEAR_STEP_DEG)))  # within one scan step
    # 3N input with the T/R vectors built at 62 deg: the densest direction (T, |w| = sqrt(N)) wins
    # over the bends, where a plain 4-fold mean of all votes would land in between
    ax, u = FR._linear_start(X, m)
    f62 = FR._linear_frame_at(ax, u, np.radians(62.0))
    s0 = ModeScorer(sym, X @ f62.T, bonds)
    tr = [md["vector"] @ f62 for md in s0.construct_T() + s0.construct_R()]
    f3, i3 = FR.frame_linear(X, sym, m, bonds, vecs + tr, linear_rule="vector")
    assert abs(i3["phi"] - 62.0) < 1e-6
    assert np.allclose(np.abs(f3 @ f62.T).max(axis=1), 1, atol=1e-9)


def test_eigh_rule_restores_the_submitted_frame():
    raw, _ = M.load_inputs("C6H6", "normal")
    FR.set_frame_rule("eigh")
    try:
        sc, _ = M.build_scorer_and_final(raw, "normal")
    finally:
        FR.set_frame_rule("symmetry")
    ref = ModeScorer(raw["atoms"], raw["coords"], raw["bonds"])
    ref.MIT()
    assert np.allclose(sc.coords, ref.coords, atol=1e-12)
