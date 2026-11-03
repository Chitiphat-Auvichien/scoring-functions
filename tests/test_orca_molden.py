# [EDIT 2026-11-02 | by Claude (AI assistant) for Than (ThanchonBK)]
# Purpose: new tests for the ORCA .hess / Molden input path (src/orca_molden.py,
#          QCFileParser, the .bonds file, and the lossless intermediate cache).
# Previous version: none -- this whole file is new.
# [END EDIT NOTE]
"""ORCA .hess and Molden input to the pipeline (src/orca_molden.py, QCFileParser).

The same molecule must score the same whichever program wrote it. Two kinds
of check, matching the web app's suite this module was ported with:

  * real files (tests/data/cclib/: divinylbenzene, B3LYP/STO-3G, from
    Gaussian 16, ORCA 6 and a Molden writer -- see its README) run through
    main.load_inputs() in a throwaway data/ tree, cache included;
  * synthetic .hess/.molden files written from data/logs/HOCl.log and
    C10H8.log, which must reproduce the Gaussian path's scores exactly.

Run from ``scoring-functions/``:   python -m pytest tests/test_orca_molden.py
"""
import contextlib
import io
import os
import shutil
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
CCLIB = os.path.join(ROOT, "tests", "data", "cclib")
DATA = os.path.join(ROOT, "data")

import main                                                          # noqa: E402
from src.orca_molden import (BOHR_TO_ANGSTROM, ParseError,            # noqa: E402
                             parse_frequency_file, parse_molden,
                             parse_orca_hess, read_bond_file)
from src.classifier import classify_all_modes, classify_to_rows     # noqa: E402
from src.parser import GaussianParser, IntermediateIO                # noqa: E402
from src.utils import atomicMass                                     # noqa: E402


def _read(name):
    with open(os.path.join(CCLIB, name)) as f:
        return f.read()


def _scores_v(raw):
    """(frequency, V_S, S/B label) per real vibration, frequency-sorted -- the
    rows the CLI writes to data/results/<mol>_normal.csv."""
    scorer, final = main.build_scorer_and_final(raw, "normal")
    rows = classify_to_rows(classify_all_modes(scorer, final))
    vib = [r for r in rows if str(r["Mode"]).startswith("Vib")]
    return sorted(((r["Freq"], r["V_Stretch"], r["vib_label"]) for r in vib), key=lambda t: t[0])


@pytest.fixture
def data_tree(tmp_path):
    """A throwaway data/ tree; load_inputs() resolves 'data' relative to it."""
    for d in ("logs", "gjf"):
        (tmp_path / "data" / d).mkdir(parents=True)

    def add(stem, freq_file, bonds_file):
        ext = os.path.splitext(freq_file)[1]
        shutil.copy(os.path.join(CCLIB, freq_file), tmp_path / "data" / "logs" / f"{stem}{ext}")
        shutil.copy(os.path.join(CCLIB, bonds_file), tmp_path / "data" / "gjf" / f"{stem}.bonds")
    return tmp_path, add


def _load(tmp, stem):
    with contextlib.redirect_stdout(io.StringIO()):
        raw, _ = main.load_inputs(stem, "normal", str(tmp / "data"))
    return raw


# ----------------------------------------------------------------------
# Real files through load_inputs()
# ----------------------------------------------------------------------
@pytest.mark.parametrize("freq_file, bonds_file, max_dv", [
    ("dvb_ir_orca6.hess", "dvb_orca6.bonds", 0.01),   # same method as the g16 run
    ("dvb_ir.molden", "dvb_molden.bonds", 0.03),      # a separate calculation
])
def test_real_file_scores_like_gaussian(data_tree, freq_file, bonds_file, max_dv):
    tmp, add = data_tree
    add("g16", "dvb_ir_g16.out", "dvb_g16.bonds")
    add("other", freq_file, bonds_file)
    g = _load(tmp, "g16")
    o = _load(tmp, "other")
    assert len(o["modes"]) == len(g["modes"]) == 54
    assert len(o["bonds"]) == len(g["bonds"]) == 20
    sg, so = _scores_v(g), _scores_v(o)
    assert [lab for _, _, lab in so] == [lab for _, _, lab in sg]
    assert max(abs(a[1] - b[1]) for a, b in zip(sg, so)) < max_dv


def test_cached_run_equals_first_run(data_tree):
    # The intermediate cache must not truncate ORCA's extra digits: T/R move
    # on a 1e-6 A geometry change for symmetric molecules.
    tmp, add = data_tree
    add("dvb", "dvb_ir_orca6.hess", "dvb_orca6.bonds")
    first = _load(tmp, "dvb")
    cached = _load(tmp, "dvb")
    assert np.allclose(first["coords"], cached["coords"], rtol=0, atol=1e-13)
    for a, b in zip(first["modes"], cached["modes"]):
        assert np.allclose(a["vector"], b["vector"], rtol=0, atol=1e-13)
        assert a["reduced_mass"] is None and b["reduced_mass"] is None
    for (fa, va, la), (fb, vb, lb) in zip(_scores_v(first), _scores_v(cached)):
        assert fa == fb and la == lb and abs(va - vb) < 1e-12


def test_missing_bonds_file_gives_no_bonds(data_tree):
    tmp, add = data_tree
    add("dvb", "dvb_ir_orca6.hess", "dvb_orca6.bonds")
    os.remove(tmp / "data" / "gjf" / "dvb.bonds")
    raw = _load(tmp, "dvb")
    assert raw["bonds"] == []
    with pytest.raises(ValueError, match="bonds list"):
        main.build_scorer_and_final(raw, "normal")


def test_gaussian_intermediate_format_unchanged():
    # Gaussian data needs exactly 6 dp (geometry) / 5 dp (vectors): the
    # adaptive writer must reproduce the historical fixed-width layout.
    raw = GaussianParser(os.path.join(DATA, "logs", "HOCl.log")).parse()
    buf = os.path.join(os.path.dirname(CCLIB), "_tmp_HOCl_normal_data.txt")
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            IntermediateIO.save(raw, buf)
        with open(buf) as f:
            lines = f.read().splitlines()
    finally:
        if os.path.exists(buf):
            os.remove(buf)
    coord = [ln for ln in lines if ln.startswith("      1          8")][0]
    assert coord == f"{1:>7d}{8:>11d}{0:>12d}" + "".join(
        f"{v:>{w}.6f}" for v, w in zip(raw["coords"][0], (16, 12, 12)))
    vec = [ln for ln in lines if ln.startswith("   1     1     8")][0]
    assert len(vec) == 16 + 10 * 3 and vec.count(".") == 3


# ----------------------------------------------------------------------
# Synthetic files from real Gaussian logs
# ----------------------------------------------------------------------
def _gaussian(stem):
    return GaussianParser(os.path.join(DATA, "logs", f"{stem}.log")).parse()


def _unit(v, masses=None):
    x = np.asarray(v, float)
    if masses is not None:
        x = x * np.sqrt(masses)[:, None]
    return x / np.linalg.norm(x)


def _hess_text(raw, ncols=5, fmt="{:11.6f}", n_zero=None):
    n = len(raw["atoms"])
    nc = 3 * n
    nz = nc - len(raw["modes"]) if n_zero is None else n_zero
    vecs = [np.zeros((n, 3))] * nz + [_unit(m["vector"]) for m in raw["modes"]]
    freqs = [0.0] * nz + [m["frequency"] for m in raw["modes"]]
    mat = np.array([x.ravel() for x in vecs]).T
    out = ["$orca_hessian_file", "", "$vibrational_frequencies", str(len(freqs))]
    out += [f"{k:5d}   {f:14.6f}" for k, f in enumerate(freqs)]
    out += ["", "$normal_modes", f"{nc} {mat.shape[1]}"]
    for c0 in range(0, mat.shape[1], ncols):
        cols = range(c0, min(c0 + ncols, mat.shape[1]))
        out.append("        " + "".join(f"{c:11d}" for c in cols))
        for r in range(nc):
            out.append(f"{r:6d}  " + "".join(fmt.format(mat[r, c]) for c in cols))
    out += ["", "#", "# The atoms: label  mass x y z (in bohrs)", "#", "$atoms", str(n)]
    for a, xyz in zip(raw["atoms"], raw["coords"]):
        b = [c / BOHR_TO_ANGSTROM for c in xyz]
        out.append(f" {a:2s} {atomicMass[a.lower()]:10.4f} {b[0]:15.9f} {b[1]:15.9f} {b[2]:15.9f}")
    return "\n".join(out + ["", "$end", ""])


def _molden_text(raw, mass_weighted=False):
    masses = np.array([atomicMass[a.lower()] for a in raw["atoms"]]) if mass_weighted else None
    out = ["[Molden Format]", "[FREQ]"] + [f"{m['frequency']:12.4f}" for m in raw["modes"]]
    out += ["[FR-COORD]"] + [f"{a} " + " ".join(f"{c / BOHR_TO_ANGSTROM:16.10f}" for c in xyz)
                             for a, xyz in zip(raw["atoms"], raw["coords"])]
    out += ["[FR-NORM-COORD]"]
    for k, m in enumerate(raw["modes"], 1):
        out.append(f"vibration {k}")
        out += [" ".join(f"{c:14.8f}" for c in row) for row in _unit(m["vector"], masses)]
    return "\n".join(out) + "\n"


def _same_v_and_labels(raw, g):
    g = dict(g, bonds=raw["bonds"])
    a, b = _scores_v(raw), _scores_v(g)
    assert [x[2] for x in a] == [x[2] for x in b]
    assert max(abs(x[1] - y[1]) for x, y in zip(a, b)) < 2e-4


@pytest.mark.parametrize("stem", ["HOCl", "C10H8"])
@pytest.mark.parametrize("fmt", ["hess", "hess-E6", "molden", "molden-mw"])
def test_synthetic_matches_gaussian(stem, fmt):
    raw = _gaussian(stem)
    text = {"hess": lambda: _hess_text(raw),
            "hess-E6": lambda: _hess_text(raw, ncols=6, fmt="{:17.9E}"),
            "molden": lambda: _molden_text(raw),
            "molden-mw": lambda: _molden_text(raw, mass_weighted=True)}[fmt]()
    g = parse_frequency_file(text)
    assert g["atoms"] == raw["atoms"]
    assert np.allclose(g["coords"], raw["coords"], atol=1e-6)
    _same_v_and_labels(raw, g)


def test_mass_weighted_vectors_are_converted():
    raw = _gaussian("HOCl")
    g = parse_molden(_molden_text(raw, mass_weighted=True))
    for got, m in zip(g["modes"], raw["modes"]):
        assert np.allclose(got["vector"], _unit(m["vector"]), atol=1e-6)


def test_linear_cartesian_is_left_alone():
    # Symmetry makes linear-molecule modes orthogonal under BOTH metrics.
    raw = _gaussian("CO2")
    g = parse_molden(_molden_text(raw))
    for got, m in zip(g["modes"], raw["modes"]):
        assert np.allclose(got["vector"], _unit(m["vector"]), atol=1e-7)
    assert len(parse_orca_hess(_hess_text(raw))["modes"]) == 4


def test_zero_frequency_vibration_after_externals_is_kept():
    text = _hess_text(_gaussian("HOCl"))
    first_vib = text.split("    6   ")[1].split("\n")[0].strip()
    assert len(parse_orca_hess(text.replace(first_vib, "0.000000", 1))["modes"]) == 3


@pytest.mark.parametrize("text, match", [
    ("$orca_hessian_file\n$atoms\n1\n H 1.0 0 0 0\n$end\n", r"\$vibrational_frequencies"),
    ("[Molden Format]\n[Atoms] Angs\nH 1 1 0 0 0\n", r"\[FREQ\]"),
    ("just some text\n", "Not recognised"),
])
def test_missing_sections_are_refused(text, match):
    with pytest.raises(ParseError, match=match):
        parse_frequency_file(text)


@pytest.mark.parametrize("bad", ["nan", "inf", "abc"])
def test_junk_numbers_are_value_errors(bad):
    raw = _gaussian("HOCl")
    text = _hess_text(raw)
    target = text.split("$atoms")[1].split("\n")[2].split()[2]   # first x coordinate
    with pytest.raises(ValueError):
        parse_orca_hess(text.replace(target, bad, 1))


def test_wrong_number_of_zero_modes_is_refused():
    with pytest.raises(ParseError):
        parse_orca_hess(_hess_text(_gaussian("HOCl"), n_zero=4))


def test_non_orthogonal_vectors_are_refused():
    raw = _gaussian("HOCl")
    text = _molden_text(raw).split("[FR-NORM-COORD]")[0] + "[FR-NORM-COORD]\n"
    for k in range(3):
        text += f"vibration {k + 1}\n1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n"
    with pytest.raises(ParseError, match="not mutually orthogonal"):
        parse_molden(text)


# ----------------------------------------------------------------------
# .bonds files
# ----------------------------------------------------------------------
def test_bond_file_formats():
    assert read_bond_file("# HOCl\n1 2\n1 3   # O-Cl\n", 3) == [(0, 1), (0, 2)]
    assert read_bond_file("1 2 1.0 3 1.0\n2\n3\n", 3) == [(0, 1), (0, 2)]
    assert read_bond_file("\r\n1 2\r\n", 3) == [(0, 1)]


def test_bond_file_errors():
    with pytest.raises(ParseError, match="references atom 4"):
        read_bond_file("1 4\n", 3)
    with pytest.raises(ParseError, match="no bonds"):
        read_bond_file("# nothing\n", 3)
    with pytest.raises(ParseError, match="neither an atom index"):
        read_bond_file("1 x\n", 3)
