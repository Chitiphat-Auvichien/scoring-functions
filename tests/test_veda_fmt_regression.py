"""Regression tests for ped/'s VEDA4-.fmt-building tooling (2026-08-02
reorg: ped/12_export_veda_fmt.py -> ped/build_veda_fmt.py + ped/blocks.py +
ped/molecule.py + ped/hessian_reconstruct.py + ped/hessian_fchk.py).

Run from ``Github/scoring-functions/``:
    py -m pytest tests/test_veda_fmt_regression.py
    py tests/test_veda_fmt_regression.py       (standalone; no pytest needed)

Covers:
  - blocks.fortran_d against byte-exact values pulled directly from
    ped/reference/C6H6.fmt (the VEDA4-validated baseline -- these are the
    same numbers ped/archive_python_ped/12_export_veda_fmt.py's docstring
    describes as byte-exact-verified against A1-91.FMT).
  - molecule.resolve_molecule's dual-extension (.com/.gjf) gjf lookup.
  - Framework 1 end-to-end: round-trip discrepancy < 2.0 cm^-1, and the
    freshly built C6H6 .fmt's Hessian floats match ped/reference/C6H6.fmt
    numerically (nothing about the math changed in the reorg, only the
    code organization).
  - hessian_fchk's scalar/array/.fchk-Hessian parsing against a synthetic,
    hand-built fixture only -- there is no real Gaussian .fchk anywhere in
    this repo to validate against (see hessian_fchk.py's module docstring).
"""
import os
import re
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "ped"))

import blocks  # noqa: E402
import hessian_fchk  # noqa: E402
import hessian_reconstruct  # noqa: E402
from molecule import resolve_molecule  # noqa: E402
from src.utils import atomicMass  # noqa: E402

REFERENCE_FMT = os.path.join(ROOT, "ped", "reference", "C6H6.fmt")


def _extract_hessian_floats(fmt_path):
    """Pull every Fortran D-notation float out of the 'Force constants in
    Cartesian coordinates:' block of a .fmt file, in file order."""
    with open(fmt_path) as f:
        text = f.read()
    idx = text.index("Force constants in Cartesian coordinates:")
    tail = text[idx:]
    matches = re.findall(r"[-]?\d\.\d{6}D[+-]\d{2}", tail)
    return np.array([float(m.replace("D", "E")) for m in matches])


# ---------------------------------------------------------------------------
# blocks.fortran_d -- byte-exact against ped/reference/C6H6.fmt
# ---------------------------------------------------------------------------

def test_fortran_d_matches_reference_fmt_values():
    """Pull a few real (value, formatted-string) pairs straight out of the
    already-VEDA4-validated ped/reference/C6H6.fmt and confirm
    blocks.fortran_d reproduces the exact same string -- these are the same
    reference values ped/archive_python_ped/12_export_veda_fmt.py's
    docstring already verified byte-exact against A1-91.FMT."""
    with open(REFERENCE_FMT) as f:
        text = f.read()
    idx = text.index("Force constants in Cartesian coordinates:")
    tail = text[idx:]
    tokens = re.findall(r"[-]?\d\.\d{6}D[+-]\d{2}", tail)
    assert len(tokens) > 100  # sanity: the Hessian block was actually found

    # Deterministic sample: first, middle, and last formatted value.
    sample_idxs = [0, len(tokens) // 2, len(tokens) - 1]
    for i in sample_idxs:
        token = tokens[i]
        value = float(token.replace("D", "E"))
        formatted = blocks.fortran_d(value).strip()
        assert formatted == token, f"index {i}: {formatted!r} != {token!r}"


def test_fortran_d_zero():
    assert blocks.fortran_d(0.0).strip() == "0.000000D+00"


def test_fortran_d_rounding_edge_case():
    """A mantissa that rounds up to 1.0 must re-normalize into the next
    decade (e.g. 0.9999995 -> 1.000000D+00 would be wrong; must renormalize
    to 0.100000D+01)."""
    formatted = blocks.fortran_d(0.9999999).strip()
    mantissa_str, exp_str = formatted.split("D")
    mantissa = float(mantissa_str)
    assert 0.1 <= mantissa < 1.0, f"mantissa {mantissa} not in [0.1, 1.0)"
    # Round-trip: reparse and confirm it's numerically close to the input.
    assert abs(float(formatted.replace("D", "E")) - 0.9999999) < 1e-6


# ---------------------------------------------------------------------------
# molecule.resolve_molecule -- dual-extension (.com/.gjf) gjf lookup
# ---------------------------------------------------------------------------

def test_resolve_molecule_com_extension():
    mol = resolve_molecule("H2O", ROOT)
    assert mol.log_path == os.path.join(ROOT, "data", "logs", "H2O.log")
    assert os.path.isfile(mol.log_path)
    assert mol.gjf_path == os.path.join(ROOT, "data", "gjf", "H2O.com")
    assert os.path.isfile(mol.gjf_path)


def test_resolve_molecule_gjf_extension():
    mol = resolve_molecule("PH3", ROOT)
    assert mol.log_path == os.path.join(ROOT, "data", "logs", "PH3.log")
    assert os.path.isfile(mol.log_path)
    assert mol.gjf_path == os.path.join(ROOT, "data", "gjf", "PH3.gjf")
    assert os.path.isfile(mol.gjf_path)


def test_resolve_molecule_missing_log_raises():
    raised = False
    try:
        resolve_molecule("NoSuchMolecule123", ROOT)
    except FileNotFoundError:
        raised = True
    assert raised


def test_resolve_molecule_no_fchk_present():
    """No .fchk exists anywhere in this repo (see hessian_fchk.py's
    docstring) -- fchk_path must resolve to None for a real molecule."""
    mol = resolve_molecule("C6H6", ROOT)
    assert mol.fchk_path is None


# ---------------------------------------------------------------------------
# Framework 1 (reconstruct) end-to-end regression, C6H6
# ---------------------------------------------------------------------------

def test_framework1_c6h6_round_trip_and_matches_reference():
    mol = resolve_molecule("C6H6", ROOT)
    data = hessian_reconstruct.load_geometry_and_modes(mol.log_path)
    symbols = data["atoms"]
    coords_ang = data["coords"]
    modes = data["modes"]

    vibfreq = np.array([m["frequency"] for m in modes])
    reduced_mass = np.array([m["reduced_mass"] for m in modes])
    l_raw = np.array([m["vector"] for m in modes])
    masses = np.array([atomicMass[s.lower()] for s in symbols])

    H_AU, max_diff = hessian_reconstruct.build_hessian(
        mol, symbols, coords_ang, vibfreq, masses, reduced_mass, l_raw)

    assert max_diff < 2.0, f"round-trip discrepancy {max_diff} >= 2.0 cm^-1"

    with open(mol.log_path) as f:
        log_lines = f.readlines()
    geom_block = blocks.extract_geometry_block(log_lines, coords_ang)
    freq_block = blocks.extract_standard_freq_block(log_lines, len(symbols), len(vibfreq))
    hessian_text = blocks.format_cartesian_hessian_block(H_AU)
    fresh_fmt_text = blocks.assemble_fmt(geom_block, freq_block, hessian_text)

    fresh_floats = np.array(
        [float(m.replace("D", "E")) for m in
         re.findall(r"[-]?\d\.\d{6}D[+-]\d{2}",
                     fresh_fmt_text[fresh_fmt_text.index(
                         "Force constants in Cartesian coordinates:"):])])
    reference_floats = _extract_hessian_floats(REFERENCE_FMT)

    assert len(fresh_floats) == len(reference_floats)
    assert np.allclose(fresh_floats, reference_floats, atol=1e-6), (
        "Freshly reconstructed Hessian floats diverge from "
        "ped/reference/C6H6.fmt -- max abs diff = "
        f"{np.max(np.abs(fresh_floats - reference_floats))}")


# ---------------------------------------------------------------------------
# hessian_fchk -- synthetic fixture only (no real .fchk exists in this repo)
# ---------------------------------------------------------------------------

def _synthetic_fchk_lines():
    """A hand-built, minimal fake .fchk excerpt for a fictitious 2-atom
    (H, O) 'molecule' -- exercises _parse_fchk_scalar/_parse_fchk_array/
    parse_fchk_hessian's shape logic only; the numbers are arbitrary."""
    lines = [
        "Fake Gaussian fchk fixture\n",
        "SP        RHF                                                     STO-3G\n",
        "Number of atoms                           I               2\n",
        "Atomic numbers                            I   N=           2\n",
        "           1           8\n",
        "Real atomic weights                       R   N=           2\n",
        "  1.00782504E+00  1.59949146E+01\n",
        "Cartesian Force Constants                 R   N=          21\n",
    ]
    flat = [float(v) for v in range(1, 22)]  # 1.0 .. 21.0
    for i in range(0, 21, 5):
        chunk = flat[i:i + 5]
        lines.append("".join(f"  {v:.6f}D+00" for v in chunk) + "\n")
    return lines


def test_parse_fchk_scalar():
    lines = _synthetic_fchk_lines()
    assert hessian_fchk._parse_fchk_scalar(lines, "Number of atoms") == 2.0


def test_parse_fchk_scalar_not_found_raises():
    lines = _synthetic_fchk_lines()
    raised = False
    try:
        hessian_fchk._parse_fchk_scalar(lines, "Total Energy")
    except ValueError:
        raised = True
    assert raised


def test_parse_fchk_array():
    lines = _synthetic_fchk_lines()
    values = hessian_fchk._parse_fchk_array(lines, "Atomic numbers")
    assert values == [1.0, 8.0]

    hessian_vals = hessian_fchk._parse_fchk_array(lines, "Cartesian Force Constants")
    assert len(hessian_vals) == 21
    assert hessian_vals == [float(v) for v in range(1, 22)]


def test_parse_fchk_hessian_shape_and_symmetry():
    lines = _synthetic_fchk_lines()
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        fchk_path = os.path.join(d, "fake.fchk")
        with open(fchk_path, "w") as f:
            f.writelines(lines)

        H = hessian_fchk.parse_fchk_hessian(fchk_path, natoms=2, symbols=["H", "O"])
        assert H.shape == (6, 6)
        assert np.allclose(H, H.T)
        # Lower-triangle row-major fill: row0=[1], row1=[2,3], row2=[4,5,6]...
        assert H[0, 0] == 1.0
        assert H[1, 0] == H[0, 1] == 2.0
        assert H[1, 1] == 3.0
        assert H[5, 5] == 21.0


def test_parse_fchk_hessian_natoms_mismatch_raises():
    lines = _synthetic_fchk_lines()
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        fchk_path = os.path.join(d, "fake.fchk")
        with open(fchk_path, "w") as f:
            f.writelines(lines)

        raised = False
        try:
            hessian_fchk.parse_fchk_hessian(fchk_path, natoms=3)
        except ValueError:
            raised = True
        assert raised


def test_parse_fchk_hessian_symbol_mismatch_raises():
    lines = _synthetic_fchk_lines()
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        fchk_path = os.path.join(d, "fake.fchk")
        with open(fchk_path, "w") as f:
            f.writelines(lines)

        raised = False
        try:
            hessian_fchk.parse_fchk_hessian(fchk_path, natoms=2, symbols=["H", "H"])
        except ValueError:
            raised = True
        assert raised


def test_parse_fchk_masses():
    lines = _synthetic_fchk_lines()
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        fchk_path = os.path.join(d, "fake.fchk")
        with open(fchk_path, "w") as f:
            f.writelines(lines)

        masses = hessian_fchk.parse_fchk_masses(fchk_path, natoms=2)
        assert np.allclose(masses, [1.00782504, 15.9949146])


def test_parse_fchk_masses_natoms_mismatch_raises():
    lines = _synthetic_fchk_lines()
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        fchk_path = os.path.join(d, "fake.fchk")
        with open(fchk_path, "w") as f:
            f.writelines(lines)

        raised = False
        try:
            hessian_fchk.parse_fchk_masses(fchk_path, natoms=3)
        except ValueError:
            raised = True
        assert raised


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    passed = 0
    for fn in tests:
        try:
            fn()
            print(f"PASS  {fn.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"FAIL  {fn.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            print(f"ERROR {fn.__name__}: {type(e).__name__}: {e}")
    print(f"\n{passed}/{len(tests)} passed")
    sys.exit(0 if passed == len(tests) else 1)
