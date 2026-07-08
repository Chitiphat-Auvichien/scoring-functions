"""Regression tests for src/parser.py's Gaussian-direct rework (2026-07-08):
mu (reduced mass) / k (force constant) / irrep extraction in GaussianParser,
and the new Gaussian-direct IntermediateIO normal-mode format.

Run from ``Github/scoring-functions/``:
    py -m pytest tests/test_parser.py
    py tests/test_parser.py       (standalone; no pytest needed)

Uses real on-disk .log files already in data/logs/ (no synthetic fixtures),
matching the convention the rest of tests/ uses:
  - data/logs/C6H6_MP2_3-21G_D6h.log -- clean HP-block case (real irreps).
  - data/logs/AlCl3.log -- HP block with placeholder irreps ('?A', '?B',
    and one real symbol 'A2"') -- must be stored as-is, not validated.
  - data/logs/TeH2.log -- has BOTH an HP block and a (redundant) Standard
    block back-to-back for the same 3 modes; GaussianParser's existing
    force_hp_only logic must still only keep the HP block's 3 modes (not
    6), and mu/k/irrep extraction must work correctly for the HP path here.
  - A direct, isolated Standard-block extraction check (bypassing
    force_hp_only) confirms the Standard branch's mu/k/irrep offsets too.
"""
import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from src.parser import GaussianParser, IntermediateIO   # noqa: E402

TOL = 1e-3


def _parse(mol):
    gp = GaussianParser(os.path.join(ROOT, "data", "logs", f"{mol}.log"))
    return gp.parse(parse_modes=True)


# ---------------------------------------------------------------------------
# HP block: mu/k/irrep extraction
# ---------------------------------------------------------------------------

def test_benzene_hp_block_mu_k_irrep():
    """Clean HP-block case: first 6 modes' mu/k/irrep match the raw .log
    (data/logs/C6H6_MP2_3-21G_D6h.log:1789-1793) exactly."""
    raw = _parse("C6H6_MP2_3-21G_D6h")
    modes = raw["modes"]
    assert len(modes) == 30  # 3*12-6

    expected = [
        (411.9847, 2.9803, 0.2980, "E2U"),
        (411.9847, 2.9803, 0.2980, "E2U"),
        (645.8926, 6.2210, 1.5291, "E2G"),
        (645.8926, 6.2210, 1.5291, "E2G"),
        (678.1652, 7.5391, 2.0429, "B2G"),
    ]
    for m, (freq, mu, k, irrep) in zip(modes[:5], expected):
        assert abs(m["frequency"] - freq) < TOL
        assert abs(m["reduced_mass"] - mu) < TOL
        assert abs(m["force_constant"] - k) < TOL
        assert m["irrep"] == irrep

    # Every mode in this clean file must have real (non-None) mu/k/irrep.
    for m in modes:
        assert m["reduced_mass"] is not None
        assert m["force_constant"] is not None
        assert m["irrep"] is not None


def test_alcl3_placeholder_irrep_stored_as_is():
    """AlCl3.log:1164 has placeholder irreps ('?A', '?B') alongside one real
    symbol ('A2"') in the SAME row -- all must be stored verbatim, never
    validated against a point-group table."""
    raw = _parse("AlCl3")
    modes = raw["modes"]
    assert len(modes) == 6  # 3*4-6

    expected_irreps = ["?A", "?A", 'A2"', "?B", "?A", "?A"]
    assert [m["irrep"] for m in modes] == expected_irreps

    expected = [
        (139.7371, 33.7431, 0.3882),
        (139.7371, 33.7431, 0.3882),
        (189.5614, 28.3041, 0.5992),
        (351.2037, 34.9689, 2.5413),
        (586.2406, 29.1616, 5.9049),
        (586.2406, 29.1616, 5.9049),
    ]
    for m, (freq, mu, k) in zip(modes, expected):
        assert abs(m["frequency"] - freq) < TOL
        assert abs(m["reduced_mass"] - mu) < TOL
        assert abs(m["force_constant"] - k) < TOL


def test_teh2_hp_and_standard_back_to_back_keeps_only_hp():
    """TeH2.log has both an HP block (line ~1221) and a redundant Standard
    block (line ~1241) for the SAME 3 modes -- force_hp_only must still keep
    only the HP block's 3 modes (not 6), with correct mu/k/irrep."""
    raw = _parse("TeH2")
    modes = raw["modes"]
    assert len(modes) == 3  # 3*3-6

    expected = [
        (892.9572, 1.0160, 0.4773, "A1"),
        (2062.3274, 1.0150, 2.5436, "A1"),
        (2070.5516, 1.0157, 2.5657, "B2"),
    ]
    for m, (freq, mu, k, irrep) in zip(modes, expected):
        assert abs(m["frequency"] - freq) < TOL
        assert abs(m["reduced_mass"] - mu) < TOL
        assert abs(m["force_constant"] - k) < TOL
        assert m["irrep"] == irrep


def test_teh2_standard_block_mu_k_irrep_direct():
    """Isolate the Standard (2-dash) block branch directly: force_hp_only=False
    on TeH2's Standard 'Frequencies --' line (data/logs/TeH2.log:1243) must
    extract the SAME mu/k/irrep as the HP block above -- confirming the
    Standard branch's offsets (line_idx-1/+1/+2) independently of the HP path."""
    gp = GaussianParser(os.path.join(ROOT, "data", "logs", "TeH2.log"))
    gp._parse_standard_orientation()
    gp.natoms = len(gp.coordinates)

    standard_line_idx = None
    for i, line in enumerate(gp.lines):
        # The Standard block's own "Frequencies --" line (2-dash, no 3rd
        # dash) -- distinguish from the HP block's "Frequencies ---" by
        # requiring the line NOT contain the 3-dash form.
        if "Frequencies --" in line and "Frequencies ---" not in line:
            standard_line_idx = i
            break
    assert standard_line_idx is not None, "no Standard-style Frequencies line found in TeH2.log"

    gp.modes = []
    gp._parse_block(standard_line_idx, force_hp_only=False)
    modes = gp.modes
    assert len(modes) == 3

    expected = [
        (892.9572, 1.0160, 0.4773, "A1"),
        (2062.3274, 1.0150, 2.5436, "A1"),
        (2070.5516, 1.0157, 2.5657, "B2"),
    ]
    for m, (freq, mu, k, irrep) in zip(modes, expected):
        assert abs(m["frequency"] - freq) < TOL
        assert abs(m["reduced_mass"] - mu) < TOL
        assert abs(m["force_constant"] - k) < TOL
        assert m["irrep"] == irrep


# ---------------------------------------------------------------------------
# IntermediateIO -- new Gaussian-direct normal-mode format round-trip
# ---------------------------------------------------------------------------

def _roundtrip(mol, tmp_path):
    raw = _parse(mol)
    IntermediateIO.save(raw, tmp_path)
    loaded = IntermediateIO.load(tmp_path)
    return raw, loaded


def test_intermediate_roundtrip_water(tmp_path):
    tmp_file = os.path.join(str(tmp_path), "H2O-MP2-321G_normal_data.txt")
    raw, loaded = _roundtrip("H2O-MP2-321G", tmp_file)

    assert loaded["atoms"] == raw["atoms"]
    assert sorted(loaded["bonds"]) == sorted(raw["bonds"])
    assert np.allclose(loaded["coords"], raw["coords"], atol=1e-5)
    assert len(loaded["modes"]) == len(raw["modes"])

    for m1, m2 in zip(raw["modes"], loaded["modes"]):
        assert abs(m1["frequency"] - m2["frequency"]) < 1e-3
        assert abs(m1["reduced_mass"] - m2["reduced_mass"]) < 1e-3
        assert abs(m1["force_constant"] - m2["force_constant"]) < 1e-3
        assert m1["irrep"] == m2["irrep"]
        assert np.allclose(m1["vector"], m2["vector"], atol=1e-5)
        assert m2["is_emit"] is False


def test_intermediate_roundtrip_benzene_multi_block(tmp_path):
    """30 modes -> 6 frequency blocks of 5 -- exercises the block-boundary
    logic (5-modes-per-block) end to end, not just a single block."""
    tmp_file = os.path.join(str(tmp_path), "C6H6_MP2_3-21G_D6h_normal_data.txt")
    raw, loaded = _roundtrip("C6H6_MP2_3-21G_D6h", tmp_file)

    assert len(loaded["modes"]) == 30
    for m1, m2 in zip(raw["modes"], loaded["modes"]):
        assert abs(m1["frequency"] - m2["frequency"]) < 1e-3
        assert abs(m1["reduced_mass"] - m2["reduced_mass"]) < 1e-3
        assert abs(m1["force_constant"] - m2["force_constant"]) < 1e-3
        assert m1["irrep"] == m2["irrep"]
        assert np.allclose(m1["vector"], m2["vector"], atol=1e-5)


def test_intermediate_header_line_matches_natoms_linear_flag(tmp_path):
    tmp_file = os.path.join(str(tmp_path), "H2O-MP2-321G_normal_data.txt")
    raw, _ = _roundtrip("H2O-MP2-321G", tmp_file)
    with open(tmp_file) as f:
        header = f.readline().split()
    assert int(header[0]) == 3          # natoms
    assert int(header[1]) == 0          # water is nonlinear
    assert header[2] == "H2O-MP2-321G"  # molecule name derived from filename


def test_intermediate_load_fails_loud_on_tampered_linear_flag(tmp_path):
    """load() must fail loud if the header's LINEAR flag disagrees with the
    actual mode count -- do not silently trust a stale/hand-edited flag."""
    tmp_file = os.path.join(str(tmp_path), "H2O-MP2-321G_normal_data.txt")
    raw = _parse("H2O-MP2-321G")
    IntermediateIO.save(raw, tmp_file)

    with open(tmp_file) as f:
        lines = f.readlines()
    parts = lines[0].split()
    parts[1] = "1"  # water is nonlinear (0) -- flip to a wrong value
    lines[0] = " ".join(parts) + "\n"
    with open(tmp_file, "w") as f:
        f.writelines(lines)

    raised = False
    try:
        IntermediateIO.load(tmp_file)
    except ValueError:
        raised = True
    assert raised, "load() should raise ValueError on a tampered LINEAR flag"


def test_intermediate_load_accepts_bare_neighbor_pair(tmp_path):
    """A hand-edited bare 'atom neighbor' bond line (CLAUDE.md's documented
    manual bond-repair workflow, e.g. '1 2') must be accepted, defaulting
    bond order to 1.0, alongside the full 'atom neighbor order' form."""
    tmp_file = os.path.join(str(tmp_path), "H2O-MP2-321G_normal_data.txt")
    raw = _parse("H2O-MP2-321G")
    IntermediateIO.save(raw, tmp_file)

    with open(tmp_file) as f:
        lines = f.readlines()
    # Line 1 (0-based index 1, right after the header) is atom 1's bond line
    # -- currently the full form (e.g. ' 1 2 1.0 3 1.0'). Rewrite it as a
    # bare pair to just the first bond, dropping the order token and the
    # second neighbor, to simulate a human's simplified manual edit.
    assert lines[1].strip().startswith("1 ")
    lines[1] = " 1 2\n"  # bare pair: atom 1 bonded to atom 2 only, no order
    with open(tmp_file, "w") as f:
        f.writelines(lines)

    loaded = IntermediateIO.load(tmp_file)
    assert (0, 1) in loaded["bonds"]


if __name__ == "__main__":
    import tempfile
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    passed = 0
    for fn in tests:
        try:
            if "tmp_path" in fn.__code__.co_varnames[: fn.__code__.co_argcount]:
                with tempfile.TemporaryDirectory() as d:
                    fn(d)
            else:
                fn()
            print(f"PASS  {fn.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"FAIL  {fn.__name__}: {e}")
        except Exception as e:  # noqa: BLE001
            print(f"ERROR {fn.__name__}: {type(e).__name__}: {e}")
    print(f"\n{passed}/{len(tests)} passed")
    sys.exit(0 if passed == len(tests) else 1)
