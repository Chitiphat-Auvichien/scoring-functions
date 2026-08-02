"""Framework 2: parse the Cartesian Hessian directly out of a Gaussian
formatted checkpoint (.fchk) file, already in Hartree/Bohr^2 -- no
reconstruction from printed normal modes needed (contrast
hessian_reconstruct.py's Framework 1, which recovers the Hessian from
Gaussian's printed frequencies/eigenvectors instead).

WARNING: this module has NEVER been integration-tested against a real
Gaussian-produced .fchk file -- none exists anywhere in this repo (no
data/fchk/ directory, and no .log job here was ever run with %chk set).
Only a synthetic hand-built fixture (tests/test_veda_fmt_regression.py)
exercises the parsing logic below. Treat Framework 2 with suspicion on its
first real use, until it's been run against a real Gaussian .fchk and
independently cross-checked (e.g. against Framework 1's reconstruction for
the same job, or against VEDA4's own recomputed frequencies).

To obtain a usable .fchk: add `%chk=<name>.chk` to the Gaussian route of the
freq job, then, after the job finishes, run Gaussian's `formchk <name>.chk
<name>.fchk` utility to convert the binary checkpoint into the text .fchk
format this module parses. (A less portable alternative is requesting a
Hessian-printing route option directly in the .log itself, avoiding the
.fchk detour entirely -- but that output format isn't what this module
parses.)
"""
import os
import re
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)                        # sibling ped/ modules (blocks.py)
sys.path.insert(0, os.path.join(_HERE, '..'))     # repo root (src/)
from src.utils import get_atomic_number  # noqa: E402

import blocks  # noqa: E402

_N_EQUALS_RE = re.compile(r'N=\s*(\d+)')


def _parse_fchk_scalar(lines, header_prefix):
    """A .fchk scalar line looks like
    'Number of atoms                           I               12'
    -- prefix-match the header, take the last whitespace-separated token."""
    for line in lines:
        if line.startswith(header_prefix):
            parts = line.split()
            if not parts:
                raise ValueError(f"'{header_prefix}' scalar line has no tokens: {line!r}")
            return float(parts[-1])
    raise ValueError(f"'{header_prefix}' scalar not found in .fchk.")


def _parse_fchk_array(lines, header_prefix):
    """A .fchk array section looks like
    'Cartesian Force Constants                  R   N=          66'
    followed by N values, 5 per line, general/E-notation (defensively
    'D'->'E'-replaced in case of Fortran D-notation) -- greedily collect N
    whitespace-separated numeric tokens starting on the line after the
    header."""
    for i, line in enumerate(lines):
        if line.startswith(header_prefix):
            m = _N_EQUALS_RE.search(line)
            if not m:
                raise ValueError(f"'{header_prefix}' array header missing "
                                  f"'N=' count: {line!r}")
            n = int(m.group(1))
            values = []
            j = i + 1
            while len(values) < n and j < len(lines):
                for tok in lines[j].split():
                    values.append(float(tok.replace('D', 'E').replace('d', 'e')))
                    if len(values) == n:
                        break
                j += 1
            if len(values) != n:
                raise ValueError(f"'{header_prefix}' array expected {n} "
                                  f"values, found {len(values)}.")
            return values
    raise ValueError(f"'{header_prefix}' array not found in .fchk.")


def parse_fchk_hessian(fchk_path, natoms, symbols=None):
    """Locate 'Cartesian Force Constants' and reshape the flat lower
    triangle (row-major: row i has i+1 values, row1=[col1],
    row2=[col1,col2], ..., row-3N=[col1..col3N]) into a full symmetric
    (3*natoms, 3*natoms) array, already in Hartree/Bohr^2 -- no unit
    conversion needed.

    Cross-checks the .fchk's own 'Number of atoms' scalar against `natoms`,
    and (if `symbols` given) its 'Atomic numbers' array against
    get_atomic_number(symbols) -- raises ValueError loudly on any mismatch,
    catching the case where the .log and .fchk are from different jobs."""
    with open(fchk_path, 'r') as f:
        lines = f.readlines()

    n_atoms_fchk = int(round(_parse_fchk_scalar(lines, "Number of atoms")))
    if n_atoms_fchk != natoms:
        raise ValueError(
            f".fchk 'Number of atoms' ({n_atoms_fchk}) does not match the "
            f"expected natoms ({natoms}) -- the .log and .fchk are likely "
            f"from different jobs: {fchk_path}")

    if symbols is not None:
        atomic_numbers_fchk = [int(round(v)) for v in
                                _parse_fchk_array(lines, "Atomic numbers")]
        expected_numbers = [get_atomic_number(s) for s in symbols]
        if atomic_numbers_fchk != expected_numbers:
            raise ValueError(
                f".fchk 'Atomic numbers' {atomic_numbers_fchk} do not match "
                f"the expected {expected_numbers} derived from the .log's "
                f"symbols {symbols} -- the .log and .fchk are likely from "
                f"different jobs: {fchk_path}")

    flat = _parse_fchk_array(lines, "Cartesian Force Constants")
    n = 3 * natoms
    expected_len = n * (n + 1) // 2
    if len(flat) != expected_len:
        raise ValueError(
            f".fchk 'Cartesian Force Constants' has {len(flat)} values, "
            f"expected {expected_len} for a {n}x{n} lower triangle "
            f"(natoms={natoms}).")

    H = np.zeros((n, n))
    idx = 0
    for i in range(n):
        for j in range(i + 1):
            H[i, j] = flat[idx]
            H[j, i] = flat[idx]
            idx += 1
    return H


def build_hessian(mol, symbols, vibfreq, masses):
    """Orchestrate Framework 2: parse .fchk -> validate against the .log's
    own frequencies (the only sanity check available until real .fchk data
    exists). Returns (H_AU, max_diff_cm1)."""
    if not mol.fchk_path:
        raise FileNotFoundError(
            f"No .fchk resolved for molecule {mol.name!r}; Framework 2 "
            "(fchk) requires one -- see this module's docstring for how to "
            "produce one.")
    print(f"Parsing Hessian for {mol.name} from {mol.fchk_path} "
          "(Framework 2: parse from .fchk)")
    H_AU = parse_fchk_hessian(mol.fchk_path, len(symbols), symbols=symbols)
    max_diff = blocks.round_trip_validate_au(H_AU, masses, vibfreq)
    return H_AU, max_diff
