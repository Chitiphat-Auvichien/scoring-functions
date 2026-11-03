# [EDIT 2026-11-02 | by Claude (AI assistant) for Than (ThanchonBK)]
# Purpose: new module -- ORCA .hess / Molden readers, so the pipeline accepts
#          non-Gaussian programs (JCC reviewer 1). Ported from the vscore web app.
# Previous version: none -- this whole file is new.
# [END EDIT NOTE]
"""ORCA .hess and Molden readers -- the non-Gaussian inputs to the pipeline.

Each returns ``{'atoms', 'coords', 'modes'}`` in exactly the shape
GaussianParser produces: element symbols, coordinates in Angstrom, and 3N-6
(3N-5) modes whose vectors are unit-norm Cartesian displacements. Bonds are
not read here -- neither format carries them -- see ``read_bond_file``.

Ported from the web app's ``app/core/parsers.py`` (vscore-webapp), which
carries the validation record: checked against real ORCA output by
rediagonalising each file's own Hessian, against Gaussian 16 for the same
molecule and level of theory, and fuzzed with 20 000 corrupted files.

    * ``parse_orca_hess``  ``$atoms`` / ``$vibrational_frequencies`` /
      ``$normal_modes`` of an ORCA ``.hess`` (bohr -> Angstrom)
    * ``parse_molden``     ``[FREQ]`` / ``[FR-COORD]`` / ``[FR-NORM-COORD]`` of a
      Molden file -- Molpro ``put,molden`` and any other writer
    * the zero-frequency translations/rotations a 3N listing carries are
      dropped, and mass-weighted vectors are detected and converted
      (``_cartesian_modes``), so the scorer never sees a different convention.

Errors are ``ParseError``, a ``ValueError`` subclass, matching the rest of
``src/parser.py``, which raises ValueError for malformed input.
"""

import re

import numpy as np

from .utils import atomicMass, get_symbol


class ParseError(ValueError):
    """Malformed or unsupported ORCA/Molden input; the message says why."""


def _is_int(s):
    return s.lstrip("+-").isdigit()


# Neither format carries connectivity, so -- exactly as for a Gaussian .log --
# the bonds come from a separate file (``read_bond_file``) and are never
# inferred from distances.
BOHR_TO_ANGSTROM = 0.529177210903

# A non-orthogonality above this, under the best-fitting metric, means the
# vectors are not normal modes of one Hessian (or were mis-parsed). Gaussian's
# 5-dp hpmodes already reach 7.5e-3 on heavy-atom molecules (PCl3, SeBr2), and
# 4-dp writers more; unrelated vectors overlap at O(0.1-1).
_ORTHO_TOL = 5e-2
# How much better the mass-weighted reading must fit before vectors are
# converted; anything less is ambiguous and stays Cartesian.
_DECISIVE = 3.0


def parse_orca_hess(text):
    """Geometry and normal modes from an ORCA ``.hess`` file.

    Reads ``$atoms`` (symbol, mass, x y z in bohr), ``$vibrational_frequencies``
    and ``$normal_modes`` (a 3N x 3N matrix, one mode per column, printed in
    column blocks). ORCA documents these columns as Cartesian displacements,
    normalised -- the same convention as Gaussian -- which
    ``_cartesian_modes`` confirms from the vectors themselves rather than
    trusting. The Hessian block is not needed and is not parsed.

    Returns ``{'atoms', 'coords', 'modes'}``; no 'bonds', see above.
    """
    sections = _dollar_sections(text)
    for req in ("atoms", "vibrational_frequencies", "normal_modes"):
        if req not in sections:
            raise ParseError(
                f"No ${req} block found. This does not look like an ORCA .hess file "
                "from a frequency calculation.")

    rows = sections["atoms"]
    try:
        natoms = int(rows[0].split()[0])
    except (IndexError, ValueError):
        raise ParseError("The $atoms block does not start with an atom count.")
    if natoms < 1:
        raise ParseError("The $atoms block declares no atoms.")
    atoms, coords, masses = [], [], []
    for n, line in enumerate(rows[1:natoms + 1], 1):
        parts = line.split()
        if len(parts) < 5:
            raise ParseError(f"$atoms line {n} should read 'symbol mass x y z': {line!r}")
        atoms.append(_element(parts[0], f"$atoms line {n}"))
        mass, *xyz = _numbers(parts[1:5], f"$atoms line {n}")
        if mass <= 0:
            raise ParseError(f"$atoms line {n}: mass {mass} is not positive.")
        masses.append(mass)
        coords.append([v * BOHR_TO_ANGSTROM for v in xyz])
    if len(atoms) != natoms:
        raise ParseError(f"$atoms declares {natoms} atoms but lists {len(atoms)}.")

    freqs = []
    for line in sections["vibrational_frequencies"][1:]:
        parts = line.split()
        if len(parts) >= 2:
            freqs.append(_numbers(parts[1:2], "$vibrational_frequencies")[0])
    ncoord = 3 * natoms
    if len(freqs) != ncoord:
        raise ParseError(
            f"$vibrational_frequencies lists {len(freqs)} values, expected 3N = {ncoord}.")

    matrix = _orca_column_blocks(sections["normal_modes"], ncoord)
    modes = [{"frequency": freqs[k], "vector": matrix[:, k].reshape(natoms, 3),
              "is_emit": False, "reduced_mass": None, "force_constant": None,
              "irrep": None} for k in range(ncoord)]
    modes = _drop_projected_modes(modes, natoms, np.array(coords))
    _cartesian_modes(modes, np.array(masses), "ORCA .hess")
    return {"atoms": atoms, "coords": np.array(coords), "modes": modes}


def _dollar_sections(text):
    """'$name' blocks of an ORCA .hess -> {name: [non-blank lines]}.

    ORCA puts '#' comment lines between and inside blocks (e.g. '# The atoms:
    label  mass x y z'); they are dropped.
    """
    sections, current = {}, None
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("#"):
            continue
        if s.startswith("$"):
            current = s[1:].split()[0].lower() if len(s) > 1 else None
            if current == "end":
                break
            sections[current] = []
        elif current and s:
            sections[current].append(s)
    return sections


def _orca_column_blocks(lines, ncoord):
    """The '$normal_modes' matrix: a 'rows cols' line, then blocks of columns,
    each a header of column indices followed by 'row v v v ...' lines."""
    try:
        nrow, ncol = (int(v) for v in lines[0].split()[:2])
    except (IndexError, ValueError):
        raise ParseError("$normal_modes does not start with its 'rows cols' dimensions.")
    if nrow != ncoord or ncol != ncoord:
        raise ParseError(
            f"$normal_modes is {nrow} x {ncol}, but the geometry needs {ncoord} x {ncoord}.")

    matrix = np.full((ncoord, ncoord), np.nan)
    cols = []
    for line in lines[1:]:
        parts = line.split()
        if all(_is_int(p) for p in parts):          # a column-index header
            cols = [int(p) for p in parts]
            if len(set(cols)) != len(cols) or not all(0 <= c < ncoord for c in cols):
                raise ParseError(f"$normal_modes has an invalid column header: {line!r}")
            continue
        if not _is_int(parts[0]):
            raise ParseError(f"$normal_modes has an unreadable line: {line!r}")
        row = int(parts[0])
        vals = _numbers(parts[1:], "$normal_modes")
        if len(vals) != len(cols) or not 0 <= row < ncoord:
            raise ParseError(f"$normal_modes line does not match its block header: {line!r}")
        matrix[row, cols] = vals
    if np.isnan(matrix).any():
        raise ParseError("$normal_modes is incomplete -- the file looks truncated.")
    return matrix


_MOLDEN_SECTION_RE = re.compile(r"^\s*\[([A-Za-z0-9_\- ]+)\]\s*(.*)$")


def parse_molden(text):
    """Geometry and normal modes from the vibrational part of a Molden file.

    Reads ``[FREQ]``, ``[FR-COORD]`` (symbol x y z, bohr) and
    ``[FR-NORM-COORD]`` (one 'vibration n' header per mode, then N rows of
    dx dy dz). This covers Molpro (``put,molden``) and any other program that
    writes the format. Writers disagree on whether the vectors are mass-weighted,
    and the format does not say, so the convention is detected from the vectors
    (``_cartesian_modes``) and converted to Cartesian if needed. The [Atoms],
    [GTO] and [MO] sections are not read: [FR-COORD] is the geometry the modes
    belong to.

    Returns ``{'atoms', 'coords', 'modes'}``; no 'bonds'.
    """
    sections, current = {}, None
    for line in text.splitlines():
        m = _MOLDEN_SECTION_RE.match(line)
        if m:
            current = m.group(1).strip().upper()
            sections[current] = []
            continue
        if current and line.strip():
            sections[current].append(line.strip())

    if not text.lstrip().upper().startswith("[MOLDEN FORMAT]"):
        raise ParseError("This does not look like a Molden file (no '[Molden Format]' header).")
    for req in ("FREQ", "FR-COORD", "FR-NORM-COORD"):
        if req not in sections:
            raise ParseError(
                f"No [{req}] section. The Molden file has no vibrational data -- write "
                "it after a frequency calculation (Molpro: 'put,molden' after 'frequencies').")

    atoms, coords = [], []
    for n, line in enumerate(sections["FR-COORD"], 1):
        parts = line.split()
        if len(parts) < 4:
            raise ParseError(f"[FR-COORD] line {n} should read 'symbol x y z': {line!r}")
        atoms.append(_element(parts[0], f"[FR-COORD] line {n}"))
        coords.append([v * BOHR_TO_ANGSTROM
                       for v in _numbers(parts[1:4], f"[FR-COORD] line {n}")])
    natoms = len(atoms)
    if natoms == 0:
        raise ParseError("[FR-COORD] is empty.")

    freqs = [_numbers(line.split()[:1], "[FREQ]")[0] for line in sections["FREQ"]]

    vectors, rows = [], None
    for line in sections["FR-NORM-COORD"]:
        if line.lower().startswith("vibration"):
            if rows is not None:
                vectors.append(rows)
            rows = []
            continue
        if rows is None:
            raise ParseError(f"[FR-NORM-COORD] row before any 'vibration' header: {line!r}")
        parts = line.split()
        if len(parts) < 3:
            raise ParseError(f"[FR-NORM-COORD] row is not three numbers: {line!r}")
        rows.append(_numbers(parts[:3], "[FR-NORM-COORD]"))
    if rows is not None:
        vectors.append(rows)

    if len(vectors) != len(freqs):
        raise ParseError(
            f"[FREQ] lists {len(freqs)} frequencies but [FR-NORM-COORD] holds "
            f"{len(vectors)} vibrations.")
    bad = [k + 1 for k, v in enumerate(vectors) if len(v) != natoms]
    if bad:
        raise ParseError(
            f"Vibration {bad[0]} has the wrong number of rows for {natoms} atoms in "
            "[FR-COORD].")

    modes = [{"frequency": f, "vector": np.array(v), "is_emit": False,
              "reduced_mass": None, "force_constant": None, "irrep": None}
             for f, v in zip(freqs, vectors)]
    modes = _drop_projected_modes(modes, natoms, np.array(coords))
    masses = np.array([atomicMass[a.lower()] for a in atoms])
    _cartesian_modes(modes, masses, "Molden file")
    return {"atoms": atoms, "coords": np.array(coords), "modes": modes}


def _numbers(tokens, where):
    """Finite floats, or a ParseError naming the offending token.

    NaN and inf are refused, not just non-numbers: a NaN coordinate makes the
    linearity SVD fail to converge, and a NaN frequency cannot be serialised to
    the JSON the result page is built from -- both would surface as a 500.
    """
    out = []
    for tok in tokens:
        try:
            v = float(tok)
        except ValueError:
            raise ParseError(f"{where}: {tok!r} is not a number.")
        if not np.isfinite(v):
            raise ParseError(f"{where}: {tok!r} is not a finite number.")
        out.append(v)
    return out


def _element(tok, where):
    """Element symbol from a token like 'C', 'cl', 'C1' or '6'."""
    if tok.isdigit():
        sym = get_symbol(int(tok))
        if sym != "X":
            return sym
    else:
        sym = re.sub(r"[^A-Za-z]", "", tok).capitalize()
        if sym.lower() in atomicMass:
            return sym
    raise ParseError(
        f"{where}: unrecognised element {tok!r}; no atomic mass on file. Element "
        "identity affects all seven scores, so this cannot be defaulted.")


def _is_linear(coords, tol=1e-3):
    centred = coords - coords.mean(axis=0)
    if len(coords) < 3:
        return True
    return np.linalg.svd(centred, compute_uv=False)[1] < tol


def _drop_projected_modes(modes, natoms, coords):
    """Remove the projected-out translations/rotations a 3N listing carries.

    ORCA (and some Molden writers) list all 3N modes, with the 5-6 external
    ones set to exactly 0 cm-1 -- as all-zero vectors in a .hess, as real T/R
    vectors in symmetry-adapted runs. Those are dropped so the result has the
    same 3N-6 / 3N-5 vibrations a Gaussian log does. A set with no zeros is
    passed through untouched (a 3N set with unprojected T/R, or vibrations
    only), and the pipeline's own mode-set detection decides what it is.
    """
    zero = [k for k, m in enumerate(modes)
            if abs(m["frequency"]) < 1e-3 or not np.any(m["vector"])]
    if not zero:
        return modes
    n_ext = 5 if _is_linear(coords) else 6
    if len(zero) > n_ext and zero[:n_ext] == list(range(n_ext)):
        # ORCA always lists the external modes first, then any imaginary ones,
        # then the rest; a genuine vibration that prints as 0.00 can follow.
        # Drop by position and keep it -- unless its vector is empty too.
        zero = zero[:n_ext]
        if any(not np.any(m["vector"]) for k, m in enumerate(modes) if k not in zero):
            raise ParseError(
                "A mode after the translations and rotations has an all-zero "
                "displacement vector, so it cannot be scored.")
    if len(zero) != n_ext:
        raise ParseError(
            f"Found {len(zero)} zero-frequency modes; a {'linear' if n_ext == 5 else 'nonlinear'} "
            f"molecule has {n_ext} translations and rotations. Either the Hessian was not "
            "projected cleanly or the file is mis-parsed, so these are not dropped by guesswork.")
    keep = [m for k, m in enumerate(modes) if k not in zero]
    if len(keep) != 3 * natoms - n_ext:
        raise ParseError(
            f"After removing {n_ext} zero modes, {len(keep)} remain; expected "
            f"{3 * natoms - n_ext}.")
    return keep


def _cartesian_modes(modes, masses, what):
    """Make every vector a unit-norm Cartesian displacement, in place.

    s[V_S] is invariant to a mode's overall scale but not to the RELATIVE
    displacement of its atoms, so a mass-weighted vector (q = M^1/2 x) scores
    differently from the Cartesian one Gaussian prints. Which one a file holds
    is read off the vectors with two independent tests:

      * orthogonality -- normal modes are orthogonal as mass-weighted vectors
        (Q^T Q diagonal), and as Cartesian displacements only under the mass
        metric (X^T M X diagonal);
      * zero net momentum (Eckart) -- a vibration has sum_i m_i x_i = 0 as a
        Cartesian displacement, sum_i sqrt(m_i) q_i = 0 when mass-weighted.

    Either test alone can be blind: in a linear molecule or a set with one mode
    per irrep, symmetry makes the vectors orthogonal under both metrics, and
    with near-equal masses both tests barely differ. So vectors are converted
    only when the mass-weighted reading fits DECISIVELY better on the two tests
    combined; otherwise they are taken as Cartesian, the convention both ORCA
    and the Molden format document. Vectors orthogonal under neither metric are
    not normal modes of one Hessian and are refused.

    """
    X = np.array([np.asarray(m["vector"], float).ravel() for m in modes]).T
    if not np.isfinite(X).all():
        raise ParseError(f"The {what} holds non-numeric (NaN/inf) displacements.")
    convert = False
    if len(modes) >= 2 and np.ptp(masses) > 0:
        m3 = np.repeat(masses, 3)
        # numpy 2.x on macOS Accelerate raises spurious FP warnings from matmul
        # on finite input; the inputs were just checked finite, so silence them.
        with np.errstate(all="ignore"):
            off_cart = _offdiag(X.T @ (m3[:, None] * X))
            off_mw = _offdiag(X.T @ X)
        if min(off_cart, off_mw) > _ORTHO_TOL:
            raise ParseError(
                f"The {what}'s mode vectors are not mutually orthogonal under either the "
                f"Cartesian or the mass-weighted convention (worst overlaps {off_cart:.3f} / "
                f"{off_mw:.3f}). They are not normal modes of one Hessian, or the file was "
                "mis-parsed.")
        mw_fit = off_mw + _net_momentum(X, np.sqrt(masses))
        cart_fit = off_cart + _net_momentum(X, masses)
        convert = mw_fit * _DECISIVE < cart_fit
    if convert:
        for m in modes:
            m["vector"] = np.asarray(m["vector"], float) / np.sqrt(masses)[:, None]
    _normalise(modes)


def _net_momentum(X, weights):
    """Median over modes of |sum_i w_i x_i| / sum_i w_i |x_i| -- ~0 when the
    weighting matches the vectors' convention. The median, not the max, so the
    rigid translations a 3N set may carry do not swamp the vibrations."""
    n = len(weights)
    per_mode = []
    for k in range(X.shape[1]):
        x = X[:, k].reshape(n, 3) * weights[:, None]
        scale = np.abs(x).sum()
        per_mode.append(np.linalg.norm(x.sum(axis=0)) / scale if scale > 0 else 0.0)
    return float(np.median(per_mode))


def _offdiag(G):
    d = np.sqrt(np.clip(np.diag(G), 1e-300, None))
    C = G / np.outer(d, d)
    np.fill_diagonal(C, 0.0)
    return float(np.abs(C).max())


def _normalise(modes):
    for m in modes:
        v = np.asarray(m["vector"], float)
        n = np.linalg.norm(v)
        m["vector"] = v / n if n > 0 else v


def parse_frequency_file(text):
    """ORCA .hess or Molden, told apart by content rather than file extension."""
    if text.lstrip().upper().startswith("[MOLDEN FORMAT]"):
        return parse_molden(text)
    if "$normal_modes" in text or "$orca_hessian_file" in text:
        return parse_orca_hess(text)
    raise ParseError(
        "Not recognised as an ORCA .hess file (no $normal_modes block) or a Molden "
        "file (no '[Molden Format]' header). For a Gaussian .log use GaussianParser.")


_COMMENT_RE = re.compile(r"[#!].*$")


def read_bond_file(text, natoms):
    """Bond list for an ORCA/Molden molecule, from a plain ``.bonds`` file.

    One line per atom in the gjf ``geom=connectivity`` row format --
    ``center partner [order] partner [order] ...``, 1-based -- so a bare
    ``i j`` per line works too. A token with a decimal point is a bond order
    and is ignored (they play no part in the V-score); write orders as
    ``2.0``, since a bare ``2`` reads as an atom index. ``#`` or ``!`` starts a
    comment. Returns 0-based, de-duplicated, sorted pairs.
    """
    bonds = []
    for n, line in enumerate(text.splitlines(), 1):
        parts = _COMMENT_RE.sub("", line).split()
        if not parts or parts[0].startswith("["):
            continue
        for tok in parts:
            if "." not in tok and not _is_int(tok):
                raise ParseError(f"bond file line {n}: {tok!r} is neither an atom index nor a bond order.")
        idx = [int(t) - 1 for t in parts if "." not in t]
        center, partners = idx[0], idx[1:]
        for j in partners:
            if j == center:
                continue
            for a in (center, j):
                if not 0 <= a < natoms:
                    raise ParseError(
                        f"bond file line {n} references atom {a + 1}, but the geometry has "
                        f"{natoms} atoms.")
            pair = tuple(sorted((center, j)))
            if pair not in bonds:
                bonds.append(pair)
    if not bonds:
        raise ParseError("The bond file lists no bonds.")
    return sorted(bonds)
