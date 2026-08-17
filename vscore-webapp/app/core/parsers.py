"""parsers.py -- readers for the three input paths, all producing one shape.

Every path converges on the same four values, which are exactly what
``ModeScorer`` needs and nothing more (verified against the reference
implementation in scoring-functions/):

    atoms   list[str]        element symbol per atom      -> atomic mass
    coords  (N,3) float      equilibrium geometry         -> COM, principal axes, bonds
    bonds   list[(i,j)]      0-based index pairs          -> s[V_S] ONLY
    modes   list[dict]       {'frequency','vector',...}   -> the mode itself

Frequency, reduced mass, force constant and irrep ride along as row metadata;
the scorer never reads them. The Cartesian Hessian -- the largest block in a
VEDA .fmt -- is not needed at all and is not parsed.

Three readers:
  * ``parse_gaussian_log``  geometry + modes from a .log (NO connectivity --
    Gaussian does not print it in the output; it lives in the input deck)
  * ``parse_connectivity``  bonds from a .com/.gjf deck
  * ``parse_vsc``           all four from one .vsc file

and one writer, ``write_vsc``, so a run can be replayed from a single small
file. The writer stores the ORIGINAL geometry, not the principal-axis-aligned
one, so a .vsc is a faithful record of what was uploaded.
"""

from __future__ import annotations

import re

import numpy as np

from .utils import get_symbol, get_atomic_number, atomicMass

VSC_VERSION = "1.0"

_DASH_RE = re.compile(r"^-+$")


class ParseError(ValueError):
    """Raised with a message intended to be shown directly to the user."""


# ======================================================================
# Shared Gaussian frequency-block helpers (ported from scoring-functions
# src/parser.py so the two agree line-for-line on block layout).
# ======================================================================
def _values_after_dash_marker(line):
    parts = line.split()
    for i, tok in enumerate(parts):
        if _DASH_RE.match(tok):
            return [float(x) for x in parts[i + 1:]]
    raise ValueError(f"no dash-marker token found in line: {line!r}")


def _parse_freq_block_header(lines, line_idx):
    freqs = _values_after_dash_marker(lines[line_idx])
    n = len(freqs)

    irreps = lines[line_idx - 1].split() if line_idx - 1 >= 0 else []
    if len(irreps) != n:
        irreps = [None] * n

    def _optional(offset):
        try:
            vals = _values_after_dash_marker(lines[line_idx + offset])
        except (ValueError, IndexError):
            return [None] * n
        return vals if len(vals) == n else [None] * n

    return freqs, irreps, _optional(1), _optional(2)


def _parse_hp_mode_vectors(lines, data_start, natoms, num_modes):
    temp = np.zeros((num_modes, natoms, 3))
    cur, count = data_start, 0
    while count < 3 * natoms and cur < len(lines):
        parts = lines[cur].split()
        if len(parts) < 3:
            break
        try:
            c_idx = int(parts[0]) - 1
            a_idx = int(parts[1]) - 1
            if a_idx >= natoms:
                cur += 1
                continue
            vals = [float(v) for v in parts[3:]]
            for m in range(min(len(vals), num_modes)):
                temp[m, a_idx, c_idx] = vals[m]
            count += 1
            cur += 1
        except ValueError:
            break
    return temp, cur


# ======================================================================
# Path 1a: Gaussian .log  ->  atoms, coords, modes
# ======================================================================
def parse_gaussian_log(text):
    """Geometry and normal modes from Gaussian output text.

    Returns ``{'atoms', 'coords', 'modes'}``. Deliberately no 'bonds' key --
    a Gaussian .log carries no connectivity, and inventing one here is exactly
    the silent-corruption failure the reference implementation refuses to make.
    """
    lines = text.splitlines()

    starts = [i for i, l in enumerate(lines) if "Standard orientation" in l]
    if not starts:
        starts = [i for i, l in enumerate(lines) if "Input orientation" in l]
    if not starts:
        raise ParseError(
            "No 'Standard orientation' or 'Input orientation' block found. "
            "This does not look like a Gaussian output file.")

    geom_line = starts[-1] + 5
    coords, symbols = [], []
    i = geom_line
    while i < len(lines) and "-----" not in lines[i]:
        parts = lines[i].split()
        if len(parts) < 6:
            break
        try:
            symbols.append(get_symbol(int(parts[1])))
            coords.append([float(parts[3]), float(parts[4]), float(parts[5])])
        except ValueError:
            pass
        i += 1

    natoms = len(coords)
    if natoms == 0:
        raise ParseError("Found an orientation block but could not read any atoms from it.")

    modes = _parse_log_modes(lines, geom_line, natoms)
    if not modes:
        raise ParseError(
            "No vibrational modes found. This looks like a Gaussian job without a "
            "frequency calculation -- the scorer needs displacement vectors, so run "
            "the job with 'freq' (ideally 'freq=hpmodes') and upload that output.")

    expected = {3 * natoms - 6, 3 * natoms - 5}
    if len(modes) not in expected:
        raise ParseError(
            f"Parsed {len(modes)} vibrational modes for {natoms} atoms, but expected "
            f"3N-6={3*natoms-6} (nonlinear) or 3N-5={3*natoms-5} (linear). "
            "The frequency block was mis-parsed or the file is truncated.")

    return {"atoms": symbols, "coords": np.array(coords), "modes": modes}


def _parse_log_modes(lines, geom_line, natoms):
    has_hp = any("Coord Atom Element:" in l for l in lines)
    freq_lines = [i for i, l in enumerate(lines)
                  if "Frequencies --" in l and i > geom_line]
    if not freq_lines:
        freq_lines = [i for i, l in enumerate(lines) if "Frequencies --" in l]

    modes = []
    for start in freq_lines:
        _parse_log_block(lines, start, natoms, modes, force_hp_only=has_hp)
    return modes


def _parse_log_block(lines, line_idx, natoms, modes, force_hp_only):
    try:
        freqs, irreps, mus, ks = _parse_freq_block_header(lines, line_idx)
    except (ValueError, IndexError):
        return

    n_block = len(freqs)
    data_start, is_hp = -1, False
    for i in range(line_idx + 1, min(len(lines), line_idx + 30)):
        line = lines[i].strip()
        if "Coord Atom Element:" in line:
            data_start, is_hp = i + 1, True
            break
        if "Atom" in line and "AN" in line:
            data_start, is_hp = i + 1, False
            break
        parts = line.split()
        if len(parts) > 3 and parts[0].isdigit() and parts[1].isdigit() and parts[2].isdigit():
            try:
                float(parts[3])
                data_start, is_hp = i, True
                break
            except ValueError:
                pass

    if data_start == -1:
        return
    # When HP blocks exist, the Standard block repeats the same frequencies at
    # lower precision -- taking both would double-count every mode.
    if force_hp_only and not is_hp:
        return

    if is_hp:
        vecs, _ = _parse_hp_mode_vectors(lines, data_start, natoms, n_block)
        for m in range(n_block):
            modes.append({"frequency": freqs[m], "vector": vecs[m], "is_emit": False,
                          "reduced_mass": mus[m], "force_constant": ks[m],
                          "irrep": irreps[m]})
        return

    block = [[] for _ in range(n_block)]
    cur, count = data_start, 0
    while count < natoms and cur < len(lines):
        parts = lines[cur].split()
        try:
            raw = [float(v) for v in parts[2:]]
            for m in range(n_block):
                c = m * 3
                if c + 3 <= len(raw):
                    block[m].append(raw[c:c + 3])
            count += 1
            cur += 1
        except (ValueError, IndexError):
            break

    for m in range(n_block):
        if len(block[m]) == natoms:
            modes.append({"frequency": freqs[m], "vector": np.array(block[m]),
                          "is_emit": False, "reduced_mass": mus[m],
                          "force_constant": ks[m], "irrep": irreps[m]})


# ======================================================================
# Path 1b: Gaussian .com / .gjf  ->  bonds
# ======================================================================
def parse_connectivity(text, natoms=None):
    """Bond list from a Gaussian input deck's geom=connectivity section.

    Format is ``center  partner [order]  partner [order] ...``; bond orders are
    tolerated and ignored (they play no part in the V-score). Returns 0-based,
    de-duplicated, sorted index pairs.
    """
    lines = [l.rstrip() for l in text.splitlines()]

    # The charge/multiplicity line ("0 1") marks the start of the molecule
    # specification; connectivity follows the blank line after the geometry.
    idx, found = 0, False
    for i, line in enumerate(lines):
        parts = line.split()
        if len(parts) == 2 and _is_int(parts[0]) and _is_int(parts[1]):
            idx, found = i + 1, True
            break
    if not found:
        raise ParseError(
            "No charge/multiplicity line (e.g. '0 1') found in the input deck, so the "
            "connectivity section could not be located. Is this really a .com/.gjf file?")

    while idx < len(lines) and lines[idx].strip():
        idx += 1
    idx += 1

    bonds = []
    while idx < len(lines):
        line = lines[idx].strip()
        if not line:
            break
        parts = line.split()
        if parts and parts[0].isdigit():
            center = int(parts[0]) - 1
            for k in range(1, len(parts), 2):
                try:
                    neighbour = int(parts[k]) - 1
                except ValueError:
                    continue
                if neighbour == center:
                    continue
                pair = tuple(sorted((center, neighbour)))
                if pair not in bonds:
                    bonds.append(pair)
        idx += 1

    if not bonds:
        raise ParseError(
            "The input deck has no connectivity section. Re-run Gaussian with "
            "'geom=connectivity' in the route line, or upload a .vsc file with an "
            "explicit [CONNECTIVITY] block. Without bonds the V-score is undefined -- "
            "it would silently evaluate to 0.000 and label every mode as bending.")

    if natoms is not None:
        bad = [p for p in bonds if p[0] >= natoms or p[1] >= natoms or p[0] < 0]
        if bad:
            raise ParseError(
                f"Connectivity references atom index {max(max(p) for p in bad) + 1}, but "
                f"the geometry has only {natoms} atoms. The .com and .log do not describe "
                "the same molecule.")

    return sorted(bonds)


def _is_int(s):
    return s.lstrip("+-").isdigit()


# ======================================================================
# Path 2: the all-in-one .vsc file
# ======================================================================
_SECTION_RE = re.compile(r"^\s*\[([A-Za-z_]+)\]\s*(.*)$")


def _strip(line):
    """Drop a trailing '#' comment and surrounding whitespace."""
    return line.split("#", 1)[0].strip()


def parse_vsc(text):
    """Read a .vsc file. Returns ``{'atoms','coords','bonds','modes','title'}``."""
    raw = text.splitlines()
    title = ""
    sections, current = {}, None

    for line in raw:
        if line.lstrip().startswith("#"):
            head = line.lstrip("# \t")
            if head.upper().startswith("TITLE"):
                title = head[5:].strip()
            continue
        m = _SECTION_RE.match(line)
        if m:
            current = m.group(1).upper()
            sections[current] = {"arg": m.group(2).strip(), "lines": []}
            continue
        if current:
            body = _strip(line)
            if body:
                sections[current]["lines"].append(body)

    for req in ("GEOMETRY", "MODES"):
        if req not in sections:
            raise ParseError(
                f"Missing required [{req}] section. A .vsc file needs [GEOMETRY], "
                "[CONNECTIVITY] and [MODES].")

    atoms, coords = _parse_vsc_geometry(sections["GEOMETRY"]["lines"])
    natoms = len(atoms)

    if "CONNECTIVITY" not in sections:
        raise ParseError(
            "Missing [CONNECTIVITY] section. Without bonds the V-score is undefined "
            "(it would evaluate to 0.000 and label every mode as bending), so this "
            "file is rejected rather than scored.")
    bonds = _parse_vsc_connectivity(sections["CONNECTIVITY"]["lines"], natoms)

    modes = _parse_vsc_modes(sections["MODES"]["lines"], natoms,
                             sections["MODES"]["arg"])

    return {"atoms": atoms, "coords": coords, "bonds": bonds,
            "modes": modes, "title": title}


def _parse_vsc_geometry(lines):
    atoms, coords = [], []
    for n, line in enumerate(lines, 1):
        parts = line.split()
        if len(parts) < 4:
            raise ParseError(
                f"[GEOMETRY] line {n} has {len(parts)} fields, expected at least 4 "
                f"('index element x y z' or 'element x y z'): {line!r}")
        # Tolerate a leading serial index.
        if len(parts) >= 5 and parts[0].isdigit():
            parts = parts[1:]
        tok = parts[0]
        if tok.isdigit():
            sym = get_symbol(int(tok))
            if sym == "X":
                raise ParseError(f"[GEOMETRY] line {n}: unknown atomic number {tok}.")
        else:
            sym = tok.capitalize()
            if sym.lower() not in atomicMass:
                raise ParseError(
                    f"[GEOMETRY] line {n}: unrecognised element {tok!r}; no atomic "
                    "mass on file. Element identity affects all seven scores, not "
                    "just the V-score, so this cannot be defaulted.")
        try:
            xyz = [float(v) for v in parts[1:4]]
        except ValueError:
            raise ParseError(f"[GEOMETRY] line {n}: coordinates are not numbers: {line!r}")
        atoms.append(sym)
        coords.append(xyz)
    if not atoms:
        raise ParseError("[GEOMETRY] section is empty.")
    return atoms, np.array(coords)


def _parse_vsc_connectivity(lines, natoms):
    bonds = []
    for n, line in enumerate(lines, 1):
        parts = line.split()
        if not parts:
            continue
        try:
            center = int(parts[0]) - 1
        except ValueError:
            raise ParseError(f"[CONNECTIVITY] line {n}: first field must be an atom index.")
        for k in range(1, len(parts)):
            tok = parts[k]
            if "." in tok:          # a bond order -- ignored
                continue
            try:
                neighbour = int(tok) - 1
            except ValueError:
                raise ParseError(
                    f"[CONNECTIVITY] line {n}: {tok!r} is neither an atom index nor a "
                    "bond order.")
            if neighbour == center:
                continue
            for a in (center, neighbour):
                if not 0 <= a < natoms:
                    raise ParseError(
                        f"[CONNECTIVITY] line {n} references atom {a + 1}, but the "
                        f"geometry has {natoms} atoms.")
            pair = tuple(sorted((center, neighbour)))
            if pair not in bonds:
                bonds.append(pair)
    if not bonds:
        raise ParseError("[CONNECTIVITY] section is empty -- at least one bond is required.")
    return sorted(bonds)


_MODE_HEAD_RE = re.compile(r"^mode\s+(\d+)\s*(.*)$", re.IGNORECASE)


def _parse_vsc_modes(lines, natoms, declared):
    modes, cur, rows = [], None, []

    def flush():
        if cur is None:
            return
        if len(rows) != natoms:
            raise ParseError(
                f"[MODES] mode {cur['n']} has {len(rows)} displacement rows but the "
                f"geometry has {natoms} atoms.")
        modes.append({"frequency": cur["freq"], "vector": np.array(rows),
                      "is_emit": False, "reduced_mass": cur["mu"],
                      "force_constant": cur["k"], "irrep": cur["irrep"]})

    for line in lines:
        head = _MODE_HEAD_RE.match(line)
        if head:
            flush()
            rows = []
            meta = _parse_kv(head.group(2))
            # Every one of freq/mu/k/irrep is optional -- the scorer reads none
            # of them, they are row metadata only. Absent freq stays None rather
            # than defaulting to 0.0, so the UI shows a dash instead of an
            # invented 0.00 cm-1.
            cur = {"n": int(head.group(1)),
                   "freq": meta.get("freq"), "mu": meta.get("mu"),
                   "k": meta.get("k"), "irrep": meta.get("irrep")}
            continue
        if cur is None:
            raise ParseError(
                f"[MODES] displacement row before any 'mode N' header: {line!r}")
        parts = line.split()
        if len(parts) >= 4 and parts[0].isdigit():
            parts = parts[1:]
        if len(parts) < 3:
            raise ParseError(f"[MODES] mode {cur['n']}: expected 'dx dy dz', got {line!r}")
        try:
            rows.append([float(v) for v in parts[:3]])
        except ValueError:
            raise ParseError(
                f"[MODES] mode {cur['n']}: displacements are not numbers: {line!r}")
    flush()

    if not modes:
        raise ParseError("[MODES] section contains no modes.")
    if declared:
        try:
            want = int(declared.split()[0])
        except (ValueError, IndexError):
            want = None
        if want is not None and want != len(modes):
            raise ParseError(
                f"[MODES] header declares {want} modes but {len(modes)} were found.")
    return modes


def _parse_kv(s):
    """'freq=667.64 irrep=A'' -> {'freq': 667.64, 'irrep': "A'"}"""
    out = {}
    for tok in s.split():
        if "=" not in tok:
            continue
        key, val = tok.split("=", 1)
        key = key.lower()
        if key in ("freq", "mu", "k"):
            try:
                out[key] = float(val)
            except ValueError:
                pass
        elif key == "irrep":
            out["irrep"] = val
    return out


# ======================================================================
# The .vsc writer
# ======================================================================
def write_vsc(atoms, coords, bonds, modes, title="", source=""):
    """Serialise a parsed molecule to .vsc text.

    Callers must pass the ORIGINAL frame, not a principal-axis-aligned one.
    MIT() is NOT idempotent: its sign-fix heuristic negates the rotation when
    the heaviest atom's projected coordinates sum negative, so re-aligning an
    already-aligned geometry can flip axes and invert the T/R scores. Passing
    the frame the data arrived in makes a .vsc round trip exact.
    """
    out = [f"#VSCORE {VSC_VERSION}"]
    if title:
        out.append(f"#TITLE   {title}")
    if source:
        out.append(f"#SOURCE  {source}")
    out.append("")

    out.append("[GEOMETRY] Angstrom")
    for i, (sym, xyz) in enumerate(zip(atoms, coords), 1):
        out.append(f"{i:5d}  {sym:<2s} {xyz[0]:14.6f} {xyz[1]:14.6f} {xyz[2]:14.6f}")
    out.append("")

    out.append("[CONNECTIVITY]")
    partners = {i: [] for i in range(len(atoms))}
    for a, b in bonds:
        partners[a].append(b)
    for i in range(len(atoms)):
        row = "  ".join(str(j + 1) for j in sorted(partners[i]))
        out.append(f"{i + 1:5d}  {row}".rstrip())
    out.append("")

    out.append(f"[MODES] {len(modes)}")
    for n, m in enumerate(modes, 1):
        # Emit only the metadata that exists. All of it is optional on the way
        # in, so writing "freq=0.0000" for a mode that never had a frequency
        # would fabricate one on the way out.
        meta = []
        if m.get("frequency") is not None:
            meta.append(f"freq={m['frequency']:.4f}")
        if m.get("reduced_mass") is not None:
            meta.append(f"mu={m['reduced_mass']:.4f}")
        if m.get("force_constant") is not None:
            meta.append(f"k={m['force_constant']:.4f}")
        if m.get("irrep"):
            meta.append(f"irrep={m['irrep']}")
        out.append(f"  mode {n}" + ("   " + "   ".join(meta) if meta else ""))
        for i, d in enumerate(np.asarray(m["vector"]), 1):
            out.append(f"{i:5d} {d[0]:12.5f} {d[1]:12.5f} {d[2]:12.5f}")
    out.append("")
    return "\n".join(out)


# ======================================================================
# Precision guard
# ======================================================================
def displacement_precision(modes, sample=4000):
    """Smallest number of decimal places that reproduces every displacement.

    Standard Gaussian output prints 2 dp; ``freq=hpmodes`` prints 5. This
    matters because ``Tscore()`` unit-normalises EACH atom's displacement
    before summing, so an atom printed as '0.00 0.00 0.00' is either masked to
    zero or promoted to a full unit vector depending on rounding noise.
    Measured on real molecules, 2 dp input moves individual T-scores by up to
    0.31 on a [-1,1] scale (V and the S/B/SB labels survive, since the
    stretch/bend gap is 0.73 wide).
    """
    vals = []
    for m in modes:
        vals.extend(np.asarray(m["vector"]).ravel().tolist())
        if len(vals) > sample:
            break
    vals = np.array(vals[:sample])
    for dp in range(1, 7):
        if np.allclose(vals, np.round(vals, dp), atol=0, rtol=0):
            return dp
    return 6
