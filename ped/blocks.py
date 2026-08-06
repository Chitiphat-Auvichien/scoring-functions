"""Build the three text blocks a VEDA4-readable .fmt file needs, and assemble/
write them. Generalized (molecule-agnostic, no module-level globals) from the
archived ped/archive_python_ped/12_export_veda_fmt.py -- the block-extraction,
Fortran D-notation formatting, and atomic-unit round-trip-validation logic
below is unchanged from that script; only the parameterization changed.

VEDA4 has no custom Hessian/coordinate input format of its own -- it only
ingests a Gaussian-log-shaped text excerpt (a ".fmt" file, per
veda4/Veda_use.doc) containing exactly three blocks (reference:
veda4/examp_1/examp_1/A1-91.FMT):
  1. one geometry "...orientation:" table
  2. the Gaussian "Standard" (3-modes-per-group, "Atom  AN" header)
     harmonic-frequency/normal-coordinate block -- NOT the "HP"
     (5-per-group, "Coord Atom Element:") block the rest of this repo uses
  3. a "Force constants in Cartesian coordinates:" lower-triangle Hessian
     dump in Fortran D+00 notation
"""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from src.parser import _values_after_dash_marker  # noqa: E402

C_CM_PER_S = 2.99792458e10
LABEL_WIDTH = 15  # " Frequencies --" etc, measured from A1-91.FMT


# --- Block 1: geometry, verbatim from the log's LAST orientation block --

def extract_geometry_block(lines, expected_coords):
    # Match src.parser's GaussianParser._parse_standard_orientation: prefer
    # "Standard orientation", fall back to "Input orientation" only if no
    # Standard block exists. A bare "orientation:" substring match is too
    # broad -- it also catches non-tabular sections like a freq job's
    # "Dipole orientation:" (5 atom-number/x/y/z rows, no '---' header/
    # footer), which the dash-counting loop below then runs off the end of
    # the file trying to close.
    starts = [i for i, l in enumerate(lines) if "Standard orientation" in l]
    if not starts:
        starts = [i for i, l in enumerate(lines) if "Input orientation" in l]
    if not starts:
        raise ValueError("No 'Standard orientation'/'Input orientation' "
                          "block found in log.")
    start = starts[-1]

    end = start + 1
    dash_count = 0
    while dash_count < 3:
        if lines[end].strip().startswith('---'):
            dash_count += 1
        end += 1
    block = lines[start:end]

    # Cross-check against the caller's own parsed coordinates: confirms this
    # is the SAME geometry block the caller used elsewhere, not a stray
    # duplicate (a log can have several repeated orientation blocks from an
    # opt+freq job).
    coords = []
    for l in block:
        parts = l.split()
        if len(parts) >= 6 and parts[0].isdigit():
            try:
                coords.append([float(parts[-3]), float(parts[-2]), float(parts[-1])])
            except ValueError:
                pass
    coords = np.array(coords)
    if coords.shape != expected_coords.shape:
        raise ValueError(f"Geometry block parsed {coords.shape[0]} atoms, "
                          f"expected {expected_coords.shape[0]}.")
    max_diff = np.max(np.abs(coords - expected_coords))
    print(f"Geometry cross-check vs parsed coordinates: max|delta| = {max_diff:.2e} "
          f"({'OK' if max_diff < 1e-4 else 'WARNING: not the same geometry block!'})")
    return block


# --- Block 2: Standard-format frequency block, verbatim -----------------

def extract_standard_freq_block(lines, natoms, nvib):
    preamble_idxs = [i for i, l in enumerate(lines) if "Harmonic frequencies (cm**-1)" in l]
    standard_start = None
    for i in preamble_idxs:
        window = "".join(lines[i:i + 20])
        if "Coord Atom Element:" not in window and "Atom" in window and "AN" in window:
            standard_start = i
            break
    if standard_start is None:
        raise ValueError("Could not locate the Standard-format ('Atom  AN' "
                          "header) frequency block in the log.")

    j = standard_start
    while not lines[j].split() or not all(tok.isdigit() for tok in lines[j].split()):
        j += 1
    block_end = j

    modes_seen = 0
    while modes_seen < nvib:
        block_end += 1  # mode-index row
        block_end += 1  # irrep row
        n_this = len(_values_after_dash_marker(lines[block_end]))  # Frequencies row
        block_end += 1
        block_end += 3  # Red. masses, Frc consts, IR Inten rows
        block_end += 1  # 'Atom  AN ...' header row
        block_end += natoms  # per-atom displacement rows
        modes_seen += n_this

    if modes_seen != nvib:
        raise ValueError(f"Standard frequency block covers {modes_seen} modes, "
                          f"expected exactly {nvib}.")
    return lines[standard_start:block_end]


def insert_dummy_raman(freq_block_lines):
    """A1-91.FMT's working reference has 'Raman Activ --'/'Depolar --' rows
    after 'IR Inten'; a freq-only (no raman/polar) log has neither.
    Synthesize zero-filled rows matching the real line's measured field
    width, in case VEDA4's fixed-offset reader expects them."""
    out = []
    for line in freq_block_lines:
        out.append(line)
        if line.strip().startswith("IR Inten"):
            n_vals = len(_values_after_dash_marker(line))
            data_width = len(line.rstrip("\n")) - LABEL_WIDTH
            w = data_width // n_vals
            zero_field = f"{0.0:>{w}.4f}"
            out.append(" Raman Activ --" + zero_field * n_vals + "\n")
            out.append(" Depolar     --" + zero_field * n_vals + "\n")
    return out


# --- Block 3: Cartesian force constants ----------------------------------

def fortran_d(x, width=14, decimals=6):
    """Format x as Gaussian's '0.dddddd D+-ee' Fortran D-notation,
    right-justified in `width` chars (mantissa always in [0.1, 1.0), unlike
    Python's default [1.0, 10.0) -- verified byte-exact against several
    A1-91.FMT reference values including a zero and a rounding-edge case)."""
    if x == 0.0:
        mantissa, exp = 0.0, 0
    else:
        exp = math.floor(math.log10(abs(x))) + 1
        mantissa = x / 10.0 ** exp
        r = round(mantissa, decimals)
        if abs(r) >= 1.0:
            exp += 1
            mantissa = x / 10.0 ** exp
        elif abs(r) < 0.1 and r != 0.0:
            exp -= 1
            mantissa = x / 10.0 ** exp
    sign = '+' if exp >= 0 else '-'
    return f"{mantissa:.{decimals}f}D{sign}{abs(exp):02d}".rjust(width)


def format_cartesian_hessian_block(H, block_size=5):
    n = H.shape[0]
    out = []
    for c_start in range(1, n + 1, block_size):
        c_end = min(c_start + block_size - 1, n)
        out.append("".join(f"{c:>14d}" for c in range(c_start, c_end + 1)))
        for r in range(c_start, n + 1):
            row_end = min(r, c_end)
            vals = "".join(fortran_d(H[r - 1, c - 1]) for c in range(c_start, row_end + 1))
            out.append(f"{r:>4d}{vals}")
    return "\n".join(out) + "\n"


# --- Atomic-unit round-trip validation -----------------------------------

def round_trip_validate_au(H_AU, masses_amu, vibfreq_cm1, tol_cm1=2.0):
    """Independently re-validate the frequency round-trip *in atomic units*
    before trusting a Hessian enough to write a .fmt file with it.

    H_AU: (3N,3N) Cartesian Hessian in Hartree/Bohr^2.
    masses_amu: (N,) atomic masses in amu.
    vibfreq_cm1: (Nvib,) the vibrational frequencies this Hessian is supposed
        to reproduce.

    Returns max_diff (cm^-1). Raises ValueError (refusing to let the caller
    proceed to writing an unverified Hessian) if max_diff >= tol_cm1.
    """
    import scipy.constants as sc

    vibfreq_cm1 = np.asarray(vibfreq_cm1)
    masses_amu = np.asarray(masses_amu)
    natoms = len(masses_amu)
    n_modes = len(vibfreq_cm1)
    n_tr = 3 * natoms - n_modes  # T/R null-space size: 5 (linear) or 6 (nonlinear)
    if n_tr not in (5, 6):
        raise ValueError(
            f"Unexpected translation+rotation null-space size {n_tr} "
            f"(natoms={natoms}, n_modes={n_modes}); expected 5 (linear) or "
            "6 (nonlinear) -- masses_amu/vibfreq_cm1 likely mismatched.")

    u = sc.u
    me = sc.value('electron mass')
    t_au = sc.value('atomic unit of time')

    masses_au = masses_amu * (u / me)
    M_au = np.repeat(masses_au, 3)
    Minv_sqrt_au = 1.0 / np.sqrt(M_au)
    H_AU_sym = 0.5 * (H_AU + H_AU.T)
    D_au = (Minv_sqrt_au[:, None] * H_AU_sym) * Minv_sqrt_au[None, :]
    evals_au = np.sort(np.linalg.eigvalsh(D_au))

    near_zero_au = evals_au[:n_tr]
    vib_au = evals_au[n_tr:]
    omega_au = np.sign(vib_au) * np.sqrt(np.abs(vib_au))
    freq_recovered = np.sort(omega_au / t_au / (2 * np.pi * C_CM_PER_S))
    freq_original = np.sort(vibfreq_cm1)
    max_diff = float(np.max(np.abs(freq_recovered - freq_original)))

    print(f"Atomic-unit null-space check (should be ~0, T/R modes): {np.round(near_zero_au, 8)}")
    print(f"Atomic-unit frequency round-trip max discrepancy: {max_diff:.6f} cm^-1 "
          f"({'OK, proceeding' if max_diff < tol_cm1 else 'WARNING'})")
    if max_diff >= tol_cm1:
        raise ValueError(
            f"Atomic-unit frequency round-trip failed (max diff {max_diff:.3f} "
            f"cm^-1 >= tolerance {tol_cm1}) -- unit conversion or the Hessian "
            "itself is likely wrong; refusing to write a .fmt file with an "
            "unverified Hessian.")
    return max_diff


# --- Assemble and write ---------------------------------------------------

def assemble_fmt(geom_block, freq_lines, hessian_text):
    text = ("".join(geom_block) + "\n\n" + "".join(freq_lines) +
            "\n Force constants in Cartesian coordinates:\n" + hessian_text)
    n_orient = text.count("orientation")
    if n_orient != 1:
        raise ValueError(
            f"Expected exactly 1 'orientation' occurrence in the assembled "
            f".fmt text, found {n_orient} -- VEDA4 reads the FIRST "
            "orientation-labeled block (per Veda_use.doc); a stray "
            "duplicate risks silently feeding it the wrong geometry.")
    return text


def write_fmt_outputs(molecule_name, framework_tag, geom_block, freq_block,
                       hessian_text, output_dir, dummy_raman="both",
                       veda_convenience_dir=None):
    """Write the .fmt variant(s) for `molecule_name`/`framework_tag`
    ("reconstruct" or "fchk") to `output_dir`, best-effort convenience-copying
    to `veda_convenience_dir` if it exists (never fatal if it doesn't).

    dummy_raman: "both" (default) writes both the plain and
        _with_dummy_raman variants; "yes"/"no" writes only one.

    Returns {"written": {filename: path}, "convenience": {filename: path},
             "convenience_skipped_reason": str or None}.
    """
    if dummy_raman not in ("both", "yes", "no"):
        raise ValueError(f"dummy_raman must be 'both'/'yes'/'no', got {dummy_raman!r}")

    os.makedirs(output_dir, exist_ok=True)

    variants = {}
    if dummy_raman in ("both", "no"):
        name = f"{molecule_name}_{framework_tag}.fmt"
        variants[name] = assemble_fmt(geom_block, freq_block, hessian_text)
    if dummy_raman in ("both", "yes"):
        name = f"{molecule_name}_{framework_tag}_with_dummy_raman.fmt"
        variants[name] = assemble_fmt(geom_block, insert_dummy_raman(freq_block), hessian_text)

    veda_dir_exists = veda_convenience_dir is not None and os.path.isdir(veda_convenience_dir)
    convenience_skipped_reason = None
    if veda_convenience_dir is not None and not veda_dir_exists:
        convenience_skipped_reason = f"{veda_convenience_dir} not found -- skipped convenience copy"

    written = {}
    convenience = {}
    for name, text in variants.items():
        path = os.path.join(output_dir, name)
        with open(path, 'w') as f:
            f.write(text)
        written[name] = path
        print(f"Wrote {path} ({len(text)} bytes, {text.count(chr(10))} lines)")

        if veda_dir_exists:
            veda_name = name.replace(molecule_name, f"{molecule_name}_reconstructed", 1)
            veda_path = os.path.join(veda_convenience_dir, veda_name)
            try:
                with open(veda_path, 'w') as f:
                    f.write(text)
                convenience[name] = veda_path
                print(f"  also copied to {veda_path} (untracked convenience copy)")
            except OSError as e:
                print(f"  WARNING: could not copy to {veda_path}: {e}")

    if convenience_skipped_reason:
        print(f"  ({convenience_skipped_reason})")

    return {
        "written": written,
        "convenience": convenience,
        "convenience_skipped_reason": convenience_skipped_reason,
    }
