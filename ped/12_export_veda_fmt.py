"""
Step 12: Assemble a VEDA4-compatible .fmt input file for benzene, so the
real VEDA4.exe program (veda4/veda4e1.exe) can be run on this job directly,
instead of only this repo's from-scratch Python PED reimplementation.

VEDA4 has no custom Hessian/coordinate input format of its own -- it only
ingests a Gaussian-log-shaped text excerpt (a ".fmt" file, per
veda4/Veda_use.doc) containing exactly three blocks (reference:
veda4/examp_1/examp_1/A1-91.FMT):
  1. one geometry "...orientation:" table
  2. the Gaussian "Standard" (3-modes-per-group, "Atom  AN" header)
     harmonic-frequency/normal-coordinate block -- NOT the "HP"
     (5-per-group, "Coord Atom Element:") block the rest of ped/ uses
  3. a "Force constants in Cartesian coordinates:" lower-triangle Hessian
     dump in Fortran D+00 notation

data/logs/C6H6.log already contains authentic, Gaussian-computed versions
of blocks 1 and 2 verbatim -- they're sliced out below by text anchor, not
recomputed. Only block 3 is missing from the log (this MP2/3-21G job never
printed it, and no .chk/formchk exists for this legacy run) -- it is
synthesized here from H_cart.npy, the Cartesian Hessian
02_reconstruct_hessian.py already reconstructed purely from Gaussian's own
printed frequencies/eigenvectors.

H_cart.npy is in amu*s^-2, not Gaussian's Hartree/Bohr^2 convention: in
H = M L Lambda L^T M, L satisfies L^T M L = I with plain amu masses, so L
carries units of mass^-1/2 only -- no length unit at all (rescaling l_raw
cancels out of that normalization). This script converts units and
independently re-validates the frequency round-trip *in atomic units*
before trusting the result -- the amu/(rad/s) round-trip already done in
step 02 never exercises this conversion.

Requires: H_cart.npy, opt_coords_ang.npy, opt_symbols.txt, vibfreq.npy
          (01, 02), and data/logs/C6H6.log itself (for verbatim block 1/2
          text)
Output: C6H6.fmt, C6H6_with_dummy_raman.fmt (ped/, tracked). C6H6.log's
        Standard block has no Raman Activ/Depolar rows (freq-only job) but
        A1-91.FMT's working reference does -- the _with_dummy_raman variant
        inserts synthetic 0.0000 rows in case VEDA4's fixed-offset reader
        expects them; try both in the GUI. Also best-effort copied to
        veda4/C6H6_reconstructed*.fmt (outside the git repo, convenience
        only).
"""
import os
import sys
import math
import numpy as np
import scipy.constants as sc

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from src.utils import atomicMass
from src.parser import _values_after_dash_marker

HERE = os.path.dirname(__file__)
LOG_PATH = os.path.join(HERE, '..', 'data', 'logs', 'C6H6.log')
VEDA_DIR = os.path.join(HERE, '..', '..', '..', 'veda4')
C_CM_PER_S = 2.99792458e10
LABEL_WIDTH = 15  # " Frequencies --" etc, measured from A1-91.FMT


# --- Block 1: geometry, verbatim from the log's LAST orientation block --

def extract_geometry_block(lines, expected_coords):
    starts = [i for i, l in enumerate(lines) if "orientation:" in l]
    if not starts:
        raise ValueError("No '...orientation:' block found in log.")
    start = starts[-1]

    end = start + 1
    dash_count = 0
    while dash_count < 3:
        if lines[end].strip().startswith('---'):
            dash_count += 1
        end += 1
    block = lines[start:end]

    # Cross-check against opt_coords_ang.npy: confirms this is the SAME
    # geometry block step 01's GaussianParser used, not a stray duplicate
    # (the log has 5 repeated orientation blocks from the opt+freq job).
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
    print(f"Geometry cross-check vs opt_coords_ang.npy: max|delta| = {max_diff:.2e} "
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
    after 'IR Inten'; C6H6.log (freq-only, no raman/polar) has neither.
    Synthesize zero-filled rows matching the real line's measured field
    width, in case VEDA4's fixed-offset reader expects them (untested --
    see plan's open-risk list)."""
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


# --- Block 3: Cartesian force constants, synthesized from H_cart.npy ----

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


# --- Load pipeline data ---------------------------------------------------

symbols = [l.strip() for l in open(os.path.join(HERE, 'opt_symbols.txt'))]
natm = len(symbols)
H_cart = np.load(os.path.join(HERE, 'H_cart.npy'))                   # (36,36) amu*s^-2
vibfreq = np.load(os.path.join(HERE, 'vibfreq.npy'))                  # (30,) cm^-1
opt_coords_ang = np.load(os.path.join(HERE, 'opt_coords_ang.npy'))    # (12,3)
masses = np.array([atomicMass[s.lower()] for s in symbols])           # (12,) amu
Nvib = len(vibfreq)

with open(LOG_PATH, 'r') as f:
    log_lines = f.readlines()

geom_block = extract_geometry_block(log_lines, opt_coords_ang)
freq_block = extract_standard_freq_block(log_lines, natm, Nvib)
print(f"Extracted geometry block ({len(geom_block)} lines) and Standard "
      f"frequency block ({len(freq_block)} lines) verbatim from {os.path.basename(LOG_PATH)}.")

# --- Unit conversion: amu*s^-2 -> Hartree/Bohr^2 -------------------------
# H_cart carries no length unit (see module docstring), so only a mass/time
# conversion is needed: 1 amu/s^2 = u*a0^2/Eh Hartree/Bohr^2.
u = sc.u
a0 = sc.value('Bohr radius')
Eh = sc.value('Hartree energy')
FACTOR_AMU_S2_TO_HARTREE_BOHR2 = u * a0 ** 2 / Eh

H_AU = 0.5 * (H_cart + H_cart.T) * FACTOR_AMU_S2_TO_HARTREE_BOHR2
diag = np.diag(H_AU)
print(f"\nH_AU diagonal range: [{diag.min():.4f}, {diag.max():.4f}] Hartree/Bohr^2 "
      f"(A1-91.FMT reference range: ~0.03-0.7)")

# --- Independent atomic-unit round-trip validation, BEFORE writing anything -
me = sc.value('electron mass')
t_au = sc.value('atomic unit of time')

masses_au = masses * (u / me)
M_au = np.repeat(masses_au, 3)
Minv_sqrt_au = 1.0 / np.sqrt(M_au)
D_au = (Minv_sqrt_au[:, None] * H_AU) * Minv_sqrt_au[None, :]
evals_au = np.sort(np.linalg.eigvalsh(D_au))

near_zero_au = evals_au[:6]
vib_au = evals_au[6:]
omega_au = np.sign(vib_au) * np.sqrt(np.abs(vib_au))
freq_recovered = np.sort(omega_au / t_au / (2 * np.pi * C_CM_PER_S))
freq_original = np.sort(vibfreq)
max_diff = np.max(np.abs(freq_recovered - freq_original))

print(f"\nAtomic-unit null-space check (should be ~0, T/R modes): {np.round(near_zero_au, 8)}")
print(f"Atomic-unit frequency round-trip max discrepancy: {max_diff:.6f} cm^-1 "
      f"({'OK, proceeding' if max_diff < 2.0 else 'WARNING'})")
if max_diff >= 2.0:
    raise ValueError(
        f"Atomic-unit frequency round-trip failed (max diff {max_diff:.3f} "
        "cm^-1) -- unit conversion is likely wrong; refusing to write a "
        ".fmt file with an unverified Hessian.")

hessian_text = format_cartesian_hessian_block(H_AU)


# --- Assemble and write ---------------------------------------------------

def assemble_fmt(freq_lines):
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


outputs = {
    'C6H6.fmt': assemble_fmt(freq_block),
    'C6H6_with_dummy_raman.fmt': assemble_fmt(insert_dummy_raman(freq_block)),
}

print()
for name, text in outputs.items():
    path = os.path.join(HERE, name)
    with open(path, 'w') as f:
        f.write(text)
    print(f"Wrote {path} ({len(text)} bytes, {text.count(chr(10))} lines)")

    if os.path.isdir(VEDA_DIR):
        veda_name = name.replace('C6H6', 'C6H6_reconstructed')
        veda_path = os.path.join(VEDA_DIR, veda_name)
        try:
            with open(veda_path, 'w') as f:
                f.write(text)
            print(f"  also copied to {veda_path} (untracked convenience copy)")
        except OSError as e:
            print(f"  WARNING: could not copy to {veda_path}: {e}")
    else:
        print(f"  ({VEDA_DIR} not found -- skipped convenience copy)")

print("\nNext step (manual, cannot be automated): open C6H6.fmt (or the "
      "_with_dummy_raman variant) in veda4e1.exe and check whether VEDA4 "
      "accepts it and reproduces these 30 frequencies in its own "
      "recomputation. See the plan's open-risks list for what to try if it "
      "doesn't parse cleanly.")
