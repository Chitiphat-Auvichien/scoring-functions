#!/usr/bin/env python3
"""Merge real VEDA4 PED/TED output into this repo's mode-scoring/classification
output, producing a table of V-score (V_Stretch) vs. real %stretch/%bend
character per vibrational mode.

This is downstream of ped/build_veda_fmt.py: it does NOT run VEDA4 or build
.fmt files -- it consumes the .ved/.vdf pair a user has already produced by
manually opening a .fmt file in the real veda4e1.exe GUI and dropped into
data/ved/<mol>.ved + data/ved/<mol>.vdf (mirrors the existing per-molecule-
basename convention of data/logs/, data/gjf/, data/fchk/).

Sources combined:
  - data/results/<mol>_normal_classified.csv -- this repo's own T/R/V scores
    and S/B/M classification label (src.classifier.classify_all_modes,
    written by `python main.py -m <mol> --classify --mode normal`).
  - data/ved/<mol>.ved + .vdf -- real VEDA4 TED matrix + per-coordinate type
    definitions (STRE/BEND/TORS/OUT/...), produced by veda4e1.exe.

Category rule (explicit project decision): STRE -> stretch; every other
VEDA4 coordinate type (BEND, TORS, OUT, LIN, ...) -> bend.

Sign convention (revised 2026-08-07): the .ved file contains TWO separate
matrices, a "PED: sign = direction" table and a distinct "TED: sum = 100"
table (different values, not a transform of each other). This module uses
the PED table -- VEDA/VEDA4's own literature (Jamroz's papers, and every
published table built from VEDA output) is framed entirely around "PED
analysis"; TED is an internal VEDA4-only supplementary quantity, not what
the field reports as "%PED". Values are taken as abs() before summing per
category: the raw signed PED table does NOT reliably sum to 100 per row for
coupled/(quasi-)degenerate modes (e.g. CH4.ved's 1456.69 cm^-1 triple has a
row that sums to -98 signed), an artifact of the arbitrary rotation freedom
within a degenerate eigenspace -- taking abs() first (matching how %PED is
conventionally reported in the literature) recovers a ~100 sum instead
(that same row abs-sums to 100). This supersedes an earlier version of this
module that used the TED table -- see ped/archive_python_ped/13_compare_veda4.py
for the original PED-vs-TED table discovery (still accurate on the raw file
structure, just not on which table to report).

VEDA4's own mode-row order need not match data/results/<mol>_normal_classified.csv's
"Vib N" row order (e.g. VEDA4 often prints descending frequency, Gaussian's
parse order here is typically ascending) -- modes are matched by nearest
recomputed frequency, not by row position (see match_modes_by_frequency).

Usage:
    python ped/merge_ped_scores.py --molecule CH4
    python ped/merge_ped_scores.py --molecules CH4 H2O C6H6 --combined-output data/results/combined_ped_vs_scores.csv
    python ped/merge_ped_scores.py --all --combined-output data/results/combined_ped_vs_scores.csv
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _REPO_ROOT)

# Frequency-match warning tolerance (cm^-1). Reuses the same 2.0 cm^-1
# convention ped/build_veda_fmt.py / tests/test_veda_fmt_regression.py use
# for the Hessian round-trip discrepancy check -- VEDA4 recomputes its own
# frequencies from the imported Hessian/F-matrix, so a small (< 2 cm^-1,
# usually << 1 cm^-1) gap from Gaussian's original frequency is expected and
# is NOT itself a mismatch; this only warns (doesn't fail) when a match is
# looser than that.
DEFAULT_FREQ_TOL_CM1 = 2.0

_COMBINED_COLUMNS = [
    "Molecule", "Mode", "Freq", "V_Stretch", "label",
    "PED_Stretch_pct", "PED_Bend_pct", "VEDA_Freq", "Freq_Residual_cm-1",
]


# ---------------------------------------------------------------------------
# .ved parsing (TED matrix + frequency list)
# ---------------------------------------------------------------------------

def _read_matrix_block(lines, start_idx, n):
    """Read an n x n PED/TED data block starting at lines[start_idx] (the
    first data row, i.e. already past the '1 2 3 ... n' column-header line).

    Each logical row is printed as `<row_idx> <n coordinate values> <row_idx>`
    (VEDA4 repeats the row index as a trailing column) and, for molecules
    with many more coordinates than CH4/C6H6, MAY wrap across more than one
    physical text line -- tokens are accumulated across as many physical
    lines as needed until a full row (n+2 tokens) is collected, so this is
    robust to that wrapping even though it hasn't been observed for the two
    real molecules validated so far (CH4: 9x9, C6H6: 30x30, both single-
    line-per-row).

    Returns (matrix (n, n) ndarray, next_line_idx).
    """
    rows = []
    i = start_idx
    for expected_row in range(1, n + 1):
        tokens = []
        while len(tokens) < n + 2:
            if i >= len(lines):
                raise ValueError(
                    f"Ran out of lines while reading TED/PED row {expected_row} "
                    f"(collected {max(len(tokens) - 1, 0)} of {n} coordinate "
                    "values before hitting EOF).")
            tokens.extend(lines[i].split())
            i += 1
        if len(tokens) != n + 2:
            raise ValueError(
                f"Row {expected_row}: collected {len(tokens)} tokens, expected "
                f"exactly {n + 2} (row index + {n} coordinate values + repeated "
                "row index) -- matrix wrapping did not land on a row boundary.")
        row_idx = int(float(tokens[0]))
        trailing_idx = int(float(tokens[-1]))
        if row_idx != expected_row or trailing_idx != expected_row:
            raise ValueError(
                f"Row {expected_row}: expected leading/trailing row-index token "
                f"{expected_row}, got {row_idx}/{trailing_idx} -- matrix parse "
                "is out of sync.")
        rows.append([float(x) for x in tokens[1:-1]])
    return np.array(rows), i


def parse_ved(ved_path):
    """Parse a VEDA4 .ved file. Returns (freqs (n,), ped_matrix (n, n)).

    Only the PED block is parsed -- see module docstring for why TED is
    deliberately not used. `freqs` is in VEDA4's own row order (row i of
    ped_matrix corresponds to freqs[i]); it need not match the order of the
    log/classified.csv this molecule's frequencies came from originally.
    """
    with open(ved_path) as f:
        lines = f.readlines()

    def _find(pred, what):
        for idx, line in enumerate(lines):
            if pred(line):
                return idx
        raise ValueError(f"{ved_path}: could not find {what}.")

    dfac_idx = _find(lambda l: 'diagonality factor' in l, "the 'diagonality factor' line")
    ped_hdr_idx = _find(lambda l: l.strip().startswith('PED: sign'), "the 'PED: sign = direction' header")
    # Not parsed (see module docstring), but its presence is checked as a
    # sanity check that this is a genuine, complete VEDA4 .ved file.
    _find(lambda l: l.strip().startswith('TED:'), "the 'TED: sum = 100' header")

    freq_lines = [l for l in lines[dfac_idx + 1:ped_hdr_idx] if l.strip()]
    freqs = np.array([float(x) for l in freq_lines for x in l.split()])
    n = len(freqs)
    if n == 0:
        raise ValueError(f"{ved_path}: parsed zero frequencies between the "
                          "diagonality-factor line and the PED header.")

    # ped_hdr_idx+1 is the "1 2 3 ... n" column-header line; data starts after it.
    ped_matrix, _ = _read_matrix_block(lines, ped_hdr_idx + 2, n)
    return freqs, ped_matrix


# ---------------------------------------------------------------------------
# .vdf parsing (per-coordinate type: STRE vs everything else)
# ---------------------------------------------------------------------------

_COORD_DEF_RE = re.compile(r'^[a-zA-Z]\s*(\d+)\s+(\S+)')


def parse_coord_types(vdf_path):
    """Parse the 'definitions of modes' tail section of a .vdf file. Returns
    a list of VEDA4 coordinate-type strings (e.g. 'STRE', 'BEND', 'TORS',
    'OUT'), index 0 = coordinate s1, in order.

    Only the per-mode SUMMARY lines earlier in the .vdf show the single
    dominant coordinate per normal mode -- that's NOT what this reads (it
    isn't enough to compute a full %stretch/%bend sum). This reads the
    per-INTERNAL-COORDINATE type definitions instead, which are then applied
    to every column of the TED matrix parsed by parse_ved().

    Every definition line is `<letter> <idx> <TYPE> <atoms...> f<freq> <pct>
    [f<freq> <pct> ...]` (the trailing f/pct pairs list every mode this
    coordinate contributes to, and are not needed here -- only <idx>/<TYPE>
    are read). The leading letter is USUALLY 's', but VEDA4 also prints a
    'k' prefix for some coordinates on larger molecules (confirmed on
    C10H16's 72-coordinate set, e.g. 'k 32   STRE CC   ...') -- the letter
    itself carries no category meaning for our purposes, so any single
    letter is accepted, not just 's'.
    """
    with open(vdf_path) as f:
        lines = f.readlines()

    start = None
    for idx, line in enumerate(lines):
        if line.strip().startswith('definitions of modes'):
            start = idx
            break
    if start is None:
        raise ValueError(f"{vdf_path}: could not find the 'definitions of "
                          "modes' section.")

    types = {}
    for line in lines[start + 1:]:
        m = _COORD_DEF_RE.match(line.strip())
        if m:
            types[int(m.group(1))] = m.group(2)

    if not types:
        raise ValueError(f"{vdf_path}: found the 'definitions of modes' "
                          "header but no 's <idx> <TYPE> ...' lines after it.")
    n = max(types)
    missing = [i + 1 for i in range(n) if (i + 1) not in types]
    if missing:
        raise ValueError(f"{vdf_path}: missing coordinate-type definitions "
                          f"for s{missing} (found {len(types)} of {n}).")
    return [types[i + 1] for i in range(n)]


def compute_ped_percentages(ved_path, vdf_path):
    """Combine parse_ved + parse_coord_types into per-VEDA-mode-row
    PED_Stretch_pct/PED_Bend_pct. Returns a DataFrame (one row per VEDA mode,
    VEDA's own row order) with columns veda_freq, PED_Stretch_pct,
    PED_Bend_pct.

    Percentages are the ABSOLUTE VALUE of each PED entry, summed by
    category -- the raw signed PED table does not reliably sum to ~100 per
    row for coupled/(quasi-)degenerate modes (see module docstring), so
    abs() is applied first, matching how %PED is conventionally reported
    in the VEDA literature.
    """
    freqs, ped = parse_ved(ved_path)
    types = parse_coord_types(vdf_path)
    if len(types) != ped.shape[1]:
        raise ValueError(
            f"{vdf_path} defines {len(types)} coordinate types but "
            f"{ved_path}'s PED matrix has {ped.shape[1]} columns -- these "
            "must be the same VEDA4 run's paired output files.")

    abs_ped = np.abs(ped)
    stretch_mask = np.array([t == 'STRE' for t in types])
    stretch_pct = abs_ped[:, stretch_mask].sum(axis=1)
    bend_pct = abs_ped[:, ~stretch_mask].sum(axis=1)
    return pd.DataFrame({
        'veda_freq': freqs,
        'PED_Stretch_pct': stretch_pct,
        'PED_Bend_pct': bend_pct,
    })


# ---------------------------------------------------------------------------
# Frequency matching (VEDA row order != classified.csv "Vib N" row order)
# ---------------------------------------------------------------------------

def match_modes_by_frequency(classified_freqs, veda_freqs, tol=DEFAULT_FREQ_TOL_CM1,
                              context=""):
    """Greedy nearest-frequency one-to-one matching between
    data/results/<mol>_normal_classified.csv's "Vib N" row frequencies and
    VEDA4's own recomputed mode frequencies. Both lists must be the same
    length -- this is Vib-rows-vs-VEDA-mode-rows for ONE molecule, and a
    count mismatch means the two files don't describe the same normal-mode
    set, which is fatal (mirrors main.py's build_scorer_and_final ValueError
    style: fail loud rather than silently truncate/pad).

    Exact global-optimum assignment is not attempted -- for (near-)
    degenerate frequencies (e.g. CH4's triply-degenerate 1456.69 cm^-1
    trio) any valid pairing within the group gives essentially the same
    %stretch/%bend by symmetry, so a simple greedy nearest-neighbor match is
    sufficient (see module usage docs).

    Returns a dict {classified_row_idx: (veda_row_idx, residual_cm1)}.
    Prints (does not raise) a warning for any match whose residual exceeds
    `tol`.
    """
    n_c, n_v = len(classified_freqs), len(veda_freqs)
    if n_c != n_v:
        raise ValueError(
            f"{context}Mode count mismatch: {n_c} 'Vib N' rows in "
            f"classified.csv vs {n_v} mode rows in the .ved file -- cannot "
            "match them 1:1. Check both come from the same molecule/level "
            "of theory (and that classified.csv wasn't regenerated for a "
            "different geometry after the VEDA4 run).")

    candidates = sorted(
        ((abs(cf - vf), ci, vi)
         for ci, cf in enumerate(classified_freqs)
         for vi, vf in enumerate(veda_freqs)),
        key=lambda t: t[0],
    )
    used_c, used_v = set(), set()
    assignment = {}
    for resid, ci, vi in candidates:
        if ci in used_c or vi in used_v:
            continue
        used_c.add(ci)
        used_v.add(vi)
        assignment[ci] = (vi, resid)
        if len(assignment) == n_c:
            break

    for ci, (vi, resid) in assignment.items():
        if resid > tol:
            print(f"WARNING: {context}classified freq {classified_freqs[ci]:.2f} "
                  f"cm^-1 matched to VEDA freq {veda_freqs[vi]:.2f} cm^-1, "
                  f"residual {resid:.2f} cm^-1 exceeds tolerance {tol} cm^-1.")
    return assignment


# ---------------------------------------------------------------------------
# Per-molecule merge
# ---------------------------------------------------------------------------

def _resolve_ved_paths(molecule, data_dir):
    ved_path = os.path.join(data_dir, 'ved', f'{molecule}.ved')
    vdf_path = os.path.join(data_dir, 'ved', f'{molecule}.vdf')
    missing = [p for p in (ved_path, vdf_path) if not os.path.isfile(p)]
    return ved_path, vdf_path, missing


def merge_molecule_ped(molecule, repo_root=_REPO_ROOT, freq_tol=DEFAULT_FREQ_TOL_CM1,
                        write=True):
    """Merge one molecule's data/ved/<mol>.ved+.vdf into its
    data/results/<mol>_normal_classified.csv.

    Returns (full_df, vib_rows, out_path):
      - full_df: the classified.csv DataFrame with exactly two new columns
        appended, PED_Stretch_pct and PED_Bend_pct (NaN for the Tx/Ty/Tz/
        Rx/Ry/Rz ideal-reference rows, which have no PED). This is written
        to data/results/<mol>_normal_classified_ped.csv when write=True.
        The original _classified.csv is never modified.
      - vib_rows: list of dicts (one per "Vib N" row only) with the
        _COMBINED_COLUMNS fields, for combined multi-molecule tables.
      - out_path: path the per-molecule CSV was (or would be) written to.

    Raises FileNotFoundError if the classified.csv or the .ved/.vdf pair is
    missing, ValueError if the Vib-row count doesn't match the VEDA mode
    count (see match_modes_by_frequency).
    """
    data_dir = os.path.join(repo_root, 'data')
    classified_path = os.path.join(data_dir, 'results', f'{molecule}_normal_classified.csv')
    if not os.path.isfile(classified_path):
        raise FileNotFoundError(
            f"{classified_path} not found -- run "
            f"`python main.py -m {molecule} --classify --mode normal` first.")

    ved_path, vdf_path, missing = _resolve_ved_paths(molecule, data_dir)
    if missing:
        raise FileNotFoundError(
            f"Missing VEDA4 output for {molecule!r}: {missing}. Drop the "
            f"<{molecule}>.ved + <{molecule}>.vdf pair from a VEDA4 GUI run "
            "into data/ved/ (see ped/README.md).")

    df = pd.read_csv(classified_path)
    if 'Freq' not in df.columns:
        raise ValueError(
            f"{classified_path} has no 'Freq' column -- is this really a "
            "'normal' (not 'EMIT') classified.csv? PED merge only supports "
            "normal-mode classifications.")

    is_vib = df['Mode'].astype(str).str.startswith('Vib ')
    vib_idx = df.index[is_vib]
    if len(vib_idx) == 0:
        raise ValueError(f"{classified_path}: no 'Vib N' rows found.")

    ped_df = compute_ped_percentages(ved_path, vdf_path)
    assignment = match_modes_by_frequency(
        df.loc[vib_idx, 'Freq'].to_numpy(), ped_df['veda_freq'].to_numpy(),
        tol=freq_tol, context=f"{molecule}: ")

    df['PED_Stretch_pct'] = np.nan
    df['PED_Bend_pct'] = np.nan

    vib_rows = []
    # Iterate local_ci in increasing order (i.e. in vib_idx/classified.csv row
    # order), not assignment's insertion order (which is by increasing
    # residual from the greedy match) -- keeps vib_rows/output in the same
    # row order as classified.csv.
    for local_ci in sorted(assignment):
        vi, resid = assignment[local_ci]
        orig_idx = vib_idx[local_ci]
        stretch = ped_df.loc[vi, 'PED_Stretch_pct']
        bend = ped_df.loc[vi, 'PED_Bend_pct']
        veda_freq = ped_df.loc[vi, 'veda_freq']
        df.loc[orig_idx, 'PED_Stretch_pct'] = stretch
        df.loc[orig_idx, 'PED_Bend_pct'] = bend
        vib_rows.append({
            'Molecule': molecule,
            'Mode': df.loc[orig_idx, 'Mode'],
            'Freq': df.loc[orig_idx, 'Freq'],
            'V_Stretch': df.loc[orig_idx, 'V_Stretch'],
            'label': df.loc[orig_idx, 'label'],
            'PED_Stretch_pct': stretch,
            'PED_Bend_pct': bend,
            'VEDA_Freq': veda_freq,
            'Freq_Residual_cm-1': resid,
        })

    out_path = os.path.join(data_dir, 'results', f'{molecule}_normal_classified_ped.csv')
    if write:
        df.to_csv(out_path, index=False, float_format='%.4f')

    return df, vib_rows, out_path


# ---------------------------------------------------------------------------
# Batch / combined table
# ---------------------------------------------------------------------------

def _read_roster(repo_root):
    roster_path = os.path.join(repo_root, 'data', 'mol_list_method.csv')
    roster = pd.read_csv(roster_path)
    return roster['molecule'].tolist()


def build_combined_table(molecules, repo_root=_REPO_ROOT, freq_tol=DEFAULT_FREQ_TOL_CM1,
                          skip_missing=True, write=True):
    """Run merge_molecule_ped for each molecule in `molecules`. Molecules
    missing data/ved/<mol>.ved or .vdf (or their classified.csv) are skipped
    with a printed message if skip_missing=True (used by --all/--molecules),
    or raise if skip_missing=False (used by a single explicit --molecule so
    a typo/missing-input request fails loud rather than silently vanishing).

    Returns (combined_df, per_molecule_paths).
    """
    all_rows = []
    per_molecule_paths = []
    for mol in molecules:
        try:
            _df, vib_rows, out_path = merge_molecule_ped(mol, repo_root, freq_tol, write=write)
        except FileNotFoundError as e:
            if skip_missing:
                print(f"Skipping {mol!r}: {e}")
                continue
            raise
        per_molecule_paths.append(out_path)
        all_rows.extend(vib_rows)

    combined_df = pd.DataFrame(all_rows, columns=_COMBINED_COLUMNS)
    return combined_df, per_molecule_paths


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument('--molecule', help="Single molecule (e.g. CH4). Fails loud "
                        "if its classified.csv or .ved/.vdf pair is missing.")
    group.add_argument('--molecules', nargs='+', help="Multiple molecules "
                        "(space-separated). Skips (with a message) any that "
                        "are missing inputs, rather than aborting the batch.")
    group.add_argument('--all', action='store_true', help="Read the roster "
                        "from data/mol_list_method.csv's 'molecule' column; "
                        "skips (with a message) any missing data/ved/<mol>.ved "
                        "or .vdf.")
    p.add_argument('--combined-output', default=None,
                    help="Also write one combined table across all requested "
                         "molecules to this path (e.g. "
                         "data/results/combined_ped_vs_scores.csv).")
    p.add_argument('--freq-tol', type=float, default=DEFAULT_FREQ_TOL_CM1,
                    help=f"Warn (not fail) if a matched Vib/VEDA frequency "
                         f"pair's residual exceeds this many cm^-1 (default "
                         f"{DEFAULT_FREQ_TOL_CM1}).")
    return p.parse_args(argv)


def main(argv=None):
    args = _parse_args(argv)

    if args.all:
        molecules = _read_roster(_REPO_ROOT)
        skip_missing = True
    elif args.molecules:
        molecules = args.molecules
        skip_missing = True
    else:
        molecules = [args.molecule]
        skip_missing = False

    try:
        combined_df, per_molecule_paths = build_combined_table(
            molecules, _REPO_ROOT, freq_tol=args.freq_tol, skip_missing=skip_missing)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")
        return 1

    for p in per_molecule_paths:
        print(f"Wrote {p}")

    if args.combined_output:
        combined_df.to_csv(args.combined_output, index=False, float_format='%.4f')
        print(f"Wrote {len(combined_df)}-row combined table -> {args.combined_output}")
    elif len(molecules) > 1:
        print("(no --combined-output given; combined table not written)")

    if len(per_molecule_paths) == 0:
        print("No molecules were successfully merged.")
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
