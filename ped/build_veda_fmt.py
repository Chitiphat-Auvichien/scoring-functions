#!/usr/bin/env python3
"""Build a VEDA4-readable .fmt file for any molecule with a Gaussian
freq log in data/logs/. This is the one script a user runs -- it replaces
the archived, benzene-only ped/archive_python_ped/12_export_veda_fmt.py.

Two frameworks can supply the Cartesian Hessian (block 3 of the .fmt file);
either way, block 1 (geometry) and block 2 (frequency block) ALWAYS come
verbatim from the .log -- a .fchk has no such text blocks, it is an
*additional* input for Framework 2, not a replacement for the log:
  - reconstruct (Framework 1): recover the Hessian purely from the log's own
    printed frequencies/eigenvectors. Needs only data/logs/<mol>.log.
  - fchk (Framework 2): parse the Hessian directly out of a Gaussian
    formatted checkpoint. Needs data/logs/<mol>.log AND a .fchk (see
    ped/hessian_fchk.py's docstring for how to produce one).

Usage:
    python ped/build_veda_fmt.py --molecule C6H6
    python ped/build_veda_fmt.py --molecule C6H6 --framework reconstruct
    python ped/build_veda_fmt.py --molecule C6H6 --framework fchk --fchk-path path/to/C6H6.fchk
"""
import argparse
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_HERE)
sys.path.insert(0, _HERE)
sys.path.insert(0, _REPO_ROOT)

from src.utils import atomicMass  # noqa: E402

import blocks  # noqa: E402
import hessian_fchk  # noqa: E402
import hessian_reconstruct  # noqa: E402
from molecule import resolve_molecule  # noqa: E402

_DEFAULT_OUTPUT_DIR = os.path.join(_HERE, 'output')
# Same '../../../veda4' depth ped/archive_python_ped/12_export_veda_fmt.py
# uses, computed from THIS script's own location (ped/), which is one
# directory shallower than the archived copy.
_DEFAULT_VEDA_DIR = os.path.normpath(os.path.join(_HERE, '..', '..', '..', 'veda4'))


def _parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--molecule', required=True,
                    help="Molecule name matching data/logs/<name>.log (e.g. C6H6).")
    p.add_argument('--framework', choices=['auto', 'reconstruct', 'fchk'], default='auto',
                    help="'auto' (default) uses fchk if a .fchk resolves, else "
                         "reconstruct. Explicit 'reconstruct'/'fchk' error loudly "
                         "if their required file is missing.")
    p.add_argument('--fchk-path', default=None,
                    help="Override the .fchk search (data/fchk/<name>.fchk, "
                         "then data/logs/<name>.fchk).")
    p.add_argument('--dummy-raman', choices=['both', 'yes', 'no'], default='both',
                    help="Which .fmt variant(s) to write (default: both).")
    p.add_argument('--output-dir', default=_DEFAULT_OUTPUT_DIR,
                    help=f"Default: {_DEFAULT_OUTPUT_DIR}")
    p.add_argument('--veda-dir', default=_DEFAULT_VEDA_DIR,
                    help=f"Best-effort convenience-copy destination if it "
                         f"exists. Default: {_DEFAULT_VEDA_DIR}")
    return p.parse_args(argv)


def main(argv=None):
    args = _parse_args(argv)

    mol = resolve_molecule(args.molecule, _REPO_ROOT, fchk_path_override=args.fchk_path)

    if args.framework == 'auto':
        framework = 'fchk' if mol.fchk_path else 'reconstruct'
    elif args.framework == 'fchk':
        if not mol.fchk_path:
            raise FileNotFoundError(
                f"--framework fchk requested but no .fchk resolved for "
                f"{mol.name!r} (looked in data/fchk/{mol.name}.fchk, "
                f"data/logs/{mol.name}.fchk). See ped/hessian_fchk.py's "
                "docstring for how to produce one, or pass --fchk-path.")
        framework = 'fchk'
    else:
        framework = 'reconstruct'

    # Geometry/frequency block text ALWAYS comes from the .log, regardless
    # of which framework supplies the Hessian.
    with open(mol.log_path, 'r') as f:
        log_lines = f.readlines()

    data = hessian_reconstruct.load_geometry_and_modes(mol.log_path)
    symbols = data["atoms"]
    coords_ang = data["coords"]
    modes = data["modes"]
    natoms = len(symbols)

    vibfreq = np.array([m["frequency"] for m in modes])
    reduced_mass = np.array([m["reduced_mass"] for m in modes])
    l_raw = np.array([m["vector"] for m in modes])
    masses = np.array([atomicMass[s.lower()] for s in symbols])
    Nvib = len(vibfreq)

    print(f"=== Building .fmt for {mol.name} (framework: {framework}) ===")
    if framework == 'reconstruct':
        H_AU, max_diff = hessian_reconstruct.build_hessian(
            mol, symbols, coords_ang, vibfreq, masses, reduced_mass, l_raw)
    else:
        H_AU, max_diff = hessian_fchk.build_hessian(mol, symbols, vibfreq, masses)

    geom_block = blocks.extract_geometry_block(log_lines, coords_ang)
    freq_block = blocks.extract_standard_freq_block(log_lines, natoms, Nvib)
    print(f"Extracted geometry block ({len(geom_block)} lines) and Standard "
          f"frequency block ({len(freq_block)} lines) verbatim from "
          f"{os.path.basename(mol.log_path)}.")

    hessian_text = blocks.format_cartesian_hessian_block(H_AU)

    result = blocks.write_fmt_outputs(
        mol.name, framework, geom_block, freq_block, hessian_text,
        args.output_dir, dummy_raman=args.dummy_raman,
        veda_convenience_dir=args.veda_dir)

    print(f"\n=== Summary ===")
    print(f"Molecule: {mol.name}")
    print(f"Framework: {framework}")
    print(f"Round-trip discrepancy: {max_diff:.6f} cm^-1 (tolerance 2.0 cm^-1)")
    print(f"Files written: {list(result['written'].values())}")
    if result['convenience']:
        print(f"Convenience copies: {list(result['convenience'].values())}")
    elif result['convenience_skipped_reason']:
        print(f"Convenience copy skipped: {result['convenience_skipped_reason']}")

    print("\nNext step (manual, cannot be automated): open the .fmt file (or "
          "the _with_dummy_raman variant) in veda4e1.exe and check whether "
          "VEDA4 accepts it and reproduces these frequencies in its own "
          "recomputation.")
    return result


if __name__ == '__main__':
    main()
