# Building VEDA4 input files

Given a Gaussian-computed geometry/frequency log for any molecule in `data/logs/`,
this folder builds a VEDA4-readable `.fmt` file, so PED (Potential Energy
Distribution) analysis comes from the **real VEDA4 program**, not from this
repo's own from-scratch calculation. `ped/build_veda_fmt.py` is the one script
you run.

## What a `.fmt` file is, and why VEDA4 needs it

VEDA4 has no custom Hessian/coordinate input format of its own -- it only ingests
a Gaussian-log-shaped text excerpt (a `.fmt` file, per `veda4/Veda_use.doc`)
containing exactly three blocks: (1) one geometry `"...orientation:"` table,
(2) the Gaussian "Standard" (3-modes-per-group, `"Atom  AN"` header)
harmonic-frequency/normal-coordinate block, and (3) a `"Force constants in
Cartesian coordinates:"` lower-triangle Hessian dump in Fortran D-notation.
VEDA4 auto-perceives connectivity from the geometry itself (confirmed from a
real VEDA4 session log), so no bond list needs to be encoded anywhere in the
`.fmt` file.

## Two frameworks for supplying the Hessian (block 3)

Blocks 1 and 2 always come verbatim from `data/logs/<mol>.log`, regardless of
framework. Only block 3 (the Hessian) differs:

| Framework | Source | Requires |
|---|---|---|
| `reconstruct` (Framework 1) | Recovers the full Cartesian Hessian purely from the log's own printed frequencies + Cartesian displacement eigenvectors (mass-weighted eigendecomposition inversion) | `data/logs/<mol>.log` only |
| `fchk` (Framework 2) | Parses the Hessian directly out of a Gaussian formatted checkpoint | `data/logs/<mol>.log` **and** `data/fchk/<mol>.fchk` (or `data/logs/<mol>.fchk`) |

**Honest caveat:** Framework 2 has **zero real test data** in this repo right
now -- no `.fchk` file exists anywhere in it. Its parsing logic is only
exercised by a synthetic, hand-built fixture in
`tests/test_veda_fmt_regression.py`. Treat it with suspicion on first real use,
until it's been run against a genuine Gaussian `.fchk` and independently
cross-checked (e.g. against Framework 1's reconstruction for the same job, or
against VEDA4's own recomputed frequencies).

## How to produce a `.fchk`

Add `%chk=<name>.chk` to the Gaussian route of the `freq` job, then, after the
job finishes, run Gaussian's `formchk <name>.chk <name>.fchk` utility to
convert the binary checkpoint into the text `.fchk` format `ped/hessian_fchk.py`
parses. (A less portable alternative is requesting a Hessian-printing route
option directly in the `.log` itself, avoiding the `.fchk` detour entirely --
but that's not what this module parses.)

## Usage

```bash
# Auto-detect: uses fchk if a .fchk resolves for the molecule, else reconstruct.
python ped/build_veda_fmt.py --molecule C6H6

# Force a specific framework (errors loudly if its required file is missing).
python ped/build_veda_fmt.py --molecule C6H6 --framework reconstruct
python ped/build_veda_fmt.py --molecule C6H6 --framework fchk

# Point at a .fchk explicitly instead of the default data/fchk/ search.
python ped/build_veda_fmt.py --molecule C6H6 --framework fchk --fchk-path /path/to/C6H6.fchk
```

Output `.fmt` files are written to `ped/output/` (gitignored, regenerable) by
default; override with `--output-dir`. A best-effort convenience copy is also
made to `--veda-dir` (default `../../../veda4` relative to `ped/`) if that
directory exists on disk -- never fatal if it doesn't.

## Verification

`tests/test_veda_fmt_regression.py` automatically checks: the round-trip
frequency discrepancy stays under 2.0 cm⁻¹, the freshly built C6H6 `.fmt`'s
Hessian floats match `ped/reference/C6H6.fmt` numerically, `blocks.fortran_d`
reproduces byte-exact reference values, and `molecule.resolve_molecule`'s
dual-extension (`.com`/`.gjf`) lookup works.

What it does **not** check automatically -- and can't -- is whether VEDA4
itself accepts the file: that requires manually opening the `.fmt` (or its
`_with_dummy_raman` variant) in `veda4e1.exe` and confirming it parses cleanly
and reproduces the reference frequencies in VEDA4's own recomputation.

## Relationship to the manuscript

The JCC manuscript's Table 6 previously reported PED percentages from this
repo's **own** redundant-internal-coordinate calculation (now archived in
`ped/archive_python_ped/`, which explicitly claimed no VEDA pass was needed).
**Project policy changed today (2026-08-02):** real VEDA4 is now the
authoritative PED source going forward. Table 6 is being updated separately to
reflect that (a parallel task, not part of this reorganization).

The archived Python calculation remains a valuable, documented cross-check: 3
of the 4 manuscript-cited modes agreed closely with real VEDA4, with one
genuine disagreement at 1056.39 cm⁻¹ (CCC-bend/CCH-bend character swapped
between the two methods -- a real rotation-ambiguity effect for a degenerate
mode pair, not an error in either method). See
`ped/archive_python_ped/ARCHIVE_NOTE.md` and
`ped/archive_python_ped/VEDA_STYLE_METHODOLOGY.md`/`PED_METHODOLOGY.md` for the
full historical record.

## Superseded pipeline

The prior from-scratch Python PED pipeline (redundant-coordinate `01`-`05`, the
abandoned VEDA-style reimplementation `06`-`11`, and the original benzene-only
VEDA4 bridge `12`-`14`) lives in `ped/archive_python_ped/` -- see that
directory's own `ARCHIVE_NOTE.md`.
