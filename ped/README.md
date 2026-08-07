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

## Merging real VEDA4 PED into this repo's mode-scoring output

Once you've opened a `.fmt` file in `veda4e1.exe` (see above) and it has
converged, VEDA4 leaves a `<something>.ved` (full numeric PED/TED matrices)
and a paired `<something>.vdf` (abbreviated summary + per-internal-coordinate
type definitions) in its own working folder -- ignore the various
`skra.vdf`/`mskra.vdf`/`wskra.vdf`/etc. duplicates VEDA4 also drops there from
earlier internal optimization cycles; only the final, named pair matters.

Copy that `.ved` + `.vdf` pair into **`data/ved/<mol>.ved`** and
**`data/ved/<mol>.vdf`** (same basename convention as `data/logs/`,
`data/gjf/`, `data/fchk/` -- one pair per molecule). `ped/merge_ped_scores.py`
then combines them with this repo's own `data/results/<mol>_normal.csv` (T/R/V
scores + S/B/M classification label, from `python main.py -m <mol> --mode
normal`) into a per-mode table of V-score vs. real %stretch/%bend character:

```bash
# Single molecule -- appends PED_Stretch_pct/PED_Bend_pct (+ per-bond-type
# columns) into data/results/<mol>_normal.csv IN PLACE (an intentional
# enrichment of that one canonical per-molecule file, not a separate `_ped`
# file). Fails loud if inputs are missing.
python ped/merge_ped_scores.py --molecule CH4

# Several molecules, plus one combined table for correlating V-score against
# real PED %stretch across all of them. Molecules missing a .ved/.vdf pair
# are skipped with a message, not fatal to the batch.
python ped/merge_ped_scores.py --molecules CH4 H2O C6H6 \
    --combined-output data/results/combined_ped_vs_scores.csv

# Every molecule in data/mol_list_method.csv's roster that has a .ved/.vdf pair.
python ped/merge_ped_scores.py --all \
    --combined-output data/results/combined_ped_vs_scores.csv

# Full per-internal-coordinate |PED| table for one molecule (one row per VEDA
# mode, one column per internal coordinate, e.g. '1_STRE_CH') -- a developer
# diagnostic, not the aggregate Stretch/Bend/bond-type merge above. Reads
# only data/ved/<mol>.ved+.vdf, no <mol>_normal.csv needed. Writes
# data/results/<mol>_full_ped_table.csv (override with --full-table-output).
python ped/merge_ped_scores.py --molecule CH4 --full-table
```

The single-molecule and roster-batch merges are also wired into `main.py`'s
flag dispatch, alongside `--emit-projection`/`--library`/`--calibrate`/
`--figures`, so ordinary use doesn't require calling this script directly
(`--full-table` is developer-only and stays `ped/merge_ped_scores.py`-only):

```bash
# Per-molecule (equivalent to --molecule above); requires -m.
python main.py -m CH4 --ped-merge

# Roster batch + combined table (equivalent to --all above); global, ignores -m.
python main.py --ped-merge-all --combined-output data/results/combined_ped_vs_scores.csv
```

**Category rule:** VEDA4's `STRE` coordinate type maps to stretch; every
other type (`BEND`, `TORS`, `OUT`, `LIN`, ...) maps to bend. **The `PED` table
is used** (`"PED: sign = direction"` in the `.ved` file), not `TED`
(`"TED: sum = 100"`) -- VEDA4/VEDA's own literature (Jamroz's papers, and
every published table built from VEDA output) is framed entirely around "PED
analysis"; TED is an internal VEDA4-only supplementary quantity, not what the
field reports as "%PED" (see `ped/merge_ped_scores.py`'s module docstring for
the full rationale, including why an earlier version of this module used TED
instead). Values are taken as **absolute value** before summing per category
-- the raw signed PED table does not reliably sum to ~100 per row for
coupled/(quasi-)degenerate modes (an artifact of the arbitrary rotation
freedom within a degenerate eigenspace); abs()-summing first recovers a ~100
sum, matching how %PED is conventionally reported in the literature.

VEDA4's own mode-row order need not match `<mol>_normal.csv`'s `Vib N` row
order (e.g. VEDA4 often prints descending frequency; this repo's parse order
is typically ascending) -- `merge_ped_scores.py` matches modes by nearest
recomputed frequency (greedy, warns above 2.0 cm⁻¹ residual -- the same
round-trip tolerance convention `build_veda_fmt.py` uses), not by row
position, and fails loud if the Vib-row count and VEDA mode-row count differ.

The first real cross-check (CH4, `fchk` framework, 2026-08-05) validated the
whole story qualitatively: the two CH-stretch frequency groups (~3086,
~3195 cm⁻¹) score `V_Stretch` ≈ 0.995-1.0 and land at `PED_Stretch_pct` ≈
98-101%; the two degenerate bending groups (~1457, ~1668 cm⁻¹) score
`V_Stretch` ≈ 0.03 or lower and land at `PED_Bend_pct` ≈ 99-101%. The same
held for C6H6's mixed-character modes (e.g. the two `1056.39` cm⁻¹ modes,
labeled `SB` by this repo's classifier, show intermediate `PED_Stretch_pct`
of 1% and 87%). `data/ved/CH4.ved`/`.vdf` and `data/ved/C6H6.ved`/`.vdf` are
committed as genuine ground-truth regression fixtures (see
`tests/test_merge_ped_scores.py`), not throwaway output.

## Relationship to the manuscript

The JCC manuscript's Table 6 previously reported PED percentages from this
repo's **own** redundant-internal-coordinate calculation (now archived in
`ped/archive_python_ped/`, which explicitly claimed no VEDA pass was needed).
**Project policy changed today (2026-08-02):** real VEDA4 is now the
authoritative PED source going forward. Table 6 has been updated accordingly
to report VEDA4's TED values (see `JCC/JCC_man_scoring/JCC_SI_PED_benzene.tex`
for the full methodology and validation).

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
