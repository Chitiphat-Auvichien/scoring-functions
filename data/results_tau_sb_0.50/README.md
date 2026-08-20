# tau_SB=0.50 output set (RETIRED as of 2026-08-20)

This folder is a frozen, fully self-consistent snapshot of canonical
`data/results/` exactly as it stood at commit `7c18326` -- the last commit
before `tau_SB`'s canonical default changed 0.50 -> 0.42 (`3427ff8`, "Change
canonical tau_SB default 0.50 -> 0.42 and regenerate CSV results"). Every one
of the 133 tracked files under `data/results/` at that commit was copied here
verbatim via `git show 7c18326:data/results/<file>` -- nothing recomputed,
nothing relabeled, byte-identical to that commit.

## Addendum: the switch turned out to be two bugs, not one

The same-day archival of this folder started because an initial audit found
the 0.50->0.42 switch (`3427ff8`) incomplete -- only 9 of the tau_SB-
sensitive CSVs had been relabeled, missing `C6H6_EMIT.csv` and every other
per-molecule EMIT/normal CSV. Digging further turned up a second, deeper bug:
`3427ff8` had only hand-edited `data/results/thresholds.json`'s `"tau_SB"`
key, never the actual source of truth (`Thresholds.tau_SB`'s dataclass
default in `src/classifier.py`, still 0.50 at that point) -- so that JSON
edit would have been silently wiped by the very next real
`--calibrate` run regardless of the file-coverage bug. Both are fixed as of
this same 2026-08-20 session (see `data/results_tau_sb_0.42/README.md` and
`IMPLEMENTATION_PLAN.md`'s Recent history for the full writeup) via a
genuine code-level default change plus a full `reproduce.py` rerun, not
another round of hand-editing.

## Why this exists

`data/results/`'s `tau_SB` default moved to 0.42 in two steps this session:
first an incomplete "cheap relabel" (`3427ff8`) that only re-derived 9 of the
tau_SB-sensitive CSVs from their already-computed `V_Stretch` column (missing,
among others, `C6H6_EMIT.csv` and every other per-molecule EMIT/normal CSV --
see the 2026-08-20 entry in `IMPLEMENTATION_PLAN.md`'s "Recent history" for
the full bug writeup), then a full `reproduce.py` rerun that made canonical
`data/results/` completely self-consistent at 0.42. Once that rerun
overwrites `data/results/` in place, the pre-change tau_SB=0.50 output is not
reconstructible from disk state alone (git history has it, but not without
digging) -- so it is archived here as a clean, dedicated snapshot, mirroring
the naming already used for `data/results_tau_sb_0.42/`.

**As of 2026-08-20, tau_SB=0.42 is canonical and tau_SB=0.50 is retired.**
See `data/results/thresholds.json`'s top-level `"tau_SB": 0.42` key and
`data/results_tau_sb_0.42/README.md`.

## What's in here

All 133 files that were tracked under `data/results/` at commit `7c18326`,
reproducing that commit's directory structure exactly (including the
`data/results/` prefix stripped, so e.g. `data/results/C6H6_EMIT.csv` at that
commit is here as `data/results_tau_sb_0.50/C6H6_EMIT.csv`). This includes
every per-molecule `*_normal.csv`, `*_EMIT.csv`/`*_EMIT_full.csv`/
`*_EMIT_full_cartesian.csv`, `*_full_ped_table.csv`, `library_scores.csv`,
the benzene diagnostic/reference CSVs, transferability confusion tables,
`thresholds.json`, CPU/sensitivity benchmarks -- the complete tau_SB=0.50
canonical results set, not just the tau_SB-sensitive subset.

## How this was built

```
git ls-tree -r 7c18326 --name-only -- data/results/   # 133 files
# for each path, `git show 7c18326:<path> > data/results_tau_sb_0.50/<path minus data/results/ prefix>`
```

No relabeling, no rescoring, no pipeline invocation -- a pure git-history
retrieval. Spot-checked `C6H6_EMIT.csv` byte-identical against
`git show 7c18326:data/results/C6H6_EMIT.csv`.

## Do not resurrect this as canonical

This folder is retired. If a future session needs to compare tau_SB=0.50 vs
0.42 behavior, read from here rather than reverting `data/results/` or
`thresholds.json` -- the 0.42 default is the locked decision (see
`IMPLEMENTATION_PLAN.md`).
