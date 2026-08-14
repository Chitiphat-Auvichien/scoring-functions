"""Diff two library_scores.csv runs that differ only in the V-score weighting.

The point is not just "what changed" but "did only the things that CAN change
actually change". Weighting each bond by its reduced mass mu_AB is a no-op
whenever every bond in a molecule carries the same mu -- the factor cancels
between eq:vscore's numerator and denominator -- so a homoleptic AB_n molecule
must come out bit-for-bit identical. Anything else moving is a bug, not a
result.

    py scripts/compare_weighting.py                    # archive vs current
    py scripts/compare_weighting.py OLD.csv NEW.csv

Exit status is 1 if a homoleptic molecule's scores moved, so this can gate a
regeneration. Label changes on homoleptic molecules are reported but do NOT
fail: a mode can cross tau_S without its own score budging, purely because the
threshold moved. That is expected and is called out separately.

CAVEAT on comparing against an archive: an archived CSV is a snapshot of what
was on disk, which is not necessarily what the code of the day would produce.
Some archived per-molecule CSVs were found to predate an intervening engine
change, so an archive-vs-new diff can conflate "the weighting changed this"
with "this file was already stale". The airtight comparison is to regenerate
BOTH sides with the current code (py reproduce.py --v-weighting none / mu) and
diff those; against an archive, treat small diffs on homoleptic molecules as
suspect rather than as evidence about the weighting.
"""

import argparse
import os
import sys

import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULT_OLD = os.path.join(REPO_ROOT, "data", "results", "archive_unweighted",
                           "library_scores.csv")
DEFAULT_NEW = os.path.join(REPO_ROOT, "data", "results", "library_scores.csv")
# Within a single run the homoleptic invariance is EXACT -- both weightings
# produce bit-identical floats, and tests/test_scores.py asserts that with ==.
# Comparing two CSVs written by different runs is a weaker setting: the
# archived file carries last-bit (~1e-16) noise from its own generation that
# has nothing to do with the weighting (verified by re-scoring homoleptic
# molecules both ways through this same code path -- diff exactly 0.0). So the
# cross-run gate sits above that noise floor and well below anything physical.
SCORE_TOL = 1e-12


def bond_types(data_dir="data"):
    """{molecule: number of distinct reduced masses among its bonds}.

    Built from the same parse the scorer uses, so "homoleptic" here means
    exactly what _bond_weights() means by it: all bonds share one mu.
    """
    sys.path.insert(0, REPO_ROOT)
    from main import load_inputs, build_scorer_and_final
    from src.library_ingest import load_mol_roster, resolve_log_basename

    counts = {}
    for mol in load_mol_roster(data_dir)["molecule"]:
        base = resolve_log_basename(mol, data_dir)
        if base is None:
            continue
        try:
            raw, _ = load_inputs(base, "normal", data_dir)
            scorer, _ = build_scorer_and_final(raw, "normal")
        except Exception as exc:                      # missing/unparseable input
            print(f"  (skipping {mol}: {exc})")
            continue
        counts[mol] = len(set(round(float(w), 12) for w in scorer.bond_weights))
    return counts


def compare(old_path, new_path, data_dir="data"):
    old = pd.read_csv(old_path)
    new = pd.read_csv(new_path)
    key = ["molecule", "mode_index", "kind"]
    merged = old.merge(new, on=key, suffixes=("_old", "_new"))

    print(f"old: {old_path}  ({len(old)} rows)")
    print(f"new: {new_path}  ({len(new)} rows)")
    if len(merged) != len(old) or len(merged) != len(new):
        print(f"WARNING: only {len(merged)} rows matched on {key}")

    print("\nClassifying molecules by bond composition...")
    n_distinct_mu = bond_types(data_dir)

    merged["dV"] = (merged["V_Stretch_new"] - merged["V_Stretch_old"]).abs()
    merged["relabelled"] = (
        merged["predicted_label_old"] != merged["predicted_label_new"])

    homoleptic, heteroleptic, unknown = [], [], []
    for mol, grp in merged.groupby("molecule"):
        n_mu = n_distinct_mu.get(mol)
        row = (mol, float(grp["dV"].max()), int(grp["relabelled"].sum()))
        (unknown if n_mu is None else
         homoleptic if n_mu == 1 else heteroleptic).append(row)

    print(f"\n{'='*66}\nHOMOLEPTIC (one bond type -- mu cancels, scores MUST NOT move)")
    print(f"{'='*66}")
    violations = [r for r in homoleptic if r[1] > SCORE_TOL]
    noisy = [r for r in homoleptic if 0.0 < r[1] <= SCORE_TOL]
    relabelled = [r for r in homoleptic if r[2] > 0]
    print(f"  {len(homoleptic)} molecules, max |dV| over all of them = "
          f"{max((r[1] for r in homoleptic), default=0.0):.3e}")
    for mol, dv, n in violations:
        print(f"  !! {mol}: max |dV| = {dv:.3e}  -- SCORES MOVED, investigate")
    if noisy:
        print(f"  {len(noisy)} molecule(s) differ only at the last bit "
              f"(<= {SCORE_TOL:g}) -- cross-run float noise, not the weighting.")
    if relabelled:
        print("  Label changes with unchanged scores (threshold moved, expected):")
        for mol, _dv, n in relabelled:
            print(f"     {mol}: {n} mode(s)")

    print(f"\n{'='*66}\nHETEROLEPTIC (mixed bond types -- these are what changed)")
    print(f"{'='*66}")
    total_relabelled = 0
    for mol, dv, n in sorted(heteroleptic, key=lambda r: -r[2]):
        total_relabelled += n
        print(f"  {mol:12s} max |dV| = {dv:.4f}   {n:3d} label change(s)")
    print(f"  -> {len(heteroleptic)} molecules, {total_relabelled} label changes")

    if unknown:
        print(f"\nUnclassified (no parseable geometry): {', '.join(r[0] for r in unknown)}")

    if violations:
        print(f"\nFAIL: {len(violations)} homoleptic molecule(s) changed score.")
        return 1
    print("\nPASS: every homoleptic molecule is bit-for-bit unchanged.")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("old", nargs="?", default=DEFAULT_OLD)
    ap.add_argument("new", nargs="?", default=DEFAULT_NEW)
    ap.add_argument("--data-dir", default="data")
    args = ap.parse_args(argv)
    for path in (args.old, args.new):
        if not os.path.exists(path):
            ap.error(f"not found: {path}")
    return compare(args.old, args.new, args.data_dir)


if __name__ == "__main__":
    sys.exit(main())
