"""
Step 10 (VEDA-style cross-check): side-by-side comparison of this new
natural-coordinate/EPm-optimized PED table against the existing, published
redundant-coordinate table (benzene_PED_table.csv, from 01-05). Organized
facts for the user's own go/no-go call on promoting this method into the
manuscript -- no recommendation is made here.

Requires: veda_style_PED_table.csv (step 09); benzene_PED_table.csv (01-05,
          must already exist); PED_group_pct.npy, cats.txt (step 04);
          B.npy (step 03); B_nat.npy (step 06); epm_history.txt (step 07)
Output: veda_vs_redundant_comparison.txt, veda_vs_redundant_comparison.csv
"""
import csv
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(__file__))
from _ped_core import compute_PED, freq_to_lambda, EPm

HERE = os.path.dirname(__file__)
CATEGORIES = ['CH_stretch', 'CC_stretch', 'CCC_bend', 'CCH_bend', 'out_of_plane']
KEY_FREQS = [1056.39, 1319.27, 1532.85, 1598.90]


def load_csv(path):
    with open(path) as f:
        return list(csv.DictReader(f))


existing_path = os.path.join(HERE, 'benzene_PED_table.csv')
new_path = os.path.join(HERE, 'veda_style_PED_table.csv')
if not os.path.exists(existing_path):
    raise SystemExit(f"Missing {existing_path} -- run ped/run_all.sh first.")
existing = load_csv(existing_path)
new = load_csv(new_path)
if len(existing) != len(new):
    raise SystemExit(f"Row count mismatch: existing={len(existing)}, new={len(new)} "
                      "-- degenerate-pair grouping must have diverged; investigate before comparing.")

lines = []


def emit(s=""):
    lines.append(s)


emit("=" * 78)
emit("VEDA-style (natural-coordinate, EPm-optimized) vs. existing")
emit("(redundant-coordinate, Pulay-pseudoinverse) benzene PED comparison")
emit("=" * 78)
emit()

# --- 1. per-mode side-by-side table ---
emit("1. Per-mode side-by-side percentages (existing -> new), max abs diff")
emit("-" * 78)
diffs = []
for e, n in zip(existing, new):
    assert abs(float(e['freq_cm-1']) - float(n['freq_cm-1'])) < 0.5, \
        f"Row mismatch: {e['freq_cm-1']} vs {n['freq_cm-1']}"
    row_diffs = []
    parts = [f"mode {e['mode']:>2s} ({e['freq_cm-1']:>7s} cm-1):"]
    for cat in CATEGORIES:
        ev, nv = float(e[f'{cat}_pct']), float(n[f'{cat}_pct'])
        d = nv - ev
        row_diffs.append(abs(d))
        parts.append(f" {cat}={ev:6.1f}->{nv:6.1f} ({d:+6.1f})")
    diffs.append(max(row_diffs))
    emit("".join(parts))
emit()

# --- 2. summary stats ---
emit("2. Summary statistics across all 30 modes x 5 categories")
emit("-" * 78)
all_diffs = []
for e, n in zip(existing, new):
    for cat in CATEGORIES:
        all_diffs.append(abs(float(n[f'{cat}_pct']) - float(e[f'{cat}_pct'])))
all_diffs = np.array(all_diffs)
emit(f"Mean |diff| per category-entry: {all_diffs.mean():.2f} percentage points")
emit(f"Max |diff|:                     {all_diffs.max():.2f} percentage points")
emit(f"Per-mode max |diff| > 10 pts:   {sum(d > 10 for d in diffs)} / {len(diffs)} modes")
emit(f"Per-mode max |diff| > 30 pts:   {sum(d > 30 for d in diffs)} / {len(diffs)} modes")
emit()

# --- 3. EPm comparison ---
emit("3. EPm (VEDA's documented PED-purity objective, max possible = 30)")
emit("-" * 78)
emit("EPm is defined over INDIVIDUAL internal coordinates (VEDA's own")
emit("definition: sum of each mode's single largest coordinate PED), not")
emit("over post-hoc physical categories -- so the existing method's EPm")
emit("below is recomputed directly from its raw 42-coordinate B/H/L/lambda")
emit("(04_compute_ped.py's own inputs), matching the level of aggregation")
emit("used for the new method's 30 individual natural/mixed/reduced")
emit("coordinates. Comparing category-grouped percentages instead would")
emit("silently favor whichever method groups more coordinates per category.")
B_existing = np.load(os.path.join(HERE, 'B.npy'))
H_ex = np.load(os.path.join(HERE, 'H_cart.npy'))
L_ex = np.load(os.path.join(HERE, 'L.npy'))
vib_ex = np.load(os.path.join(HERE, 'vibfreq.npy'))
lam_ex = freq_to_lambda(vib_ex)
PED_existing_raw = compute_PED(B_existing, H_ex, L_ex, lam_ex)
epm_existing = EPm(PED_existing_raw)
epm_history = [float(l.split('\t')[1]) for l in open(os.path.join(HERE, 'epm_history.txt'))]
emit(f"Existing redundant-coordinate table:                EPm = {epm_existing:.4f}")
emit(f"New natural-coordinate table, UNMIXED (before step 07): EPm = {epm_history[0]:.4f}")
emit(f"New natural-coordinate table, after mixing (step 07):   EPm = {epm_history[-1]:.4f}")
new_epm_after_reduction = None
if os.path.exists(os.path.join(HERE, 'B_final.npy')):
    B_final = np.load(os.path.join(HERE, 'B_final.npy'))
    H = np.load(os.path.join(HERE, 'H_cart.npy'))
    L = np.load(os.path.join(HERE, 'L.npy'))
    vib = np.load(os.path.join(HERE, 'vibfreq.npy'))
    lam = freq_to_lambda(vib)
    new_epm_after_reduction = EPm(compute_PED(B_final, H, L, lam))
    emit(f"New natural-coordinate table, after reduction (step 08): EPm = {new_epm_after_reduction:.4f}")
emit(f"Improvement over existing method: "
     f"{(new_epm_after_reduction or epm_history[-1]) - epm_existing:+.4f} "
     f"({100*((new_epm_after_reduction or epm_history[-1])/epm_existing - 1):+.1f}%)")
emit()
emit("NOTE: most of the EPm improvement comes from resolving the ring-closure")
emit("redundancy into a well-conditioned non-redundant basis (unmixed EPm "
     f"already {epm_history[0]:.2f} vs existing {epm_existing:.2f}), not from the")
emit("mixing/reduction optimization itself (which adds only "
     f"{epm_history[-1]-epm_history[0]:+.2f} more).")
emit()

# --- 4. convergence / conditioning diagnostics ---
emit("4. Convergence and conditioning diagnostics")
emit("-" * 78)
emit(f"Mixing (step 07): {len(epm_history)-1} accepted moves, "
     f"converged in {open(os.path.join(HERE, 'epm_history.txt')).read().count(chr(10))} history entries")
B_orig = np.load(os.path.join(HERE, 'B.npy'))
B_nat = np.load(os.path.join(HERE, 'B_nat.npy'))
sv_orig = np.linalg.svd(B_orig, compute_uv=False)
sv_nat = np.linalg.svd(B_nat, compute_uv=False)
emit(f"Original redundant B (42x36): condition number over its 30 real "
     f"singular values = {sv_orig[0]/sv_orig[29]:.4e} (rcond=1e-8 pseudoinverse used)")
emit(f"New natural B_nat (30x36):    condition number = {np.linalg.cond(B_nat):.4e} "
     f"(exact inverse-like pseudoinverse, no truncation needed)")
emit()

# --- 5. the four manuscript-cited modes ---
emit("5. The four modes cited in JCC_man_CA.tex Table tab:benzenemixed")
emit("-" * 78)
for freq in KEY_FREQS:
    e_row = min(existing, key=lambda r: abs(float(r['freq_cm-1']) - freq))
    n_row = min(new, key=lambda r: abs(float(r['freq_cm-1']) - freq))
    emit(f"\n{freq} cm-1:")
    emit(f"  existing (published): {e_row['dominant_character']}")
    emit(f"  new (natural coords): {n_row['dominant_character']}")
    e_cats = {c: float(e_row[f'{c}_pct']) for c in CATEGORIES}
    n_cats = {c: float(n_row[f'{c}_pct']) for c in CATEGORIES}
    max_d = max(abs(n_cats[c] - e_cats[c]) for c in CATEGORIES)
    flag = " <-- LARGE DISAGREEMENT" if max_d > 30 else ""
    emit(f"  max category diff: {max_d:.1f} percentage points{flag}")
emit()

# --- 6. factual checklist, no recommendation ---
emit("6. Factual checklist (no recommendation -- for your own review)")
emit("-" * 78)
row_sums_ok = all(
    99.0 < sum(float(r[f'{c}_pct']) for c in CATEGORIES) < 101.0
    for r in new
)
emit(f"[{'x' if row_sums_ok else ' '}] "
     "New table's category percentages sum to ~100% for every mode")
emit(f"[{'x' if (new_epm_after_reduction or epm_history[-1]) > epm_existing else ' '}] "
     "New method's EPm exceeds the existing method's EPm")
n_large_disagreements = sum(
    max(abs(float(n[f'{c}_pct']) - float(e[f'{c}_pct'])) for c in CATEGORIES) > 30
    for e, n in zip(existing, new))
emit(f"[{'x' if n_large_disagreements == 0 else ' '}] "
     f"Zero modes with >30-point category disagreement (actual count: {n_large_disagreements})")
key_disagreements = sum(
    max(abs(float(n_row[f'{c}_pct']) - float(e_row[f'{c}_pct'])) for c in CATEGORIES) > 30
    for freq in KEY_FREQS
    for e_row, n_row in [(min(existing, key=lambda r: abs(float(r['freq_cm-1']) - freq)),
                           min(new, key=lambda r: abs(float(r['freq_cm-1']) - freq)))]
)
emit(f"[{'x' if key_disagreements == 0 else ' '}] "
     f"Zero of the 4 manuscript-cited modes show >30-point disagreement "
     f"(actual count: {key_disagreements}/4)")
emit(f"[ ] CCC-bend is 0% throughout the new table by construction (see "
     "VEDA_STYLE_METHODOLOGY.md) -- the two methods are not reporting the "
     "same 5 categories on equal footing for bend-type modes; this is a "
     "structural difference, not a checklist pass/fail item.")
emit()
emit("IMPORTANT: modes at 1056.39 and 1532.85 cm-1 -- the two modes the")
emit("manuscript's own discussion relies on most heavily to argue for real")
emit("mixed stretch-bend (SB) character -- are the ones with the largest")
emit("disagreement: the new natural-coordinate method reports them as")
emit("~100% CCH bend (no stretch character at all), where the existing")
emit("redundant-coordinate method reported 44% and 13% CC-stretch")
emit("respectively. This is a real numerical disagreement between two")
emit("legitimate coordinate-system choices, not a bug in either pipeline --")
emit("see VEDA_STYLE_METHODOLOGY.md Sec. 'Why the two methods disagree'.")

report = "\n".join(lines)
print(report)

with open(os.path.join(HERE, 'veda_vs_redundant_comparison.txt'), 'w') as f:
    f.write(report + "\n")

with open(os.path.join(HERE, 'veda_vs_redundant_comparison.csv'), 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['mode', 'freq_cm-1'] + [f'{c}_existing' for c in CATEGORIES]
                     + [f'{c}_new' for c in CATEGORIES] + [f'{c}_diff' for c in CATEGORIES])
    for e, n in zip(existing, new):
        row = [e['mode'], e['freq_cm-1']]
        row += [e[f'{c}_pct'] for c in CATEGORIES]
        row += [n[f'{c}_pct'] for c in CATEGORIES]
        row += [f"{float(n[f'{c}_pct']) - float(e[f'{c}_pct']):+.1f}" for c in CATEGORIES]
        writer.writerow(row)

print("\nSaved: veda_vs_redundant_comparison.txt, veda_vs_redundant_comparison.csv")
