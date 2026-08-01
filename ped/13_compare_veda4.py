"""
Step 13: Combine this pipeline's own PED results (benzene_PED_table.csv,
redundant 42-coordinate, Pulay pseudoinverse) with the REAL VEDA4.exe
results (veda4/c6h6_reconstructed.ved, produced by running VEDA4 on the
.fmt file 12_export_veda_fmt.py built) into a single comparison table.

This is a different comparison from ped/10_compare_report.py, which only
compares this repo's own two Python pipelines (redundant vs the from-
scratch VEDA-style reimplementation, ped/06-11_*.py) against each other --
neither side of that one is the real VEDA4 program.

Category mapping. VEDA4's internal coordinates (defined in the .fmt-derived
.dd2, type-labeled in the "definitions of modes" section of .vdf) map
one-to-one onto this pipeline's 5 reporting categories:
  STRE CH        -> CH_stretch
  STRE CC        -> CC_stretch
  BEND CCC       -> CCC_bend
  BEND HCC       -> CCH_bend
  TORS HCCC/CCCC -> out-of-plane (matches this pipeline's own merge of
                     ring-torsion + CH-wag into one category, see
                     PED_METHODOLOGY.md section 4)

Sign convention. veda4/c6h6_reconstructed.ved actually contains TWO
separate 30x30 matrices, not one relabeled: a "PED: sign = direction"
table first, then a distinct "TED: sum = 100" table with different
(not just sign-flipped or |.|'d) values. The TED table's own signed
values already sum to ~100 per mode directly -- no absolute-value
post-processing needed or applied here; this script uses the TED table.

Degenerate E-symmetry pairs are matched between the two programs by rank
position in the descending-frequency list after grouping consecutive
modes with a gap <= 0.5 cm^-1 (the same threshold and grouping logic
05_final_table.py already validated against this exact spectrum) --
NOT by comparing absolute frequency values, since VEDA4 recomputes its
own frequencies from the imported Hessian and they differ slightly
(e.g. 3223.49 vs Gaussian's original 3223.17) even for the same mode.

Requires: vibfreq.npy, PED_group_pct.npy, cats.txt (01, 05, this repo's
          own pipeline) and, from the OUTSIDE-the-git-repo veda4/
          directory, c6h6_reconstructed.ved + .vdf (only exist after the
          user has actually run VEDA4 on 12_export_veda_fmt.py's output --
          this step cannot run standalone).
Output: veda4_comparison.csv, veda4_comparison.txt
"""
import os
import re
import csv
import numpy as np

HERE = os.path.dirname(__file__)
VEDA_DIR = os.path.join(HERE, '..', '..', '..', 'veda4')
VED_PATH = os.path.join(VEDA_DIR, 'c6h6_reconstructed.ved')
VDF_PATH = os.path.join(VEDA_DIR, 'c6h6_reconstructed.vdf')

CAT_ORDER = ['CH_stretch', 'CC_stretch', 'CCC_bend', 'CCH_bend', 'OOP']
CAT_HEADER = ['CHstr', 'CCstr', 'CCCbend', 'CCHbend', 'OOP']

for p in (VED_PATH, VDF_PATH):
    if not os.path.isfile(p):
        raise FileNotFoundError(
            f"{p} not found. This step requires VEDA4 to have already been "
            "run (via veda4e1.exe) on the .fmt file 12_export_veda_fmt.py "
            "produced -- it cannot run standalone.")


# --- Our own pipeline's PED (redundant 42-coordinate, Step 4/5) ---------

vibfreq = np.load(os.path.join(HERE, 'vibfreq.npy'))          # (30,) Gaussian original, raw parse order
pg = np.load(os.path.join(HERE, 'PED_group_pct.npy'))          # (6, 30) same raw order
raw_cats = [l.strip() for l in open(os.path.join(HERE, 'cats.txt'))]
raw_idx = {c: i for i, c in enumerate(raw_cats)}

our_pcts = np.vstack([
    pg[raw_idx['CH_stretch']],
    pg[raw_idx['CC_stretch']],
    pg[raw_idx['CCC_bend']],
    pg[raw_idx['CCH_bend']],
    pg[raw_idx['ring_torsion']] + pg[raw_idx['CH_wag_oop']],
])  # (5, 30)


# --- Parse VEDA4's real output -------------------------------------------

def parse_ved(path):
    """Returns (freqs, ped_matrix, ted_matrix). .ved contains TWO distinct
    30x30 matrices -- a signed 'PED: sign = direction' table, then a
    separately-computed 'TED: sum = 100' table (different values, not a
    transform of the first) -- each preceded by its own '1 2 3 ... 30'
    column-header line and followed by 30 data rows."""
    with open(path) as f:
        lines = f.readlines()

    dfac_line = next(i for i, l in enumerate(lines) if 'diagonality factor' in l)
    ped_hdr = next(i for i, l in enumerate(lines) if l.strip().startswith('PED: sign'))
    ted_hdr = next(i for i, l in enumerate(lines) if l.strip().startswith('TED:'))

    freq_lines = [l for l in lines[dfac_line + 1:ped_hdr] if l.strip()]
    freqs = np.array([float(x) for l in freq_lines for x in l.split()])

    def read_matrix(header_line_idx):
        rows = []
        i = header_line_idx + 1
        while len(rows) < len(freqs):
            nums = [float(x) for x in lines[i].split()]
            row_label, values = int(nums[0]), nums[1:-1]
            if len(values) != len(freqs):
                raise ValueError(f"{path}: row {row_label} has {len(values)} "
                                  f"coordinate values, expected {len(freqs)}.")
            rows.append(values)
            i += 1
        return np.array(rows)

    ped_matrix = read_matrix(ped_hdr + 1)
    ted_matrix = read_matrix(ted_hdr + 1)
    return freqs, ped_matrix, ted_matrix  # (30,), (30,30), (30,30)


def parse_coord_categories(path):
    with open(path) as f:
        lines = f.readlines()
    start = next(i for i, l in enumerate(lines) if l.strip().startswith('definitions of modes'))
    pat = re.compile(r'^s\s*(\d+)\s+(\S+)\s+(\S+)')
    coord_type = {}
    for l in lines[start + 1:]:
        m = pat.match(l.strip())
        if not m:
            continue
        coord_type[int(m.group(1))] = (m.group(2), m.group(3))

    def category(kind, sub):
        if kind == 'STRE' and sub == 'CH':
            return 'CH_stretch'
        if kind == 'STRE' and sub == 'CC':
            return 'CC_stretch'
        if kind == 'BEND' and sub == 'CCC':
            return 'CCC_bend'
        if kind == 'BEND' and sub in ('HCC', 'CCH'):
            return 'CCH_bend'
        if kind in ('TORS', 'OUT'):
            return 'OOP'
        raise ValueError(f"Unmapped VEDA4 coordinate type/subtype: {kind} {sub}")

    n = max(coord_type)
    return [category(*coord_type[i + 1]) for i in range(n)]


veda_freqs, _ped_matrix, veda_matrix = parse_ved(VED_PATH)  # use the TED table
coord_cats = parse_coord_categories(VDF_PATH)                # len 30

veda_pcts = np.zeros((len(CAT_ORDER), veda_matrix.shape[0]))
for ci, cat in enumerate(CAT_ORDER):
    cols = [j for j, c in enumerate(coord_cats) if c == cat]
    veda_pcts[ci] = veda_matrix[:, cols].sum(axis=1)  # signed -- TED already sums to ~100


# --- Align modes by rank within degenerate groups, not raw frequency ----

def group_degenerate(freqs, tol=0.5):
    order = np.argsort(-freqs)
    groups, i = [], 0
    while i < len(order):
        group = [order[i]]
        while i + 1 < len(order) and abs(freqs[order[i]] - freqs[order[i + 1]]) <= tol:
            i += 1
            group.append(order[i])
        groups.append(group)
        i += 1
    return groups


our_groups = group_degenerate(vibfreq)
veda_groups = group_degenerate(veda_freqs)

if len(our_groups) != 20 or len(veda_groups) != 20:
    raise ValueError(f"Expected 20 degeneracy groups each (benzene's known "
                      f"pattern); got {len(our_groups)} (ours), "
                      f"{len(veda_groups)} (VEDA4). Refusing to align by rank.")
if [len(g) for g in our_groups] != [len(g) for g in veda_groups]:
    raise ValueError("Degeneracy pattern (group sizes) differs between our "
                      "list and VEDA4's -- rank-based alignment is unsafe.")


# --- Assemble comparison rows, ascending frequency (matches existing table) -

rows = []
for og, vg in zip(our_groups, veda_groups):
    rows.append({
        'deg': len(og),
        'freq_ours': vibfreq[og].mean(),
        'freq_veda4': veda_freqs[vg].mean(),
        'ours': our_pcts[:, og].mean(axis=1),
        'veda4': veda_pcts[:, vg].mean(axis=1),
    })
rows.sort(key=lambda r: r['freq_ours'])

csv_path = os.path.join(HERE, 'veda4_comparison.csv')
with open(csv_path, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(
        ['freq_gaussian_cm-1', 'degeneracy'] +
        [f'ours_{c}_pct' for c in CAT_HEADER] +
        ['freq_veda4_cm-1'] +
        [f'veda4_{c}_pct' for c in CAT_HEADER]
    )
    for r in rows:
        w.writerow(
            [f"{r['freq_ours']:.1f}", r['deg']] +
            [f"{x:.1f}" for x in r['ours']] +
            [f"{r['freq_veda4']:.2f}"] +
            [f"{x:.1f}" for x in r['veda4']]
        )

txt_path = os.path.join(HERE, 'veda4_comparison.txt')
with open(txt_path, 'w') as f:
    header = (f"{'freq(G)':>9} {'deg':>3} | "
              f"{'ours:CHs':>8}{'CCs':>6}{'CCCb':>6}{'CCHb':>6}{'OOP':>6} | "
              f"{'freq(V4)':>9} | "
              f"{'v4:CHs':>7}{'CCs':>6}{'CCCb':>6}{'CCHb':>6}{'OOP':>6}")
    f.write(header + "\n")
    f.write("-" * len(header) + "\n")
    for r in rows:
        o, v = r['ours'], r['veda4']
        f.write(f"{r['freq_ours']:9.1f} {r['deg']:>3} | "
                f"{o[0]:8.1f}{o[1]:6.1f}{o[2]:6.1f}{o[3]:6.1f}{o[4]:6.1f} | "
                f"{r['freq_veda4']:9.2f} | "
                f"{v[0]:7.1f}{v[1]:6.1f}{v[2]:6.1f}{v[3]:6.1f}{v[4]:6.1f}\n")
    f.write(f"\n{len(rows)} unique frequency rows (30 modes, degenerate "
            "E-symmetry pairs averaged and shown once, deg=2).\n")
    f.write("'ours' = ped/benzene_PED_table.csv (redundant 42-coordinate, "
            "signed PED, this repo's own from-scratch reconstruction).\n")
    f.write("'v4' = real VEDA4.exe output (veda4/c6h6_reconstructed.ved's "
            "'TED: sum=100' table, categories summed signed, no |.| -- "
            "see this script's docstring).\n")

print(f"Wrote {csv_path}")
print(f"Wrote {txt_path}")
print(f"\n{len(rows)} comparison rows written.")
max_gap = max(abs(r['freq_ours'] - r['freq_veda4']) for r in rows)
print(f"Max |Gaussian - VEDA4| recomputed frequency gap: {max_gap:.2f} cm^-1")
