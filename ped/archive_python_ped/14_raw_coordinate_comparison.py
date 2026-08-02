"""
Step 14: Collect the RAW, per-individual-coordinate PED values from both
programs into one CSV -- no category grouping, no |value| conversion, no
degenerate-pair averaging. ped/13_compare_veda4.py's comparison summed both
sides into 5 categories first; this step is the un-collapsed version of the
same two sources, for anyone who wants to re-derive their own grouping
(or check this repo's category assignments) directly from the numbers each
program actually produced.

Sources:
  - Ours: PED_raw_pct.npy (Nint x 30, signed, %) -- step 4's per-coordinate
    PED before step 4 itself sums it into the 6 categories in
    PED_group_pct.npy. Nint=42 (redundant: 6 CC stretch, 6 CH stretch, 6
    CCC bend, 12 CCH bend, 6 ring torsion, 6 CH wag -- see
    PED_METHODOLOGY.md section 3), one row per coordinate, labeled by type
    and (0-based) atom indices from coord_labels.txt.
  - VEDA4: veda4/c6h6_reconstructed.ved contains TWO separate 30x30
    matrices, not one relabeled -- a "PED: sign = direction" table and a
    distinct "TED: sum = 100" table with different (not just sign-flipped
    or |.|'d) values, each with its own column-header line and 30 data
    rows. Both are extracted here as separate column blocks (veda4_PED_*
    and veda4_TED_*), one column per coordinate s1..s30, labeled by
    VEDA4's own STRE/BEND/TORS + CH/CC/CCC/HCC/HCCC/CCCC type tags from
    the "definitions of modes" section of c6h6_reconstructed.vdf. Nint=30
    -- non-redundant, see this repo's discussion of why VEDA4 has a
    smaller/differently-sized coordinate set per family than ours.

The two coordinate sets are NOT the same basis (different count, different
specific linear combinations within each motion-type family -- see this
repo's prior discussion) so there is no coordinate-to-coordinate
correspondence to align on. The only alignment applied here is matching
each of the 30 vibrational modes to its counterpart by RANK in descending
frequency (both lists sorted independently, then paired position-by-
position) -- because VEDA4 recomputes its own frequencies from the
imported Hessian and they differ slightly from Gaussian's originals (e.g.
3223.49 vs 3223.17 cm^-1), so matching on absolute frequency value isn't
reliable. No degenerate-pair merging is applied: all 30 modes get their
own row.

Requires: PED_raw_pct.npy, coord_labels.txt, vibfreq.npy (this repo's own
          pipeline) and, from the OUTSIDE-the-git-repo veda4/ directory,
          c6h6_reconstructed.ved + .vdf (only exist after VEDA4 has
          actually been run on 12_export_veda_fmt.py's output).
Output: raw_coordinate_comparison.csv
"""
import os
import re
import csv
import numpy as np

HERE = os.path.dirname(__file__)
VEDA_DIR = os.path.join(HERE, '..', '..', '..', '..', 'veda4')
VED_PATH = os.path.join(VEDA_DIR, 'c6h6_reconstructed.ved')
VDF_PATH = os.path.join(VEDA_DIR, 'c6h6_reconstructed.vdf')

for p in (VED_PATH, VDF_PATH):
    if not os.path.isfile(p):
        raise FileNotFoundError(
            f"{p} not found. This step requires VEDA4 to have already been "
            "run (via veda4e1.exe) on the .fmt file 12_export_veda_fmt.py "
            "produced -- it cannot run standalone.")


# --- Ours: raw per-coordinate PED (Step 4, un-grouped) -------------------

vibfreq = np.load(os.path.join(HERE, 'vibfreq.npy'))                    # (30,)
ped_raw = np.load(os.path.join(HERE, 'PED_raw_pct.npy'))                # (42, 30)
our_labels = []
for line in open(os.path.join(HERE, 'coord_labels.txt')):
    cat, atoms = line.strip().split('\t')
    idx0 = eval(atoms)                          # 0-based tuple, e.g. (0, 1)
    idx1 = '-'.join(str(a + 1) for a in idx0)   # 1-based, user-facing convention
    our_labels.append(f"ours_{cat}_{idx1}")


# --- VEDA4: raw per-coordinate PED/TED matrix -----------------------------

def parse_ved(path):
    """Returns (freqs, ped_matrix, ted_matrix) -- two distinct 30x30
    matrices, see module docstring."""
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

    return freqs, read_matrix(ped_hdr + 1), read_matrix(ted_hdr + 1)


def parse_coord_types(path):
    with open(path) as f:
        lines = f.readlines()
    start = next(i for i, l in enumerate(lines) if l.strip().startswith('definitions of modes'))
    pat = re.compile(r'^s\s*(\d+)\s+(\S+)\s+(\S+)')
    types = {}
    for l in lines[start + 1:]:
        m = pat.match(l.strip())
        if m:
            types[int(m.group(1))] = (m.group(2), m.group(3))
    n = max(types)
    return [types[i + 1] for i in range(n)]


veda_freqs, veda_ped, veda_ted = parse_ved(VED_PATH)  # (30,), (30,30), (30,30)
veda_types = parse_coord_types(VDF_PATH)               # len 30, [(kind, sub), ...]
veda_ped_labels = [f"veda4_PED_s{i+1}_{kind}_{sub}" for i, (kind, sub) in enumerate(veda_types)]
veda_ted_labels = [f"veda4_TED_s{i+1}_{kind}_{sub}" for i, (kind, sub) in enumerate(veda_types)]


# --- Align the 30 modes by rank in descending frequency only -------------

our_rank = np.argsort(-vibfreq)          # indices into vibfreq / ped_raw columns
veda_rank = np.argsort(-veda_freqs)      # indices into veda_freqs / veda_matrix rows

if len(our_rank) != len(veda_rank):
    raise ValueError(f"Mode count mismatch: {len(our_rank)} (ours) vs "
                      f"{len(veda_rank)} (VEDA4) -- cannot align by rank.")


# --- Write the combined, un-collapsed CSV ---------------------------------

csv_path = os.path.join(HERE, 'raw_coordinate_comparison.csv')
with open(csv_path, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['rank', 'freq_gaussian_cm-1', 'freq_veda4_cm-1'] +
               our_labels + veda_ped_labels + veda_ted_labels)
    for rank, (oi, vi) in enumerate(zip(our_rank, veda_rank), start=1):
        our_vals = [f"{x:.2f}" for x in ped_raw[:, oi]]
        veda_ped_vals = [f"{x:.2f}" for x in veda_ped[vi, :]]
        veda_ted_vals = [f"{x:.2f}" for x in veda_ted[vi, :]]
        w.writerow([rank, f"{vibfreq[oi]:.4f}", f"{veda_freqs[vi]:.4f}"] +
                   our_vals + veda_ped_vals + veda_ted_vals)

print(f"Wrote {csv_path}")
print(f"{len(our_rank)} rows (one per mode, unmerged), "
      f"{len(our_labels)} raw 'ours' coordinate columns + "
      f"{len(veda_ped_labels)} veda4_PED_* + {len(veda_ted_labels)} veda4_TED_* columns.")
print("\nNote: 'ours' columns are signed PED% summing to ~100% per row "
      "by construction (Pulay pseudoinverse normalization, step 4).")
print("veda4_PED_* columns are signed and generally do NOT sum to ~100% "
      "per row. veda4_TED_* columns are a SEPARATE table VEDA4 computes "
      "(not a transform of the PED one) whose signed values do sum to "
      "~100% per row directly, matching its own 'TED: sum=100' header.")
