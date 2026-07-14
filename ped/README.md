# Benzene PED (Potential Energy Distribution) from the real JCC data

Reproduces a full Wilson-Decius-Cross / Pulay redundant-internal-coordinate
normal mode analysis of benzene, using the **actual Gaussian calculation
behind the JCC manuscript's benzene results**
(`data/logs/C6H6_MP2_3-21G_D6h.log`, MP2/3-21G, D6h, `freq=hpmodes`) — not
an independently-optimized geometry/Hessian from a separate quantum
chemistry run.

No literature PED table is used as input, and no VEDA4/`.chk`/`formchk`
detour is needed either: the Gaussian log never prints the literal "Force
constants in Cartesian coordinates" block VEDA4 requires, but that block
turns out to be unnecessary. Given all 3N-6 vibrational frequencies and
their full Cartesian displacement vectors (already in the log, parsed by
the repo's own `src/parser.py`), the full Cartesian Hessian is exactly
recoverable by mass-weighted eigendecomposition inversion.

## Setup

```bash
pip install numpy
```

No other dependencies — pure numpy, reusing `src/parser.py` and
`src/utils.py` from the parent `scoring-functions` repo (added to
`sys.path` at the top of step 1).

## Run

```bash
bash run_all.sh
```

or run the five scripts individually in order — each one prints its own
sanity checks and saves `.npy`/`.txt` files consumed by the next step:

| Script | What it does | Key output |
|---|---|---|
| `01_load_gaussian.py` | Parses `data/logs/C6H6_MP2_3-21G_D6h.log` via the repo's `GaussianParser` — geometry, all 30 HP-precision vibrational modes, connectivity | `opt_coords_ang.npy`, `l_raw.npy`, `vibfreq.npy`, `reduced_mass.npy`, `bonds.txt` |
| `02_reconstruct_hessian.py` | Reconstructs the full Cartesian Hessian directly from Gaussian's own frequencies + displacement vectors (no VEDA, no `.chk`); validates by re-diagonalizing and recovering Gaussian's own frequencies | `H_cart.npy`, `L.npy` |
| `03_build_internal_coords.py` | Detects the ring + substituents generically from the parsed bonds graph; defines 42 redundant internal coordinates (6 C-C str, 6 C-H str, 6 C-C-C bend, 12 C-C-H bend, 6 ring torsions, 6 C-H wags) and builds the Wilson B-matrix by finite differences | `B.npy` |
| `04_compute_ped.py` | F_q = (B⁺)ᵀ H B⁺ (Moore-Penrose pseudoinverse handles the redundancy); PED[n,μ] = D[n,μ]·(F_qD)[n,μ]/λ_μ | `PED_group_pct.npy` |
| `05_final_table.py` | Groups categories, merges the out-of-plane ambiguity, averages degenerate pairs, prints the final table | stdout |

## Things to check / verify yourself

- **Step 2's real correctness proof is the frequency round-trip**: it
  reconstructs the Hessian from Gaussian's own eigenvectors/eigenvalues,
  then independently re-diagonalizes it and confirms the same 30
  frequencies come back out (should agree to well under 1 cm⁻¹ — a much
  tighter bound than a from-scratch independent calculation could offer,
  since this is a mathematical round-trip, not two separate calculations).
  It also prints a mass-weighted-orthonormality residual for Gaussian's
  *as-printed* eigenvectors — a nonzero but small residual there (~1e-4)
  is expected and harmless, coming from finite print precision on
  degenerate E-symmetry mode pairs, not from an error in the
  reconstruction.
- **Step 4** prints the per-mode PED column sums — should all equal
  1.00. This is the actual correctness proof of the whole PED
  calculation, not a cosmetic check.
- **The out-of-plane ambiguity is real, not a bug.** For a planar
  hexagonal ring, "ring torsion" and "C-H out-of-plane wag" are not
  orthogonal internal coordinates, so their individual PED contributions
  can come out negative or >100% for out-of-plane modes. This is a
  known, published issue in benzene PED analysis (see Jamróz, "On the
  Internal Coordinates in the PED Analysis: Bending or Torsion?").
  Script 5 merges them into one "out-of-plane" category rather than
  reporting a misleading split.
- **No frequency scale factor is applied.** There is no verified,
  citable published scale factor for MP2/3-21G (unlike, e.g., Scott &
  Radom's well-known 0.8929 for HF/6-31G(d)) — asserting one without a
  verified source would repeat a citation mistake from earlier in this
  project's development. PED percentages don't depend on scaling anyway;
  only the displayed cm⁻¹ column would change. The reported frequencies
  are Gaussian's own computed (unscaled) harmonic values.
- Mass units (amu vs. kg) and the absolute scale of λ are both provably
  irrelevant to the final PED ratio (cancels out algebraically as long
  as used consistently) — so everything here stays in Gaussian's native
  amu/Å units throughout, no SI unit conversion needed or performed.

## Things you can change to see the numbers respond

- **Internal coordinate definitions**: script 03 uses a *redundant*
  set + pseudoinverse. You could instead build a strictly
  non-redundant 30-coordinate set (removing the linear dependencies at
  each sp² carbon by hand) and compare — the PED percentages should be
  very close but not always identical, since redundant vs.
  non-redundant PED definitions aren't mathematically forced to agree
  exactly for coupled/degenerate modes.
- **Level of theory**: to compare against a different Gaussian
  calculation, point `LOG_PATH` in `01_load_gaussian.py` at a different
  `freq=hpmodes` log (same atom-ordering/connectivity assumptions in
  step 3 apply to any single-ring, one-H-per-carbon benzene job).
