# Benzene PED Analysis: Methodology, Theoretical Basis, and Numerical Validation

This document is the supporting technical record for the Potential Energy
Distribution (PED) analysis implemented in `ped/` and reported in the JCC
manuscript ("A Unified, Reference-Free Framework for Classifying the 3N
Modes of Molecular Motion"), Table `tab:benzenemixed`. It exists to answer,
in one place, three questions: what the code computes, why each step is
mathematically correct, and how the result was numerically validated.

The condensed version of this document, sized for the manuscript's
Supporting Information bundle, is `JCC/JCC_man_scoring/JCC_SI_PED_benzene.tex`.

A separate, additive cross-check of this pipeline against VEDA's published
(not reverse-engineered) methodology — a non-redundant, EPm-optimized
coordinate approach — lives in `VEDA_STYLE_METHODOLOGY.md`. It does not
modify anything described in this document; see that file for a comparison
report, including a real disagreement at two of the four modes discussed in
Section 6 below.

## 1. Purpose

The main-text classification framework in this repository scores each
normal mode's translational, rotational, and stretch/bend vibrational
character directly from Cartesian mode displacements, without reference to
internal coordinates or a force field. Table `tab:benzenemixed` in the
manuscript checks four of that framework's mixed-character benzene modes
(1056.39, 1319.27, 1532.85, and 1598.90 cm$^{-1}$) against an independent,
internal-coordinate-based method: the classical Potential Energy
Distribution (PED) analysis of Wilson, Decius, and Cross, extended to
redundant internal coordinates by Pulay and co-workers. The `ped/` pipeline
computes that independent PED reference **from scratch**, from the same
Gaussian calculation used everywhere else in the manuscript — it is not a
literature PED table. (The manuscript's reference *mode labels/numbering*
for those four rows, e.g. "$\nu_{14}$", are separately cited to
Shimanouchi's tables; the PED percentages themselves are this code's own
output.)

## 2. Data provenance

- **Source file**: `data/logs/C6H6.log` — a Gaussian MP2/3-21G calculation
  of benzene at its $D_{6h}$-symmetric equilibrium geometry, run with
  `freq=hpmodes` (high-precision mode printing). This is the exact same
  calculation behind every other benzene figure and table in the
  manuscript — not an independently re-optimized geometry or a separate
  quantum-chemistry run.
- **Connectivity**: `data/gjf/C6H6.com`, matched to the log by filename
  basename (the repository's convention — see the top-level `CLAUDE.md`).
- **Parsing**: both files are read via the repository's own
  `src/parser.py` (`GaussianParser`), the same parser used by the main
  scoring pipeline — so the PED analysis and the main framework's scores
  are computed from identical, independently-verifiable input data.

No literature PED table, VEDA/VEDA4 output, or `.chk`/`formchk` file is
used as input anywhere in this pipeline.

## 3. Method, step by step

The pipeline is five sequential scripts (`ped/01_*.py` through
`ped/05_*.py`), run via `ped/run_all.sh`. Each step's formula and its
theoretical grounding are given below in the order they run.

### Step 1 — Load geometry, normal modes, and connectivity
`ped/01_load_gaussian.py` parses `data/logs/C6H6.log`: the equilibrium
Cartesian geometry, all 30 ($=3N-6$ for $N=12$ atoms) vibrational
frequencies $\nu_k$, their raw Cartesian displacement vectors
$\boldsymbol{l}_k$ as printed by Gaussian, the reduced mass $\mu_k$ Gaussian
reports for each mode, and the bond connectivity.

### Step 2 — Reconstruct the Cartesian force-constant (Hessian) matrix
`ped/02_reconstruct_hessian.py` recovers the full $36\times36$ Cartesian
Hessian $\mathbf{H}$ directly from Gaussian's own printed eigendata,
without needing a "Force constants in Cartesian coordinates" block (which
this particular log does not print) or any external tool.

**Derivation.** The standard mass-weighted vibrational eigenproblem (Wilson,
Decius & Cross, 1980, ch. 2)[^Wil1980] is
$$
\mathbf{M}^{-1/2}\mathbf{H}\,\mathbf{M}^{-1/2} = \mathbf{L}\,\boldsymbol{\Lambda}\,\mathbf{L}^{\mathrm T},
\qquad \mathbf{L}^{\mathrm T}\mathbf{L}=\mathbf{I},
$$
with $\boldsymbol{\Lambda}=\mathrm{diag}(\omega_k^2)$,
$\omega_k = 2\pi c\,\tilde\nu_k$ ($c$ in cm/s, $\tilde\nu_k$ in cm$^{-1}$, the
standard, exact conversion). Equivalently, in terms of Cartesian
(non-mass-weighted) mode vectors $\mathbf{L}_{\text{cart}}$ satisfying
$\mathbf{L}_{\text{cart}}^{\mathrm T}\mathbf{M}\mathbf{L}_{\text{cart}}=\mathbf{I}$,
$$
\mathbf{H} = \mathbf{M}\,\mathbf{L}_{\text{cart}}\,\boldsymbol{\Lambda}\,\mathbf{L}_{\text{cart}}^{\mathrm T}\,\mathbf{M}.
$$
This is the standard eigenproblem *run in reverse*: given the eigenvectors
and eigenvalues (which Gaussian already computed and printed), $\mathbf{H}$
is recoverable exactly, with no separate force-constant calculation needed.
This reconstruction direction is not spelled out explicitly in Wilson–Decius–Cross,
so it is presented here as an algebraic consequence of their eigenproblem,
not as a separate literature-attributed result.

Gaussian's printed $\boldsymbol{l}_k$ is proportional to the correct
Cartesian eigenvector direction, but its normalization is Gaussian's own
convention, not $\mathbf{L}_{\text{cart}}^{\mathrm T}\mathbf{M}\mathbf{L}_{\text{cart}}=\mathbf{I}$.
The code rescales directly rather than assuming a specific normalization
convention:
$$
\mathbf{L}_k = \boldsymbol{l}_k \,\big/\, \sqrt{\textstyle\sum_i m_i\,l_{k,i}^2},
$$
which satisfies $\mathbf{L}^{\mathrm T}\mathbf{M}\mathbf{L}=\mathbf{I}$ by
construction regardless of Gaussian's internal convention. Mass units
(amu vs. kg) and the absolute scale of $\lambda_k$ are both algebraically
irrelevant to the final PED ratio (Step 4), so plain amu is used
throughout with no unit conversion.

**Correctness proof used.** Rather than trust the derivation alone, the
script re-diagonalizes the reconstructed $\mathbf{H}$ and confirms the
same 30 frequencies come back out — a mathematical round-trip, not an
independent calculation, but a much tighter and more direct test than a
second from-scratch quantum-chemistry run could offer. See §5 for the
numeric result.

### Step 3 — Redundant internal coordinates and the Wilson B-matrix
`ped/03_build_internal_coords.py` auto-detects the six-membered ring and
its substituents from the parsed bond graph (not a hardcoded atom order),
and defines 42 internal coordinates for benzene's 30 vibrational degrees
of freedom — deliberately redundant, in the sense of Pulay & Török (1966)[^Pulay1966]:

| Coordinate type | Count | Definition |
|---|---|---|
| C–C stretch | 6 | ring bond lengths |
| C–H stretch | 6 | substituent bond lengths |
| C–C–C bend | 6 | in-plane ring angles |
| C–C–H bend | 12 | 2 per ring carbon (both ring neighbors) |
| Ring torsion | 6 | C–C–C–C dihedrals around the ring |
| C–H out-of-plane wag | 6 | H–C–C–C dihedrals |

The Wilson $\mathbf{B}$-matrix, $B_{n,i}=\partial q_n/\partial x_i$
(internal coordinate $n$ with respect to Cartesian displacement $i$), is
built by numerical central differences at a $\pm10^{-4}$ Å step, with
dihedral-angle wraparound handled explicitly. This is a standard,
literature-established alternative to deriving analytical $B$-matrix
elements by hand for each coordinate type (Wilson, Decius & Cross, 1980).[^Wil1980]

### Step 4 — Force-constant transformation and PED
`ped/04_compute_ped.py` implements Pulay's generalized-inverse treatment of
redundant internal coordinates.[^Pulay1966] At a stationary point,
$\mathbf{H}=\mathbf{B}^{\mathrm T}\mathbf{F}_q\mathbf{B}$ exactly for a
non-redundant, square, invertible $\mathbf{B}$; for the redundant,
non-square $\mathbf{B}$ used here, the Moore–Penrose pseudoinverse
$\mathbf{B}^{+}$ replaces the ordinary inverse:
$$
\mathbf{F}_q = (\mathbf{B}^{+})^{\mathrm T}\,\mathbf{H}\,\mathbf{B}^{+}.
$$
For each normal mode $\mu$ (Cartesian eigenvector $\mathbf{L}_\mu$,
eigenvalue $\lambda_\mu=\omega_\mu^2$), the internal-coordinate projection
is $\mathbf{D}_\mu=\mathbf{B}\mathbf{L}_\mu$, and the PED of internal
coordinate $n$ in mode $\mu$ is
$$
\mathrm{PED}_{n,\mu} = \frac{D_{n,\mu}\,\big(\mathbf{F}_q\mathbf{D}_\mu\big)_n}{\lambda_\mu},
$$
the standard PED decomposition of Keresztury & Jalsovszky (1971)[^Keresztury1971]
and Fraczkiewicz & Czernuszewicz (1997),[^Fraczkiewicz1997] normalized so
$\sum_n \mathrm{PED}_{n,\mu}=1$ for every mode. That column-sum-to-one
property is the correctness criterion checked numerically at every run
(§5) — it is a genuine test of the full transformation, not a cosmetic
check, since an error anywhere in $\mathbf{B}$, $\mathbf{H}$, or the
pseudoinverse would break it.

### Step 5 — Category grouping and reporting
`ped/05_final_table.py` sums the raw per-coordinate PED into six physically
meaningful categories (CH stretch, CC stretch, CCC bend, CCH bend, and a
merged out-of-plane category — see §6), averages numerically-degenerate
$E$-symmetry mode pairs (frequency gap $\le0.5$ cm$^{-1}$, validated against
the actual spectrum in §5), and writes `benzene_PED_table.csv`/`.txt`.

## 4. Known, documented limitations

- **Ring-torsion / C–H-wag non-orthogonality.** For a planar hexagonal
  ring, "ring torsion" and "C–H out-of-plane wag" are not orthogonal
  internal coordinates, so their *individual* PED contributions can come
  out negative or exceed 100% for out-of-plane modes. This is a real,
  documented ambiguity in redundant-coordinate PED analysis, not a bug —
  see Jamróz (2014).[^Jamroz2014] Step 5 merges the two into a single
  "out-of-plane" category rather than reporting a misleading split; the
  merged category's PED is well-behaved (0–100%, as shown in §5's
  numbers).
- **No frequency scale factor applied.** No verified, citable harmonic
  scale factor exists for MP2/3-21G (unlike, e.g., Scott & Radom's
  well-established 0.8929 for HF/6-31G(d)); asserting one without a
  verified source would repeat a citation error from earlier in this
  project's development. This is a deliberate choice, not an oversight —
  PED percentages are scale-invariant regardless (they depend on
  $D_{n,\mu}(\mathbf{F}_q\mathbf{D}_\mu)_n/\lambda_\mu$ ratios, not on the
  absolute frequency), so only the displayed cm$^{-1}$ column would move
  if a factor were applied. Reported frequencies are Gaussian's own
  unscaled harmonic values throughout.
- **Redundant vs. non-redundant coordinate choice.** A strictly
  non-redundant 30-coordinate set (removing linear dependencies by hand at
  each sp$^2$ carbon) would give PED percentages that are close but not
  guaranteed identical for coupled/degenerate modes — redundant and
  non-redundant PED definitions are not mathematically forced to agree
  exactly. The redundant formulation is used here because it requires no
  manual dependency removal and is the more general, standard treatment
  (Pulay & Török, 1966).[^Pulay1966]
- **A previously-fixed bug, noted for the record.** The reduced-mass
  cross-check formula in `02_reconstruct_hessian.py` was inverted in an
  earlier version of this pipeline, producing a spurious ~35% mean
  "error" against a diagnostic that was itself correct; this was caught
  and fixed in commit `22f38ca`. It is noted here as evidence that the
  numbers in this document were arrived at with actual scrutiny, not
  assumed correct on the first pass.

## 5. Numerical validation

All figures below are from a fresh run of `bash ped/run_all.sh` against
`data/logs/C6H6.log`, plus supplementary diagnostic checks described where
noted. Full diagnostic output is reproducible by re-running the pipeline.

| Check | Result | Interpretation |
|---|---|---|
| Frequency round-trip (Step 2: reconstruct $\mathbf{H}$, re-diagonalize, compare to Gaussian's original 30 frequencies) | max discrepancy **0.0163 cm$^{-1}$** | Well under the 2 cm$^{-1}$ pass threshold; confirms the Hessian reconstruction (§3, Step 2) is correct to near machine/print precision. |
| Mass-weighted orthonormality of Gaussian's *as-printed* eigenvectors, $\max\lvert\mathbf{L}^{\mathrm T}\mathbf{M}\mathbf{L}-\mathbf{I}\rvert$ | **3.66$\times10^{-4}$** | Small residual from finite print precision on degenerate $E$-symmetry pairs, not from the reconstruction; the frequency round-trip above is the authoritative check. |
| Reduced-mass cross-check vs. Gaussian's own printed $\mu_k$ (informational; not required for correctness) | max rel. error **8.55$\times10^{-4}$** (mean 4.48$\times10^{-4}$) | Confirms the displacement-normalization convention assumed in Step 2 independently, to <0.1%. |
| T/R null-space eigenvalues of the reconstructed, re-diagonalized $\mathbf{H}$ | $\sim10^{13}$–$10^{14}$ (rad/s)$^2$, vs. $\sim10^{28}$ for real vibrational modes | ~15 orders of magnitude smaller than genuine vibrational eigenvalues — floating-point noise around an exact zero, not a defect; translations/rotations are never included in the reconstruction. |
| PED column sums, all 30 modes (Step 4's own correctness check) | **1.000–1.001** for every mode | Confirms $\mathbf{B}$, $\mathbf{H}$, and the pseudoinverse transformation are jointly self-consistent for every mode, not just the four reported in the manuscript. |
| B-matrix singular-value spectrum (42$\times$36, at production step size $10^{-4}$ Å) | 30 singular values in $[0.618,\,6.61]$; 6 near-zero, $\le3.40\times10^{-9}$; gap ratio $\approx1.8\times10^{8}$ | Confirms `rcond=1e-8` in `np.linalg.pinv` sits in a clean 8-order-of-magnitude gap, cleanly separating the 30 true internal-coordinate directions from the 6 null directions corresponding to overall translation/rotation (to which every internal coordinate is exactly invariant). Not an arbitrary tolerance choice. |
| Finite-difference step-size sensitivity ($10^{-5}$, $10^{-4}$, $10^{-3}$ Å, three orders of magnitude) | PED percentages for the four manuscript-reported modes change by $<0.001$ percentage points across all three step sizes | Confirms the numerical $\mathbf{B}$-matrix is not sitting in a truncation- or cancellation-error-sensitive regime; $10^{-4}$ Å is a safe, unremarkable choice. |
| Degenerate-mode merge threshold (0.5 cm$^{-1}$) vs. actual spectrum | Exactly 10 frequency gaps $\le0.5$ cm$^{-1}$ (all exactly 0.000000 cm$^{-1}$ as computed), matching benzene's 10 doubly-degenerate $E$-symmetry pairs; smallest *non*-degenerate gap is 2.28 cm$^{-1}$ | 0.5 cm$^{-1}$ sits in a wide, unambiguous gap — no risk of over- or under-merging for this spectrum. |
| Regression check: pipeline output before vs. after the `data/logs/C6H6_MP2_3-21G_D6h.log` → `C6H6.log` filename fix (§7) | **Byte-identical**, except one corrected filename string in the `.txt` footer | Confirms the July 2026 file rename did not silently change the underlying Gaussian data; the committed `benzene_PED_table.csv`/`.txt` from before the rename remain numerically valid. |
| Cross-check against manuscript Table `tab:benzenemixed` (4 rows: 1056.39, 1319.27, 1532.85, 1598.90 cm$^{-1}$) | **Exact match**, all four rows, all category percentages | Confirms this pipeline is in fact the source of the manuscript's reported PED percentages. |

## 6. Relationship to the manuscript

`JCC/JCC_man_scoring/JCC_man_CA.tex`, Table `tab:benzenemixed`
(§"Applications to benzene modes of motion"), reports this pipeline's PED
percentages for the four benzene modes whose classification the main
reference-free framework disagrees with a naive literature label, using
them as an independent check that the framework's disagreement reflects a
genuine mixed stretch/bend character rather than a classifier error. The
"Ref. assignment" column's mode numbering (e.g. "$\nu_{14}$") is separately
cited to Shimanouchi (1972);[^Shi1972] the "PED description" column is this
pipeline's own output, reproduced exactly as validated in §5.

## References

[^Wil1980]: Wilson, E. B.; Decius, J. C.; Cross, P. C. *Molecular
    Vibrations: The Theory of Infrared and Raman Vibrational Spectra*;
    Dover: New York, 1980.
[^Pulay1966]: Pulay, P.; Török, F. On the calculation of normal
    coordinates. *Acta Chim. Acad. Sci. Hung.* **1966**, *47*, 273.
[^Keresztury1971]: Keresztury, G.; Jalsovszky, G. An alternative
    calculation of the vibrational potential energy distribution.
    *J. Mol. Struct.* **1971**, *10*, 304.
[^Fraczkiewicz1997]: Fraczkiewicz, R.; Czernuszewicz, R. S. The internal
    coordinate normal mode analysis: vibrational potential energy
    distributions. *J. Mol. Struct.* **1997**, *435*, 109–113.
[^Jamroz2014]: Jamróz, M. H. On the Internal Coordinates in the Potential
    Energy Distribution (PED) Analysis: Bending or Torsion?
    *Enliven: Bioinformatics* **2014**, *1*(4), 006.
    doi:10.18650/2376-9416.14006.
[^Shi1972]: Shimanouchi, T. *Tables of Molecular Vibrational Frequencies,
    Consolidated Volume I*; NSRDS-NBS 39; U.S. National Bureau of
    Standards, 1972.
