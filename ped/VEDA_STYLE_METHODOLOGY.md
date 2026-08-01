# A VEDA-Style PED Cross-Check for Benzene

This document describes a **separate, additive** cross-check of the benzene
PED analysis (see `PED_METHODOLOGY.md` for the original, published pipeline),
built to test whether the redundant-coordinate/Pulay-pseudoinverse approach
already reported in `tab:benzenemixed` agrees with the *documented* algorithm
of VEDA (Jamróz's widely-used PED software). It does **not** modify `ped/01-05`,
their outputs, or the manuscript in any way — see "What this is not" below.

## 1. Why this exists

A prior session asked whether the existing PED pipeline aligns with VEDA.
Research into VEDA's published methodology (papers are paywalled; only
abstracts/summaries were accessible — see the References section) showed a
specific, meaningful difference: VEDA builds a **non-redundant** internal
coordinate set and **optimizes** it, rather than working directly with a
fixed, redundant set resolved by pseudoinverse. This document implements
that documented algorithm as a good-faith, from-the-literature
reconstruction — **not** a reproduction of VEDA4 itself, and it has never
been validated against actual VEDA4 output (VEDA4 was not run; only its
published description was used).

## 2. What is documented about VEDA vs. what is our own choice

**From VEDA's published description (citable):**
- Automatic coordinate generation of exactly (N−1) stretches + (N−2) bends +
  (N−3) torsions/out-of-plane = 3N−6, non-redundant.
- The objective function, quoted: **EPm = Σ_modes max|PED|** — the sum,
  over all vibrational modes, of that mode's single largest internal-
  coordinate PED contribution.
- An iterative "mixing" procedure: combine coordinates, keep the change if
  EPm increases.
- A final "reduction" pass eliminating components of multicomponent
  (mixed) coordinates that change EPm "in minimal degree."
- VEDA's own documentation explicitly flags that **rings break its naive
  automatic rule** and need special "ring coordinate" handling — i.e.
  benzene is a genuinely hard case for this method, not a trivial one, even
  for VEDA itself.

**Our own necessary implementation choices (not published by VEDA — see
Section 5 for exactly which choice was made and why):**
- The specific procedure used to resolve benzene's ring-closure redundancy
  into a non-redundant 30-coordinate set (Section 5.1).
- The move-generation strategy for mixing (same-family pairwise Givens
  rotations — Section 5.2), search order, and convergence criteria.
- The reduction-pass tolerance and elimination order (Section 5.3).
- The `abs()` convention inside EPm, and the largest-weight convention used
  to assign each optimized coordinate to one of the five reported PED
  categories for the final table (Section 5.4).

## 3. Method, step by step

Reuses `ped/01-03`'s already-validated output (`H_cart.npy`, `L.npy`,
`vibfreq.npy`, and the redundant `B.npy`/`coord_labels.txt`) as raw
material — nothing is re-parsed or re-derived from the Gaussian log.

### Step 6 — Natural (non-redundant) coordinates (`06_natural_coords.py`)
VEDA's own documentation flags that a naive (N−1)/(N−2)/(N−3) automatic
rule doesn't resolve a ring's "ring-closure" redundancy (going all the way
around a ring and back is a geometric constraint that isn't satisfied by
simply omitting one bond). The approach here follows the
symmetrized-linear-combination principle of natural/local-symmetry
coordinates (Pulay, Fogarasi, Pang & Boggs, *JACS* 1979, 101, 2550 — `RN12`
in `MolecularMotion.bib`; ring-redundancy treatment per the *spirit* of
Cremer & Pople, *JACS* 1975, 97, 1354 — their exact published coefficients
were not independently re-derivable here from paywalled text, so only the
general principle is used):

1. **Symmetrize** each of the 7 raw coordinate families (CC-stretch,
   CH-stretch, CCC-bend, CCH-bend split first into per-carbon
   symmetric/antisymmetric pairs, ring-torsion, CH-wag) via the real C6
   cyclic-group (Fourier) basis — a pure orthonormal basis rotation, same
   row space as the raw 42-coordinate `B.npy`, not a reduction yet
   (verified: rank before/after symmetrization is checked equal for every
   family).
2. **Rank-revealing selection** down to exactly 30 independent rows.
   A first attempt used a fixed priority order with an incremental
   per-stack SVD tolerance check; it reached 30 rows but left a
   badly-conditioned result (condition number ≈9×10⁷, several singular
   values sitting at ≈10⁻⁸ — a sign some accepted rows were only
   marginally independent). This was replaced with **max-residual-norm
   selection within each priority tier** (equivalent to modified
   Gram-Schmidt / pivoted QR, the standard robust rank-revealing
   procedure), keeping a chemically-motivated tier order (substituent-type
   coordinates first, ring-closure-suspect ones last) as a *tie-break*
   preference while letting the numerics — not an asserted-by-hand
   count — decide what's actually redundant. Result: condition number
   **21.1** (down from ≈9×10⁷), all 30 singular values in a healthy
   [0.31, 6.61] range.

**Honest, somewhat surprising finding:** with this tier ordering, the
algorithm found **zero independent CCC-bend directions** — every candidate
CCC-bend combination turned out to be a linear combination of already-kept
CC-stretch/ring-torsion directions. This is a legitimate result of
resolving the ring-closure redundancy, not a bug, but it is also **not a
unique physical truth**: which specific 12 of the 42 symmetrized
coordinates get dropped is a matter of basis choice when two candidates
compete for the same redundant subspace (here, CC-stretch and CCC-bend
candidates of matching symmetry species repeatedly competed, and the tier
ordering happened to favor keeping CC-stretch's representative each time).
A different, equally valid priority choice could instead retain CCC-bend
components and drop the corresponding CC-stretch harmonics. This has a
direct, visible consequence for Step 9's table: **the CCC-bend column is
0% for every mode**, by construction, in the tables produced here.

### Step 7 — Greedy EPm-maximizing mixing (`07_veda_mixing.py`)
Pairwise Givens rotations between two coordinate rows of the **same
family** (documented restriction: VEDA's own description is "superposition
of local modes" of compatible type; an unrestricted global search risks
producing EPm-maximizing but chemically uninterpretable composite
coordinates). For each of the 64 same-family candidate pairs, a 1D bounded
scalar optimization finds the rotation angle maximizing EPm; accepted if it
improves. Converged in 2 sweeps (10 accepted moves total, <1 second).
EPm rose from 24.82 (unmixed natural coordinates) to 27.75.

### Step 8 — Reduction pass (`08_veda_reduction.py`)
Sequentially tries dropping each mixed coordinate's smallest-magnitude
component (relative to the natural basis), accepting the drop if rank/
conditioning are preserved and EPm changes by less than `epm_tol=0.01`
(our own explicit choice, not VEDA's unpublished threshold — see Section 6
for the sensitivity check). Dropped 14 of the mixed coordinates' 50 total
components (mean nonzero components per final coordinate: 1.67 → 1.20);
EPm fell only slightly, to 27.74.

### Step 9 — Final table (`09_veda_ped_table.py`)
Reapplies the shared PED formula (`ped/_ped_core.py`, extracted from
`04_compute_ped.py` unchanged) to the final 30×36 coordinate matrix, groups
into the same five categories as `05_final_table.py`'s table (each final
coordinate assigned to whichever category carries the largest summed
squared weight in its composition — an explicit, documented convention),
and writes `veda_style_PED_table.csv`/`.txt` in the same schema as
`benzene_PED_table.csv`/`.txt`.

### Step 10 — Comparison report (`10_compare_report.py`)
Produces `veda_vs_redundant_comparison.txt`/`.csv`: a full per-mode,
per-category side-by-side table, EPm comparison, conditioning comparison,
and a factual checklist — see Section 4.

## 4. Results: how the two methods compare

**EPm** (VEDA's own purity objective, computed the same way — over
individual internal coordinates, not post-hoc categories — for both
methods, so the comparison is apples-to-apples):

| Method | EPm (max possible 30) |
|---|---|
| Existing, redundant 42-coordinate table | 6.08 |
| New, natural coordinates, unmixed | 24.82 |
| New, after mixing | 27.75 |
| New, after reduction (final) | 27.74 |

The new method's per-mode PED purity is dramatically higher — but **most
of that gain comes from resolving the ring-closure redundancy into a
well-conditioned, non-redundant basis** (unmixed EPm already 24.82), not
from the mixing/reduction optimization itself (+2.93 more). This is
itself informative: a properly-resolved non-redundant coordinate choice
matters far more than the optimization loop for this molecule.

**Agreement across all 30 modes:** mean absolute difference per
category-entry is 4.4 percentage points; 3 of 20 unique mode-lines show a
disagreement exceeding 30 percentage points in at least one category.

**The four manuscript-cited modes** (Table `tab:benzenemixed`):

| Freq (cm⁻¹) | Existing (published) | New (natural coordinates) | Max diff |
|---|---|---|---|
| 1056.39 | 44% CC stretch + 34% CCH bend + 22% CCC bend | **100% CCH bend** | 65.9 pts |
| 1319.27 | 64% CC stretch + 36% CCH bend | 64% CC stretch + 36% CCH bend | 0.0 pts |
| 1532.85 | 80% CCH bend + 13% CC stretch + 7% CCC bend | **100% CCH bend** | 19.4 pts |
| 1598.90 | 60% CC stretch + 33% CCH bend + 7% CCC bend | 60% CC stretch + 40% CCH bend | 6.8 pts |

**This is the central, decision-relevant finding.** Two of the four modes
agree closely (1319.27 exactly; 1598.90 within a few points once CCC-bend's
7% is accounted for as folded into CCH-bend). But the other two —
**1056.39 and 1532.85 cm⁻¹, which the manuscript's own discussion relies on
most heavily to argue the classifier correctly detects real mixed
stretch-bend (SB) character that a naive literature label misses** — show
the largest disagreement in the whole table: the new method reports them
as essentially pure bending (100% CCH bend, no stretch character at all),
where the existing table reports meaningful CC-stretch contributions (44%
and 13% respectively) that the manuscript's prose specifically cites as
evidence of stretch character.

This is a real numerical disagreement between two legitimate coordinate-
system choices, not a bug in either pipeline (see Section 5 for why).

## 5. Why the two methods disagree

Both are internally consistent, numerically validated (Section 6) PED
calculations. They differ because they are genuinely different coordinate
systems, and PED percentages are coordinate-system-dependent by
construction (this is a known, general property of PED analysis, not
specific to either implementation here):

- The **existing** table's 44%/13% CC-stretch contributions at these two
  modes come from the *redundant* 42-coordinate basis, where individual
  raw CC-stretch coordinates retain some independent identity even though
  they overlap with other coordinates.
- The **new** table's coordinate system was built by an algorithm that
  specifically *removed* CCC-bend as an independent direction (Section 3,
  Step 6) and favored high per-mode purity (Step 7's explicit EPm
  objective pushes toward exactly this kind of near-100%-single-category
  outcome). A method whose objective is "maximize how concentrated each
  mode's energy is in one coordinate" will systematically produce more
  monolithic-looking PED percentages than a fixed redundant basis with no
  such optimization — this is not unique to these two modes, but it is
  most visible here because these are the modes where the two coordinate
  systems' descriptions diverge most.

Neither number is "more correct" in an absolute sense without further
scrutiny — this is exactly why this is reported as a cross-check requiring
the user's own review, not resolved unilaterally here.

## 6. Numerical validation

| Check | Result |
|---|---|
| PED column sums, all stages (natural / mixed / final) | 1.0000–1.0008 for all 30 modes |
| `B_nat` condition number | 21.1 (down from ≈9×10⁷ before the max-residual-selection fix — see Step 6) |
| GF-matrix frequency cross-check (ω²=eig(G·Fq), G=B·M⁻¹·Bᵀ — independent of `PED_METHODOLOGY.md`'s Cartesian-Hessian round-trip) | max discrepancy 0.597 cm⁻¹, both for `B_nat` and `B_final` (looser than the original pipeline's 0.016 cm⁻¹ round-trip, as expected for a route involving an extra pseudoinverse/matrix-product chain, but still small relative to the 400–3200 cm⁻¹ spectrum) |
| `epm_tol` sensitivity (0.001 / 0.01 / 0.05) | EPm = 27.7495 / 27.7396 / 27.7396 — stable across the 0.01–0.05 range; 0.001 is more conservative (11 vs. 14 components dropped), not unstable |
| D6h degenerate-pair symmetry | All 10 true E-symmetry pairs have **identical purity to 4 decimal places** (max difference 0.0000). Each pair's two partners are assigned a *different* dominant coordinate label — this is expected degenerate-subspace basis-choice freedom (any orthonormal basis within a doubly-degenerate eigenspace is equally valid), not symmetry breaking by the mixing search. |
| Regression check | Confirmed zero changes to `ped/01-05`, their outputs, or `PED_METHODOLOGY.md` — `git diff --stat` before commit shows only new files added. |

## 7. What this is not

- **Not a reproduction of VEDA4.** VEDA4 was not run, decompiled, or
  otherwise inspected — only its published, peer-reviewed description was
  used (see References). Actual VEDA4 output on this exact system was
  never obtained for comparison.
- **Not a replacement for the published table.** `ped/01-05`,
  `benzene_PED_table.csv`/`.txt`, `PED_METHODOLOGY.md`, and
  `tab:benzenemixed` in the manuscript are all unchanged by this work.
  Promotion of these new numbers into the manuscript is an explicit,
  separate decision left to the user, informed by
  `veda_vs_redundant_comparison.txt`.
- **Not a claim that the new numbers are "more correct."** Section 5
  explains why the two methods disagree; neither is asserted here to be
  the superior description of these modes' physical character.

## References

- Jamróz, M. H. Vibrational Energy Distribution Analysis (VEDA): Scopes
  and limitations. *Spectrochim. Acta A* **2013**, *114*, 220–230.
  (Source of the documented EPm objective, automatic (N−1)/(N−2)/(N−3)
  coordinate generation, mixing/reduction procedure, and the ring-coordinate
  caveat. Full text was paywalled; abstract/summary content accessed via
  search engine result snippets and the developer's own software page,
  `smmg.pl/software/veda/`, not the original PDF.)
- Pulay, P.; Fogarasi, G.; Pang, F.; Boggs, J. E. Systematic ab initio
  gradient calculation of molecular geometries, force constants, and
  dipole moment derivatives. *J. Am. Chem. Soc.* **1979**, *101*,
  2550–2560. (`RN12` in `MolecularMotion.bib` — natural/local-symmetry
  coordinate principle used in Step 6.)
- Cremer, D.; Pople, J. A. A general definition of ring puckering
  coordinates. *J. Am. Chem. Soc.* **1975**, *97*, 1354–1358. (General
  principle for ring-redundancy reduction; not independently re-derived
  in detail here — see Step 6's caveat.)
- See `PED_METHODOLOGY.md` for the full reference list (Wilson, Pulay,
  Keresztury & Jalsovszky, Fraczkiewicz & Czernuszewicz, Jamróz 2014)
  underlying the original, published PED pipeline this document
  cross-checks against.
