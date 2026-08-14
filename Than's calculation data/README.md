# Than's calculation data — MP2/3-21G

Gaussian and VEDA4 results for the nine test molecules, rebuilt after the
basis-set change from **3-21G\*** to **3-21G**. Same two-category, point-group
layout as the previous 3-21G\* version of this folder.

```
gaussian/<PointGroup>/<Molecule>/   .com, .log, .fchk
ped/<PointGroup>/<Molecule>/
    input/      VEDA-ready .fmt
    veda_run/   full VEDA run directory (.ved .vdf .mpo .dd2 .log)
MANIFEST.csv    one row per molecule, including which basis it used
```

This folder holds **raw calculation data only** — no derived figures or score
tables. The scoring pipeline's default V-score weighting changed in `44295f7`
(`DEFAULT_V_WEIGHTING = "mu"`), so any figure or table is only meaningful
alongside the commit that produced it. The inputs here are independent of that
choice and can be re-scored under whichever weighting is current.

## The nine molecules

| Point group | Molecule | Basis | VEDA run |
|---|---|---|---|
| C2v | Acetone | 3-21G | `acetone_fchk` |
| C2v | H2O_C2v | 3-21G | `h2o_c2v_fchk` |
| C2v | oDifluorobenzene | 3-21G | `odifluorobenzene_fchk` |
| Cs | Cycloheptatriene | 3-21G | `cycloheptatriene_fchk` |
| Cs | FormicAcid | 3-21G | `formicacid_fchk` |
| Cs | HOCl | 3-21G | `hocl_fchk` |
| D2h | Ethylene | 3-21G | `ethylene_fchk` |
| D2h | Naphthalene | 3-21G | `naphthalene_fchk` |
| D2h | **XeF2Cl2** | **3-21G\*** | **`xef2cl2_bent2_fchk`** |

Every deck was checked against its own log: all nine agree on the basis they
actually ran.

## XeF2Cl2 — two things that make it the odd one out

### 1. It is still at 3-21G\*

Its 3-21G optimization **failed**: error termination via link 9999, 25 steps
exhausted with Maximum Force flat at ~0.100 against a 0.000450 threshold
(steps 17–25 read 0.1011, 0.0999, 0.0999, 0.1000, 0.1000, 0.1001, 0.1002,
0.1002, 0.1003 — a plateau, not slow convergence, so raising `maxcycles` would
not help). No frequencies were produced, so there was no Hessian to score or to
hand to VEDA. The d functions on the two Cl atoms appear to be what describes
the hypervalent Xe–Cl bonding; without them there is no D2h minimum to find.

Its files here are therefore regenerated from the surviving 3-21G\* run.
**Do not read the set as if the level of theory were uniform.**

Basis-function counts, for reference: the change is a genuine no-op for the
seven molecules built only from H, C, N, O and F, because `3-21G*` adds d
polarization on second-row atoms (Na–Ar) only.

| Molecule | 3-21G\* → 3-21G |
|---|---|
| Acetone, H2O_C2v, oDifluorobenzene, Cycloheptatriene, FormicAcid, Ethylene, Naphthalene | unchanged (48, 13, 80, 79, 31, 26, 106) |
| HOCl | 30 → 24 |
| XeF2Cl2 | 95 → 77 (run failed; 3-21G\* retained) |

Frequency-level effect: six molecules came back **bit-identical**; Acetone moved
in the last decimal of two modes (117.6668→117.6669, 170.4194→170.4195, i.e.
optimizer noise); only **HOCl** genuinely shifted — 687.19→667.64,
1291.14→1256.00, 3393.66→3368.54 cm⁻¹.

### 2. What `bent2` means and why it is needed

`ped/D2h/XeF2Cl2/veda_run/` is named `xef2cl2_bent2_fchk`, not `xef2cl2_fchk`.

VEDA4 cannot build internal coordinates through an exactly **180°** valence
angle: at 180° the bending plane is degenerate and the corresponding B-matrix
rows are singular, so VEDA's coordinate pre-optimization never converges.
Planar trans-XeF₂Cl₂ has **two** such angles — F2–Xe1–F3 and Cl4–Xe1–Cl5, both
exactly 180°.

`bent2` is the input produced by `ped_inputs/make_bent_fmt.py` with θ = 2°. It
rotates **F2 by +2°** and **Cl4 by −2°** about the axis normal to the molecular
plane, opening both trans angles to **178°** so VEDA can proceed.

What the distortion does and does not touch:

- The rotation is **within** the molecular plane. Every atom keeps `x = 0.000000`,
  so the molecule stays perfectly planar — it is *not* an out-of-plane pucker.
- Only **two atoms move**: F2 by 0.069 Å and Cl4 by 0.089 Å. Xe1, F3 and Cl5 do
  not move at all.
- **Bond lengths are unchanged**: Xe–F 1.974955 Å, Xe–Cl 2.540839 Å.
- Only the coordinate block of the `.fmt` differs. The frequency block and the
  entire Cartesian Hessian are copied through **byte for byte** from the
  undistorted file, so the frequencies are unaffected (matched to <0.03 cm⁻¹).

The cost: the B-matrix is built on a geometry 2° away from the one the Hessian
belongs to, so this molecule's %PED is one step less direct than the other
eight. A 5° variant (`bent5`) is kept in `input/` as a fallback if 2° is ever
too small for VEDA to converge. **The run reported here used bent2.**

## One thing to know before re-scoring

**Water here is `H2O_C2v`, which is not the library's `H2O`.** They are separate
calculations of the same molecule and differ by ~3 cm⁻¹ (3501.51 vs 3504.65,
3660.80 vs 3663.34). The VEDA %PED in this folder was run on `H2O_C2v`, but
`data/gjf` has no `H2O_C2v.com`, so `--classify` cannot build connectivity for
that geometry — any score-vs-PED join therefore reads scores from `H2O.log`
while reading %PED from `H2O_C2v`, and trips the pipeline's own 2.0 cm⁻¹ match
tolerance. This is long-standing rather than a consequence of the basis change;
the previous 3-21G\* data produced identical residuals. Adding an
`H2O_C2v.com` deck would let the join use a single consistent calculation.
