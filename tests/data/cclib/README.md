# Real-program fixtures for the ORCA / Molden readers

Divinylbenzene at B3LYP/STO-3G, from the [cclib](https://github.com/cclib/cclib)
regression data (BSD 3-Clause, see `LICENSE`). Used by `tests/test_qc_parsers.py`
to check that the same molecule scores the same whichever program wrote it.

| File | Origin |
|---|---|
| `dvb_ir_g16.out` | `data/Gaussian/basicGaussian16/dvb_ir.out`, unchanged |
| `dvb_ir.molden` | `test/io/data/molden5.7_dvb_ir.molden`, unchanged |
| `dvb_ir_orca6.hess` | **repackaged** from `data/ORCA/basicORCA6.0/dvb_ir.out`: its `CARTESIAN COORDINATES (A.U.)`, `VIBRATIONAL FREQUENCIES` and `NORMAL MODES` blocks copied verbatim into `.hess` `$atoms` / `$vibrational_frequencies` / `$normal_modes` layout (irrep label lines dropped). cclib ships no ORCA `.hess`; the `NORMAL MODES` block is the same matrix ORCA writes to `$normal_modes`. No `$hessian` block. |
| `*.bonds` | generated once from covalent radii, one per file because the atom orders differ. The app itself never infers bonds. |

Separately (not vendored, licence unknown), the reader was checked against two
genuine ORCA `.hess` files from github.com/jmcneel1/ORCAHESS4TO3 (methane, and a
43-atom molecule with one imaginary mode): rediagonalising each file's own
`$hessian` reproduces its frequencies to 0.06 cm⁻¹ and the parsed Cartesian modes
to 1e-6 or better.
