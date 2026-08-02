# Archive note (2026-08-02)

This directory holds the original from-scratch Python PED pipeline (`01`-`05`, a
redundant-internal-coordinate Wilson-Decius-Cross/Pulay calculation), the abandoned
"VEDA-style" reimplementation (`06`-`11`, never promoted), and the real-VEDA4 bridge
scripts (`12`-`14`). It was archived because project policy changed today: real VEDA4
now calculates the PED, and `ped/` at the top level holds only the tooling that
*prepares VEDA4's input* (see `ped/build_veda_fmt.py`, generalized from `12`). This
code no longer runs as the manuscript's live PED source.

What's still trustworthy here: this is the exact code that produced the manuscript's
*previous* Table 6 numbers (`01`-`05`, cross-checked internally by `03`-`04`'s PED
column sums), and the completed `13`/`14` cross-validation against real VEDA4 output
is a genuine, already-finished piece of validation work — 3 of the 4 manuscript-cited
modes agreed closely with real VEDA4 (see `veda4_comparison.csv`/`.txt`), with one
real disagreement at 1056.39 cm⁻¹ where CCC-bend/CCH-bend character comes out swapped
between the two methods, attributable to degenerate-eigenvector rotation ambiguity in
a doubly-degenerate E-symmetry mode pair, not a bug in either method.

For the current active tool, see `ped/README.md`.
