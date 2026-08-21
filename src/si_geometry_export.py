"""Cartesian-coordinate longtable fragments for the JCC Supporting
Information, closing the IMPLEMENTATION_PLAN.md "SI Cartesian-geometry
export" TODO for water and benzene.

Reuses GaussianParser's already-parsed Standard-orientation geometry
(src/parser.py:115-144, self.coordinates/self.atom_symbols) plus
get_symbol() -- no new parsing logic, this module only formats.

    python -m src.si_geometry_export       # water + benzene

Scope for this pass is water and benzene only (the two structures headlined
in the main results). CO2/gramicidin/library geometries are a deferred
follow-up -- see the split checklist item in IMPLEMENTATION_PLAN.md.
"""

import os

from .parser import GaussianParser
from .si_tables import (_esc, _fmt, _write_longtable,
                         DEFAULT_DATA_DIR, DEFAULT_SI_TABLES_DIR)

# key -> (log basename under data/logs/, display name for the caption,
#         \label{}, output filename under JCC_SI/tables/)
GEOMETRY_MOLECULES = {
    "water": ("H2O", "H$_2$O", "tab:cartwater", "tab_cartesian_water.tex"),
    "benzene": ("C6H6", "C$_6$H$_6$ (benzene)", "tab:cartbenzene", "tab_cartesian_benzene.tex"),
}


def cartesian_table(molecule_name, log_path, label, out_path):
    """Parse `log_path`'s Standard-orientation geometry and write a longtable
    fragment of (Atom, X, Y, Z) in Angstrom to `out_path`."""
    parser = GaussianParser(log_path)
    parser.parse(parse_modes=False)

    headers = ["Atom", "X (\\AA)", "Y (\\AA)", "Z (\\AA)"]
    col_spec = "lrrr"
    rows = []
    for i, (sym, xyz) in enumerate(zip(parser.atom_symbols, parser.coordinates), start=1):
        x, y, z = xyz
        rows.append([_esc(f"{sym}{i}"), _fmt(x, 6), _fmt(y, 6), _fmt(z, 6)])

    caption = (f"Cartesian coordinates (\\AA, Gaussian ``Standard orientation'') of the "
               f"optimized geometry of {molecule_name} used for all normal-mode and "
               f"EMIT-mode scoring in the main text. Atom labels match the bond "
               f"labels used in the $s^{{AB}}$ score tables.")
    return _write_longtable(rows, headers, col_spec, caption, label, out_path)


def run_all(data_dir=DEFAULT_DATA_DIR, out_dir=DEFAULT_SI_TABLES_DIR):
    paths = []
    for mol, display_name, label, out_name in GEOMETRY_MOLECULES.values():
        log_path = os.path.join(data_dir, "logs", f"{mol}.log")
        paths.append(cartesian_table(display_name, log_path, label,
                                      os.path.join(out_dir, out_name)))
    return paths


if __name__ == "__main__":
    for p in run_all():
        print(p)
