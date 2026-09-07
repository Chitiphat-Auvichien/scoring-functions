"""Source for the ``ref_label``/citation-key content that
``src/library_ingest.py`` attaches to every internal-vibration row of
``library_scores.csv``. Two author-maintained CSVs under ``data/`` supply it:

- ``data/characterised_modes.csv``: one row per (molecule, 1-based internal
  mode index). Engine-derivable columns (``freq``/``mu``/``k``) are
  regenerated from disk by ``src.library_ingest.regenerate_characterised_
  modes()``, while manually-curated literature columns are preserved --
  ``type`` (the literature ``ref_label``, stored abbreviated as
  ``"B"``/``"S"``/``"SB"`` -- translated to ``"bend"``/``"stretch"``/``"SB"``
  by ``_filtered_ref_label()`` below; ``"SB"`` is a genuine 3rd class for two
  benzene modes) and ``ref`` (a citation
  key, e.g. ``"Shi1972"``; only partially back-filled, not a bug).
- ``ref-label_citation.csv``: citation key -> bibliographic metadata (doi,
  1st author, journal, year, note).

``data/data_score.csv`` is retired: its ``type``/``ideal`` content and
frequency-agreement gate are now both served by ``characterised_modes.csv``
alone (``ideal`` itself comes structurally from ``mol_list_method.csv``'s
``mol_type`` column via ``attach_ideal_tags()``); no code path reads it.

Coverage gap (out of scope, not a bug to patch): these CSVs omit six
multi-centre molecules (``C2H2``, ``C2H4``, ``C2H6``, ``H2O2``,
``iso-C4H10``, ``n-C4H10``) that are also absent from the finalized
72-molecule roster.

Public API (consumed by ``src/library_ingest.py``):
  ``load_label_csvs(data_dir)`` -> dict of the two raw DataFrames
      (``"characterised_modes"``, ``"citations"``).
  ``build_label_lookup(tables)`` -> {(molecule, mode:int):
      {"ref_label", "ref_key"}}.
  ``get_label(lookup, molecule, mode)`` -> (ref_label, ref_key) tuple,
      (None, None) if the (molecule, mode) key is absent.
  ``build_citation_table(tables)`` -> {code: {"doi", "author", "journal",
      "year", "note"}}, a side lookup for the manuscript to cite by key
      (kept separate from ``library_scores.csv``'s locked numeric columns;
      ``ref_key`` is the only citation-related column added to that schema).
"""
import os

import pandas as pd

# Literal literature ref_label values recognized; anything else (typos,
# blanks) in the 'type' column maps to None rather than being accepted silently.
ALLOWED_REF_LABELS = ("bend", "stretch", "SB")

# characterised_modes.csv's 'type' column stores the abbreviated form
# ("B"/"S"/"SB"); ref_label (consumed by calibrate.py/figures.py) stays the
# full word internally. Public so any other direct reader of the raw CSV
# (e.g. src/figures.py's irrep-coupling figure) can normalize the same way.
TYPE_TO_REF_LABEL = {"B": "bend", "S": "stretch", "SB": "SB"}

LABEL_CSV_FILES = {
    "characterised_modes": "characterised_modes.csv",
    "citations": "ref-label_citation.csv",
}


def load_label_csvs(data_dir="data"):
    """Read the two label CSVs from `data_dir`. Raises FileNotFoundError (via
    pandas) if either is missing -- fail loud, this is the intended default
    source now, not an optional extra."""
    tables = {}
    for key, fname in LABEL_CSV_FILES.items():
        tables[key] = pd.read_csv(os.path.join(data_dir, fname))
    return tables


def _filtered_ref_label(raw_type):
    if pd.isna(raw_type):
        return None
    raw_type = str(raw_type).strip()
    ref_label = TYPE_TO_REF_LABEL.get(raw_type, raw_type)
    return ref_label if ref_label in ALLOWED_REF_LABELS else None


def build_label_lookup(tables):
    """Build {(molecule, mode:int): {"ref_label", "ref_key"}} from
    ``characterised_modes.csv`` alone."""
    cm = tables["characterised_modes"].copy()
    cm["mode"] = pd.to_numeric(cm["mode"], errors="coerce")

    lookup = {}
    for _, r in cm.iterrows():
        if pd.isna(r["mode"]):
            continue
        key = (r["molecule"], int(r["mode"]))
        ref_key = r.get("ref")
        lookup[key] = {
            "ref_label": _filtered_ref_label(r.get("type")),
            # DataFrame.iterrows() boxes a mixed-dtype row into one Series;
            # depending on pandas version this can turn a genuine None into
            # NaN when the row also holds a float column (e.g. 'mode' after
            # pd.to_numeric above) -- normalize back to None so downstream
            # consumers get a real missing-value sentinel either way.
            "ref_key": None if pd.isna(ref_key) else ref_key,
        }
    return lookup


def get_label(lookup, molecule, mode):
    """(ref_label, ref_key) for `(molecule, mode)`, or (None, None) if not
    found in the label CSVs."""
    entry = lookup.get((molecule, int(mode)))
    if entry is None:
        return None, None
    return entry["ref_label"], entry["ref_key"]


def build_citation_table(tables):
    """{code: {"doi", "author", "journal", "year", "note"}} from
    ``ref-label_citation.csv`` -- a side lookup, not merged into
    ``library_scores.csv`` beyond the ``ref_key`` column itself, so the
    manuscript can eventually resolve a key to full bibliographic detail."""
    rc = tables["citations"]
    out = {}
    for _, r in rc.iterrows():
        code = r.get("code")
        if pd.isna(code) or not str(code).strip():
            continue
        out[str(code).strip()] = {
            "doi": r.get("doi"),
            "author": r.get("1st author"),
            "journal": r.get("journal"),
            "year": r.get("year"),
            "note": r.get("note"),
        }
    return out
