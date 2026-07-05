"""New default source (2026-07-05 author decision) for the ``ref_label``/
``ideal``/citation-key content that ``src/excel_ingest.py`` attaches to every
internal-vibration row of ``library_scores.csv``, across BOTH
``source="excel"`` and ``source="gaussian"``.

Three author-maintained CSVs under ``data/`` supersede the equivalent columns
of the xlsx workbook's ``data_score``/``characterised modes`` sheets:

- ``data_score.csv`` (63 molecules): one row per (molecule, 1-based internal
  mode index) -- ``type`` (the literature ``ref_label``: now literally
  ``"bend"``/``"stretch"``/``"SB"``, a genuine 3rd class for two of benzene's
  modes, not just the classifier's own predicted mixed bucket) and ``ideal``
  (``"yes"``/``"no"``). This is a straight re-export of the xlsx
  ``data_score`` sheet's own columns plus the author's benzene relabeling --
  hand-verified to reproduce identical values to the xlsx for every molecule
  except ``C6H6`` (mode_index 19/23/24 bend->stretch, 21/22 bend->SB).
- ``characterised_modes.csv`` (63 molecules): richer per-mode literature
  detail, most importantly the ``ref`` column -- a citation key (e.g.
  ``"Shi1972"``) joined onto ``data_score.csv`` rows by (molecule, mode).
  Only ~11 of the 63 molecules currently have a citation key filled in
  (the rest are blank -- a real, partial-coverage data limitation, not a
  bug: the author has not back-filled citations for the whole library yet).
- ``ref-label_citation.csv`` (a handful of rows): citation key -> full
  bibliographic metadata (doi, 1st author, journal, year, note).

**Coverage gap and its resolution.** These three CSVs cover 63 molecules --
6 molecules the OLD xlsx-driven pipeline scored (``C2H2``, ``C2H4``,
``C2H6``, ``H2O2``, ``iso-C4H10``, ``n-C4H10``) are entirely absent from the
new CSVs (not yet migrated by the author). To honor the "byte-for-byte
identical except C6H6" regression requirement, ``build_label_lookup()``
takes an optional ``fallback_ds`` (the xlsx ``data_score`` sheet, read by
``excel_ingest.load_excel_tables``) and uses it ONLY for molecules the new
CSVs do not cover at all -- every molecule the new CSVs DO cover always wins
over the xlsx fallback, per the "full switch" instruction. Per-bond ``s_AB``
data is untouched by this module entirely; it keeps coming from its existing
source (Gaussian recompute or the xlsx ``data_mode&bond`` sheet), unrelated
to label/citation content.

Public API (consumed by ``src/excel_ingest.py``):
  ``load_label_csvs(data_dir)`` -> dict of the three raw DataFrames.
  ``build_label_lookup(tables, fallback_ds=None)`` -> {(molecule, mode:int):
      {"ref_label", "ideal", "ref_key"}}.
  ``get_label(lookup, molecule, mode)`` -> (ref_label, ideal, ref_key) tuple,
      (None, None, None) if the (molecule, mode) key is absent everywhere.
  ``build_citation_table(tables)`` -> {code: {"doi", "author", "journal",
      "year", "note"}}, a side lookup for the manuscript to cite by key
      (kept separate from ``library_scores.csv``'s locked numeric columns;
      ``ref_key`` is the only citation-related column added to that schema).
"""
import os

import pandas as pd

# Literal literature ref_label values this framework recognizes. Anything
# else in the 'type' column (typos, blanks) maps to None rather than being
# silently accepted -- same fail-quiet-but-not-fail-wrong policy the old
# xlsx-only code used for ("bend", "stretch"), now extended to the new
# genuine 3rd class "SB" (mixed stretch/bend), per the 2026-07-05 benzene
# relabeling.
ALLOWED_REF_LABELS = ("bend", "stretch", "SB")

LABEL_CSV_FILES = {
    "data_score": "data_score.csv",
    "characterised_modes": "characterised_modes.csv",
    "citations": "ref-label_citation.csv",
}


def load_label_csvs(data_dir="data"):
    """Read the three new CSVs from `data_dir`. Raises FileNotFoundError (via
    pandas) if any is missing -- fail loud, this is the intended default
    source now, not an optional extra."""
    tables = {}
    for key, fname in LABEL_CSV_FILES.items():
        tables[key] = pd.read_csv(os.path.join(data_dir, fname))
    return tables


def _filtered_ref_label(raw_type):
    if pd.isna(raw_type):
        return None
    raw_type = str(raw_type).strip()
    return raw_type if raw_type in ALLOWED_REF_LABELS else None


def _xlsx_fallback_ref_label(raw_type):
    """Old xlsx-only filter (bend/stretch only -- the xlsx sheet never had a
    literal 'SB' value, so there is nothing to extend here for the fallback
    path)."""
    if pd.isna(raw_type):
        return None
    raw_type = str(raw_type).strip()
    return raw_type if raw_type in ("bend", "stretch") else None


def build_label_lookup(tables, fallback_ds=None):
    """Build {(molecule, mode:int): {"ref_label", "ideal", "ref_key"}} with
    the new CSVs as the authoritative source and `fallback_ds` (the xlsx
    'data_score' sheet DataFrame, as returned by
    ``excel_ingest.load_excel_tables()["data_score"]``) supplying entries
    ONLY for molecules the new CSVs do not cover at all. See module
    docstring for why this fallback exists and which 6 molecules use it.
    """
    ds = tables["data_score"].copy()
    ds["mode"] = pd.to_numeric(ds["mode"], errors="coerce")
    cm = tables["characterised_modes"].copy()
    cm["mode"] = pd.to_numeric(cm["mode"], errors="coerce")

    ref_key_by_key = {}
    for _, r in cm.iterrows():
        if pd.isna(r["mode"]):
            continue
        ref_key_by_key[(r["molecule"], int(r["mode"]))] = r.get("ref")

    lookup = {}

    if fallback_ds is not None:
        fds = fallback_ds.copy()
        fds["mode"] = pd.to_numeric(fds["mode"], errors="coerce")
        new_csv_molecules = set(ds["molecule"].dropna().unique())
        for _, r in fds.iterrows():
            if pd.isna(r["mode"]) or r["molecule"] in new_csv_molecules:
                continue  # new CSVs cover this molecule -- never fall back for it
            lookup[(r["molecule"], int(r["mode"]))] = {
                "ref_label": _xlsx_fallback_ref_label(r["type"]),
                "ideal": r["ideal"],
                "ref_key": None,
            }

    for _, r in ds.iterrows():
        if pd.isna(r["mode"]):
            continue
        key = (r["molecule"], int(r["mode"]))
        lookup[key] = {
            "ref_label": _filtered_ref_label(r["type"]),
            "ideal": r["ideal"],
            "ref_key": ref_key_by_key.get(key),
        }
    return lookup


def get_label(lookup, molecule, mode):
    """(ref_label, ideal, ref_key) for `(molecule, mode)`, or (None, None,
    None) if not found in either the new CSVs or the xlsx fallback."""
    entry = lookup.get((molecule, int(mode)))
    if entry is None:
        return None, None, None
    return entry["ref_label"], entry["ideal"], entry["ref_key"]


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
