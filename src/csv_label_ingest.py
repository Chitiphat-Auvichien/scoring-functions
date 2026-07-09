"""Default source (2026-07-05 author decision, sole source since the
2026-07-07 retirement of the old precomputed-spreadsheet ingest path) for the
``ref_label``/citation-key content that ``src/library_ingest.py`` attaches to
every internal-vibration row of ``library_scores.csv``.

Two author-maintained CSVs under ``data/`` supply this content:

- ``data/characterised_modes.csv``: one row per (molecule, 1-based internal
  mode index), regenerated from a direct on-disk scan of
  ``data/logs``/``data/gjf`` by ``src.library_ingest.regenerate_characterised_
  modes()`` for its engine-derivable columns (``freq``/``μ``/``k``), while
  PRESERVING the author's manually-curated literature columns -- most
  importantly ``type`` (the literature ``ref_label``: ``"bend"``/
  ``"stretch"``/``"SB"``, the last a genuine 3rd class for two of benzene's
  modes, not just the classifier's own predicted mixed bucket) and ``ref``
  (a citation key, e.g. ``"Shi1972"``, joined onto ``data_score.csv``... no
  longer -- see below -- joined onto this SAME file's own rows). Only a
  fraction of rows currently have a citation key filled in (the rest are
  blank -- a real, partial-coverage data limitation, not a bug: the author
  has not back-filled citations for the whole library yet).
- ``ref-label_citation.csv`` (a handful of rows): citation key -> full
  bibliographic metadata (doi, 1st author, journal, year, note).

**2026-07-08: ``data/data_score.csv`` retired from this module entirely.**
Prior to this session, ``data_score.csv`` supplied the SAME ``type``/``ideal``
content as ``characterised_modes.csv``'s ``type`` column (a near-duplicate,
maintained separately) plus the frequency-agreement gate
``src.library_ingest.attach_labels()`` used. Both roles are now served by
``characterised_modes.csv`` alone (its own ``freq`` column, freshly
regenerated straight from the on-disk log, is what the gate compares
against -- see ``src.library_ingest.attach_labels()``'s docstring), and
``ideal`` itself is no longer a per-mode literature label at all: it is
sourced structurally from ``data/mol_list_method.csv``'s per-molecule
``mol_type`` column by ``src.library_ingest.attach_ideal_tags()``, entirely
outside this module. ``data/data_score.csv`` was left on disk, unread by
any code path, from 2026-07-08 until 2026-07-09, when the author deleted it
entirely (an unrelated cleanup of an already-retired file, alongside the
OH4/OF4 library-exclusion session) -- ``src.library_ingest.resync_
reference_metadata()`` was updated the same session to stop reading/writing
it, so no code path references the file at all anymore.

**Coverage gap (2026-07-07: out of scope, not a gap to patch).** These CSVs
do not cover every molecule the old workbook-driven pipeline used to score
(``C2H2``, ``C2H4``, ``C2H6``, ``H2O2``, ``iso-C4H10``, ``n-C4H10``) -- but
all six are multi-centre molecules and are also absent from
``data/mol_list_method.csv``'s finalized 72-molecule roster (confirmed),
i.e. genuinely out of scope for this paper's Gaussian-direct pipeline, not a
coverage gap that needs a fallback. The old workbook-fallback that used to
cover them (``build_label_lookup``'s old ``fallback_ds`` parameter) has been
removed entirely along with the workbook-sourced ingest path it supported
(``src/library_ingest.py``, formerly ``src/excel_ingest.py``). Per-bond
``s_AB`` data is untouched by this module entirely; it comes from the real
engine's own recompute (``src/library_ingest.py::score_geometry_molecule``),
unrelated to label/citation content.

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

# Literal literature ref_label values this framework recognizes. Anything
# else in the 'type' column (typos, blanks) maps to None rather than being
# silently accepted -- same fail-quiet-but-not-fail-wrong policy the old
# workbook-only code used for ("bend", "stretch"), extended to the genuine
# 3rd class "SB" (mixed stretch/bend), per the 2026-07-05 benzene relabeling.
ALLOWED_REF_LABELS = ("bend", "stretch", "SB")

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
    return raw_type if raw_type in ALLOWED_REF_LABELS else None


def build_label_lookup(tables):
    """Build {(molecule, mode:int): {"ref_label", "ref_key"}} from
    ``characterised_modes.csv`` alone (2026-07-08: ``data_score.csv`` is no
    longer read here -- see module docstring). See module docstring for the
    molecules these CSVs do not cover (no fallback needed -- they are
    absent from the 72-molecule roster too).
    """
    cm = tables["characterised_modes"].copy()
    cm["mode"] = pd.to_numeric(cm["mode"], errors="coerce")

    lookup = {}
    for _, r in cm.iterrows():
        if pd.isna(r["mode"]):
            continue
        key = (r["molecule"], int(r["mode"]))
        lookup[key] = {
            "ref_label": _filtered_ref_label(r.get("type")),
            "ref_key": r.get("ref"),
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
