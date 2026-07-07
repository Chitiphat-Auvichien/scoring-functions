"""Library ingest: builds data/results/library_scores.csv for every molecule
in the JCC paper's roster, ``data/mol_list_method.csv`` (72 rows: 11 "ideal"
single-centre AB_n shapes, 60 "non-ideal" substituted variants of those same
shapes, 1 "multi-centre" = benzene).

**Single, roster-driven pipeline (2026-07-07 final flip to Gaussian-direct).**
This module used to switch between two data sources -- precomputed scores
read out of a spreadsheet workbook (the interim "source=excel" path, adopted
2026-07-04 while only ~25/70 molecules had on-disk Gaussian geometry) and a
disk-driven real-engine recompute ("source=gaussian"). That dual-source
scaffolding is now retired outright, not just defaulted away from: every one
of the 72 roster molecules has a verified, on-disk ``.log``/``.gjf`` (or
``.com``) pair (Phase A of the 2026-07-07 plan: ``mol_list_method.csv``
gained a ``basename`` column mapping each canonical ``molecule`` name to its
actual on-disk file stem, e.g. ``SbH3`` -> ``SbH3_MP2_3-21G``, ``BrF3`` ->
``brf3``). There is exactly one way to build ``library_scores.csv`` now:
iterate the roster, run the real engine (Steps 1-4, ``main.load_inputs`` ->
``main.build_scorer_and_final`` -> ``src.classifier.classify_all_modes``,
``ModeScorer.score_bonds()``) on every row's basename, and label the output
with the roster's own canonical ``molecule`` name. This module never opens
the old spreadsheet workbook at all.

Pipeline
--------
1. ``load_mol_roster()`` reads ``data/mol_list_method.csv`` -- the
   authoritative (molecule, basename) registry.
2. ``check_roster_disk_consistency()`` cross-checks every roster basename
   against ``discover_geometry_molecules()``'s directory-intersection scan
   (repurposed from primary discovery, which it used to be under
   "source=gaussian", to a safety-net consistency check). Any roster row
   whose basename does NOT resolve to a real ``.log``+``.gjf``/``.com`` pair
   on disk is a real regression now that 100% coverage is the expected
   state -- ``build_library_scores()`` raises ``FileNotFoundError`` listing
   every such (molecule, basename) pair before doing any per-molecule work
   (fail fast). An on-disk basename with no roster row at all (e.g. the
   gramicidin ``1grm_MM_UFF`` companion-paper inputs -- Decision 5, out of
   scope for this paper) is not an error; it triggers a non-fatal
   ``warnings.warn`` and is simply excluded from the output.
3. ``score_geometry_molecule(base, ...)`` (unchanged since 2026-07-03) runs
   the real engine on one on-disk molecule, returning one row per external
   T/R slot and one row per internal "Vib i" mode, with ``ref_label``/
   ``ideal``/``ref_key`` left ``None`` (label-only content, attached
   separately -- see below). A molecule whose files are present but fail to
   parse/score (``FileNotFoundError``/``ValueError`` from
   ``build_scorer_and_final``, e.g. no usable bond connectivity) is a
   DIFFERENT failure mode than "missing from disk" -- that is warned about
   and skipped (excluded from the CSV), not raised, since Phase B's
   fail-loud guarantee is specifically about roster-vs-disk coverage, not
   every possible parse edge case.
4. ``attach_labels()`` joins ``ref_label``/``ideal``/``ref_key`` onto every
   internal row from ``src/csv_label_ingest.py``'s three author-maintained
   CSVs (``data/data_score.csv``, ``data/characterised_modes.csv``,
   ``data/ref-label_citation.csv``) -- see that module's docstring. The join
   is gated by a whole-molecule frequency-agreement check against
   ``data_score.csv``'s own tabulated ``freq`` column: 19 of the 72 roster
   molecules run at a fallback level of theory (``mol_list_method.csv``'s
   ``current_method`` column, e.g. B3LYP/3-21G instead of the default
   MP2/3-21G) rather than the exact calculation ``data_score.csv``'s
   independently-curated frequencies describe, so this gate still catches
   genuine mode-index mismatches between the engine's parsed frequency and
   ``data_score.csv``'s expectation for that (molecule, mode_index) -- one
   mismatched (or altogether absent) mode disqualifies the WHOLE molecule's
   label join, so a bad molecule cannot half-merge. External (T/R) rows are
   untouched by this gate (already correct, structural, geometry-only).
   Scores (``V_Stretch``, ``Tx..Rz``, ``predicted_label``, ``s_AB``, ...) are
   NEVER affected by a failed label join -- only ``ref_label``/``ideal``/
   ``ref_key`` are left null for that molecule's internal rows.

Output schema (unchanged column names so ``src/calibrate.py``/
``src/figures.py`` keep working):
  ``molecule``, ``mode_index``, ``kind``, ``freq``, ``ref_label``, ``ideal``,
  ``V_Stretch``, ``delta_b_mean``, ``s_AB``, ``rel_db``, ``has_geometry``,
  ``predicted_label``, ``predicted_annotation``, ``Tx``, ``Ty``, ``Tz``,
  ``Rx``, ``Ry``, ``Rz``, ``ref_key``. ``has_geometry`` is unconditionally
  ``True`` for every row now (every roster molecule has on-disk geometry by
  construction) -- kept in the schema rather than dropped so downstream
  consumers that still read it (e.g. ``src/calibrate.py``'s
  ``_load_geometry_pool``) do not need to change.
"""
import os
import warnings

import numpy as np
import pandas as pd

from src import csv_label_ingest

_EXTERNAL_SLOTS = ("Tx", "Ty", "Tz", "Rx", "Ry", "Rz")

# Locked output schema (src/calibrate.py and src/figures.py read these exact
# column names). `ref_key` (2026-07-05, citation key from
# src/csv_label_ingest.py, e.g. "Shi1972") is a purely-additive column
# appended at the end: existing consumers read columns by name, not
# position, so this does not disturb them.
SCHEMA_COLUMNS = [
    "molecule", "mode_index", "kind", "freq", "ref_label", "ideal",
    "V_Stretch", "delta_b_mean", "s_AB", "rel_db", "has_geometry",
    "predicted_label", "predicted_annotation",
    "Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "ref_key",
]


def load_mol_roster(data_dir="data"):
    """Read data/mol_list_method.csv: the authoritative molecule roster
    (canonical `molecule` name + on-disk `basename`). Fail loud if the file
    or either required column is missing."""
    path = os.path.join(data_dir, "mol_list_method.csv")
    df = pd.read_csv(path)
    missing_cols = {"molecule", "basename"} - set(df.columns)
    if missing_cols:
        raise ValueError(f"mol_list_method.csv missing required column(s): {missing_cols}")
    return df


def check_roster_disk_consistency(roster, data_dir="data"):
    """Cross-check roster['basename'] against discover_geometry_molecules()'s
    directory intersection. Returns (missing, orphaned):
      missing  -- [(molecule, basename), ...] roster rows whose .log+.gjf pair
                  is NOT present on disk.
      orphaned -- sorted list of on-disk basenames not referenced by any
                  roster row's basename column.
    """
    disk_bases = set(discover_geometry_molecules(data_dir))
    roster_bases = set(roster["basename"])
    missing = [(m, b) for m, b in zip(roster["molecule"], roster["basename"]) if b not in disk_bases]
    orphaned = sorted(disk_bases - roster_bases)
    return missing, orphaned


def resolve_log_basename(molecule, data_dir="data"):
    """molecule (canonical mol_list_method.csv name) -> on-disk basename, or
    None if not in the roster. Kept as its own function (name/signature
    unchanged) since src/calibrate.py calls it directly."""
    roster = load_mol_roster(data_dir)
    matches = roster.loc[roster["molecule"] == molecule, "basename"]
    return matches.iloc[0] if len(matches) else None


def discover_geometry_molecules(data_dir="data"):
    """Sorted list of basenames present in BOTH data/logs/ (.log or .out) and
    data/gjf/ (.com or .gjf) -- a directory-intersection scan. No longer the
    primary discovery mechanism (that's now the roster,
    data/mol_list_method.csv); kept as check_roster_disk_consistency()'s
    safety-net cross-check against what's actually on disk."""
    logs_dir = os.path.join(data_dir, "logs")
    gjf_dir = os.path.join(data_dir, "gjf")
    log_bases = {os.path.splitext(f)[0] for f in os.listdir(logs_dir)
                 if f.lower().endswith((".log", ".out"))}
    gjf_bases = {os.path.splitext(f)[0] for f in os.listdir(gjf_dir)
                 if f.lower().endswith((".com", ".gjf"))}
    return sorted(log_bases & gjf_bases)


def score_geometry_molecule(base, data_dir="data", thresholds=None):
    """Run the real engine (Steps 1-4) on one on-disk molecule and return a
    list of row dicts (one per external T/R slot + one per internal 'Vib i'
    mode) in the library_scores.csv schema, EXCLUDING 'molecule' (the caller
    attaches that) and with ref_label/ideal/ref_key left None (attached
    separately by attach_labels, label-only).

    Raises ValueError if no bond connectivity is available (propagated from
    build_scorer_and_final -- a molecule with a .gjf but no usable
    connectivity is a real data problem, not silently skipped here; the
    caller decides whether to skip-and-warn).
    """
    from main import load_inputs, build_scorer_and_final
    from src.classifier import classify_all_modes

    raw, _ = load_inputs(base, "normal", data_dir)
    scorer, final = build_scorer_and_final(raw, "normal")
    by_name = {m.get("label", f"Mode {i + 1}"): m for i, m in enumerate(final)}
    scored = classify_all_modes(scorer, final, thresholds)

    rows = []
    for m in scored:
        name = m["name"]
        if name in _EXTERNAL_SLOTS:
            rows.append({
                "mode_index": name, "kind": "external",
                "freq": m["frequency"],
                "ref_label": "translation" if name[0] == "T" else "rotation",
                "ideal": None,
                "V_Stretch": m["V"], "delta_b_mean": None,
                "s_AB": "", "rel_db": "",
                "has_geometry": True,
                "predicted_label": m["classification"],
                "predicted_annotation": m["annotation"],
                "Tx": m["T"]["x"], "Ty": m["T"]["y"], "Tz": m["T"]["z"],
                "Rx": m["R"]["x"], "Ry": m["R"]["y"], "Rz": m["R"]["z"],
                "ref_key": None,  # T/R rows are structural/exact -- no literature citation
            })
            continue

        # Internal ("Vib i") mode: recompute per-bond detail (s_AB + rel_db)
        # for THIS mode's own displacement directly via score_bonds() --
        # classify_all_modes() only retains 'bonds' for STRETCHING/
        # MIXED_STRETCH_BEND labels (per its own docstring/spec), but every
        # internal row here wants per-bond detail regardless of label.
        # Cheap: same scorer/atoms, no reparsing, just reloads dispVec.
        mode_index = int(name.split()[1])
        mode_vec = by_name[name]["vector"]
        scorer.calculate_scores(mode_vec)
        bonds = scorer.score_bonds()

        # i_label/j_label are ModeScorer.score_bonds()'s own atom-symbol +
        # 1-based-index labels (e.g. "C1", "H7"), matching the convention
        # src/benzene_validation.py's C-C/C-H bond-type parsing relies on.
        s_ab_str = ";".join(f"{b['i_label']}-{b['j_label']}:{b['s_AB']:.4f}" for b in bonds)
        rel_db_str = ";".join(f"{b['i_label']}-{b['j_label']}:{b['rel_db']:.4f}" for b in bonds)
        delta_b_mean = float(np.mean([abs(b["rel_db"]) for b in bonds])) if bonds else None

        rows.append({
            "mode_index": mode_index, "kind": "internal",
            "freq": m["frequency"],
            "ref_label": None, "ideal": None,  # attached later, label-only
            "V_Stretch": m["V"], "delta_b_mean": delta_b_mean,
            "s_AB": s_ab_str, "rel_db": rel_db_str,
            "has_geometry": True,
            "predicted_label": m["classification"],
            "predicted_annotation": m["annotation"],
            "Tx": m["T"]["x"], "Ty": m["T"]["y"], "Tz": m["T"]["z"],
            "Rx": m["R"]["x"], "Ry": m["R"]["y"], "Rz": m["R"]["z"],
            "ref_key": None,  # attached later, label-only
        })
    return rows


def attach_labels(df, csv_tables, label_lookup, freq_atol=0.05, freq_rtol=1e-4):
    """Join ref_label/ideal/ref_key onto `df`'s internal rows, molecule by
    molecule, gated by a whole-molecule frequency-agreement check against
    `csv_tables["data_score"]` (data/data_score.csv, see module docstring --
    this gate is NOT vestigial: it catches genuine mode-index mismatches
    between the engine's own parsed frequency and data_score.csv's
    independently-curated expectation for that (molecule, mode_index)).
    External rows are untouched (already correct, structural).

    Once a molecule passes the gate, the actual ref_label/ideal/ref_key
    VALUES written come from `label_lookup`
    (``src.csv_label_ingest.build_label_lookup()``), not from the
    data_score.csv row directly; that row is used only for the frequency
    gate.

    Returns (df, skip_report) where skip_report is a list of {'molecule',
    'n_mismatched', 'example': (mode_index, engine_freq, ds_freq)} dicts, one
    per molecule whose internal-row label join was skipped (ds_freq is None
    if no data_score.csv row exists at all for that mode index -- treated
    identically to a numeric mismatch, not silently skipped, per the
    fail-loud guarantee established in a prior session). A molecule entirely
    absent from data_score.csv (no row at all, e.g. a naming mismatch or a
    molecule the author has not yet tabulated) is left untouched with no
    skip-report entry -- not an error, just no ground truth to gate against.
    """
    ds = csv_tables["data_score"].copy()
    ds["freq"] = pd.to_numeric(ds["freq"], errors="coerce")
    ds["mode"] = pd.to_numeric(ds["mode"], errors="coerce")
    ds_molecules = set(ds["molecule"].dropna().unique())

    ds_by_key = {}
    for _, r in ds.iterrows():
        if pd.isna(r["mode"]):
            continue
        ds_by_key[(r["molecule"], int(r["mode"]))] = r

    df = df.copy()
    skip_report = []
    for mol in df["molecule"].unique():
        if mol not in ds_molecules:
            continue  # no data_score.csv counterpart at all -- not an error

        internal_idx = df.index[(df["molecule"] == mol) & (df["kind"] == "internal")]
        matched = []
        mismatches = []
        for idx in internal_idx:
            mode_index = int(df.at[idx, "mode_index"])
            engine_freq = df.at[idx, "freq"]
            ds_row = ds_by_key.get((mol, mode_index))
            if ds_row is None:
                mismatches.append((mode_index, engine_freq, None))
                continue
            ds_freq = ds_row["freq"]
            if pd.isna(ds_freq) or not np.isclose(engine_freq, ds_freq,
                                                    atol=freq_atol, rtol=freq_rtol):
                mismatches.append((mode_index, engine_freq,
                                    None if pd.isna(ds_freq) else float(ds_freq)))
                continue
            matched.append((idx, ds_row))

        if mismatches:
            skip_report.append({
                "molecule": mol, "n_mismatched": len(mismatches),
                "example": mismatches[0],
            })
            continue

        for idx, ds_row in matched:
            mode_index = int(df.at[idx, "mode_index"])
            ref_label, ideal, ref_key = csv_label_ingest.get_label(label_lookup, mol, mode_index)
            df.at[idx, "ref_label"] = ref_label
            df.at[idx, "ideal"] = ideal
            df.at[idx, "ref_key"] = ref_key

    return df, skip_report


def _build_library_scores(data_dir, thresholds, return_skip_report):
    """The single builder: iterate every roster row, score it via the real
    engine, and label-join ref_label/ideal/ref_key from the CSV label
    sources. See module docstring."""
    roster = load_mol_roster(data_dir)
    missing, orphaned = check_roster_disk_consistency(roster, data_dir)

    if missing:
        raise FileNotFoundError(
            "library_ingest: the following mol_list_method.csv roster row(s) have "
            "no matching .log+.gjf/.com pair on disk (100% coverage is the "
            f"expected state -- this is a real regression): {missing}"
        )
    for base in orphaned:
        warnings.warn(
            f"library_ingest: on-disk basename '{base}' (data/logs + data/gjf) is "
            "not referenced by any row of data/mol_list_method.csv -- excluded "
            "from library_scores.csv. Expected for out-of-scope files (e.g. the "
            "gramicidin '1grm_MM_UFF' companion-paper inputs).",
            stacklevel=3)

    csv_tables = csv_label_ingest.load_label_csvs(data_dir)
    label_lookup = csv_label_ingest.build_label_lookup(csv_tables)

    all_rows = []
    load_errors = []
    for _, r in roster.iterrows():
        molecule, base = r["molecule"], r["basename"]
        try:
            rows = score_geometry_molecule(base, data_dir, thresholds)
        except (FileNotFoundError, ValueError) as e:
            load_errors.append((molecule, base, str(e)))
            continue
        for row in rows:
            row["molecule"] = molecule
        all_rows.extend(rows)

    df = pd.DataFrame(all_rows, columns=SCHEMA_COLUMNS)
    df, skip_report = attach_labels(df, csv_tables, label_lookup)

    for molecule, base, err in load_errors:
        warnings.warn(
            f"library_ingest: skipped roster molecule '{molecule}' (basename "
            f"'{base}') -- parse/scoring failed ({err}) -- not included in "
            "library_scores.csv at all.",
            stacklevel=3)
    for entry in skip_report:
        mode_index, eng_f, ds_f = entry["example"]
        ds_f_str = f"{ds_f:.4f}" if ds_f is not None else "NO data_score.csv ROW FOUND"
        warnings.warn(
            f"library_ingest: '{entry['molecule']}' internal rows NOT label-"
            f"joined ({entry['n_mismatched']} mode(s) mismatched, e.g. mode "
            f"{mode_index}: engine {eng_f:.4f} vs data_score.csv {ds_f_str} cm-1) -- "
            "this molecule's data_score.csv row likely came from a different "
            "calculation than data/logs/ (e.g. a fallback level of theory -- see "
            "mol_list_method.csv's current_method column). Its scores "
            "(V_Stretch, Tx..Rz, predicted_label, ...) are still the real "
            "engine's own and are NOT affected; only ref_label/ideal/ref_key "
            "are left null.",
            stacklevel=3)

    if return_skip_report:
        return df, skip_report
    return df


def build_library_scores(data_dir="data", thresholds=None, return_skip_report=False):
    """Build the library_scores DataFrame: real-engine recompute for every
    one of the 72 mol_list_method.csv roster molecules. Returns a DataFrame
    with exactly SCHEMA_COLUMNS.

    Raises FileNotFoundError if any roster row's .log/.gjf pair is missing
    from disk (see _build_library_scores / check_roster_disk_consistency).
    """
    return _build_library_scores(data_dir, thresholds, return_skip_report)


def run_ingest_pipeline(data_dir="data", thresholds=None, write=True):
    """Headless entry point: build library_scores and optionally write the
    CSV."""
    df, skip_report = build_library_scores(data_dir, thresholds, return_skip_report=True)
    out_path = os.path.join(data_dir, "results", "library_scores.csv")
    if write:
        df.to_csv(out_path, index=False)
    return df, out_path, skip_report
