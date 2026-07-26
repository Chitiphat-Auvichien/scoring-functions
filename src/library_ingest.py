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
actual on-disk file stem, e.g. ``SbH3`` -> ``SbH3``, ``BrF3`` ->
``BrF3``). There is exactly one way to build ``library_scores.csv`` now:
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
3. ``score_geometry_molecule(base, ..., mol_type=...)`` runs the real engine
   on one on-disk molecule, returning one row per external T/R slot and one
   row per internal "Vib i" mode, with ``ref_label``/``ideal``/``ref_key``
   left ``None`` (label-only content, attached separately -- see below).
   ``mol_type`` (the roster's own ``mol_type`` column for this molecule)
   gates the ``d_CA`` column: only computed for ``ideal``/``non-ideal``
   molecules with a unique central/hub atom -- see ``_central_atom_index()``.
   A molecule whose files are present but fail to parse/score
   (``FileNotFoundError``/``ValueError`` from ``build_scorer_and_final``,
   e.g. no usable bond connectivity) is a DIFFERENT failure mode than
   "missing from disk" -- that is warned about and skipped (excluded from
   the CSV), not raised, since Phase B's fail-loud guarantee is specifically
   about roster-vs-disk coverage, not every possible parse edge case.
4. ``attach_labels()`` joins ``ref_label``/``ref_key`` onto every internal
   row from ``src/csv_label_ingest.py``'s CSVs (``data/characterised_modes.csv``,
   ``data/ref-label_citation.csv`` -- see that module's docstring;
   ``data/data_score.csv`` is retired from this code path entirely as of
   2026-07-08, see below). The join is gated by a whole-molecule
   frequency-agreement check against ``characterised_modes.csv``'s own
   tabulated ``freq`` column (itself regenerated straight from the on-disk
   log by ``regenerate_characterised_modes()``, so a gate failure should now
   only mean a real mode-count/index problem, not staleness -- kept as a
   safety net regardless): 19 of the 72 roster molecules run at a fallback
   level of theory (``mol_list_method.csv``'s ``current_method`` column,
   e.g. B3LYP/3-21G instead of the default MP2/3-21G), so this gate still
   catches genuine mode-index mismatches between the engine's parsed
   frequency and ``characterised_modes.csv``'s expectation for that
   (molecule, mode_index) -- one mismatched (or altogether absent) mode
   disqualifies the WHOLE molecule's label join, so a bad molecule cannot
   half-merge. External (T/R) rows are untouched by this gate (already
   correct, structural, geometry-only). Scores (``V_Stretch``, ``Tx..Rz``,
   ``predicted_label``, ``s_AB``, ...) are NEVER affected by a failed label
   join -- only ``ref_label``/``ref_key`` are left null for that molecule's
   internal rows. ``attach_ideal_tags()`` separately (unconditionally, not
   gated by this frequency check) populates ``ideal`` for every internal
   row from ``mol_list_method.csv``'s per-molecule ``mol_type`` column --
   see that function's docstring.

**``data/data_score.csv`` retirement (2026-07-08, completed 2026-07-09).**
This module (and ``src/csv_label_ingest.py``) no longer read
``data/data_score.csv`` at all -- every quantity it used to supply is now
either computed by the engine (``d_CA``, replacing its old ``|d_CA|``
column) or sourced from ``data/characterised_modes.csv``/
``data/mol_list_method.csv`` (``ref_label``/``ref_key``/frequency-gate,
``ideal``, respectively). The 2026-07-08 retirement left the file on disk
untouched as a legacy artifact still synced by ``resync_reference_metadata()``;
2026-07-09 (OH4/OF4 exclusion session) the author deleted
``data/data_score.csv`` from disk entirely as a separate, unrelated cleanup
of that already-retired file, so ``resync_reference_metadata()`` was updated
in the same session to stop touching it -- it now resyncs
``characterised_modes.csv`` only. ``regenerate_characterised_modes()`` keeps
``characterised_modes.csv`` itself in sync with the on-disk logs (a direct
disk scan, not restricted to the roster) while preserving the author's
manually-curated literature columns.

Output schema (unchanged column names so ``src/calibrate.py``/
``src/figures.py`` keep working):
  ``molecule``, ``mode_index``, ``kind``, ``freq``, ``ref_label``, ``ideal``,
  ``V_Stretch``, ``delta_b_mean``, ``s_AB``, ``rel_db``, ``has_geometry``,
  ``predicted_label``, ``predicted_annotation``, ``Tx``, ``Ty``, ``Tz``,
  ``Rx``, ``Ry``, ``Rz``, ``ref_key``, ``reduced_mass``, ``force_constant``,
  ``irrep``, ``d_CA`` (additive, 2026-07-08 -- see SCHEMA_COLUMNS' own
  comment). ``has_geometry`` is unconditionally ``True`` for every row now
  (every roster molecule has on-disk geometry by construction) -- kept in
  the schema rather than dropped so downstream consumers that still read it
  (e.g. ``src/calibrate.py``'s ``_load_geometry_pool``) do not need to
  change.
"""
import os
import warnings

import numpy as np
import pandas as pd

from src import csv_label_ingest
from src.parser import GaussianParser

_EXTERNAL_SLOTS = ("Tx", "Ty", "Tz", "Rx", "Ry", "Rz")

# Locked output schema (src/calibrate.py and src/figures.py read these exact
# column names). `ref_key` (2026-07-05, citation key from
# src/csv_label_ingest.py, e.g. "Shi1972") is a purely-additive column
# appended at the end: existing consumers read columns by name, not
# position, so this does not disturb them. `reduced_mass`/`force_constant`/
# `irrep` (2026-07-08, Gaussian-direct parser rework) are appended the same
# way -- engine-parsed metadata for internal rows, None/blank for external
# (T/R) rows since construct_T/construct_R never set these keys. Named
# `reduced_mass`/`force_constant` (NOT bare `k`) to avoid colliding with this
# same CSV's existing, differently-scoped `k` column (a scoring metric
# adjacent to `vib_scr`, unrelated to force constant).
SCHEMA_COLUMNS = [
    "molecule", "mode_index", "kind", "freq", "ref_label", "ideal",
    "V_Stretch", "delta_b_mean", "s_AB", "rel_db", "has_geometry",
    "predicted_label", "predicted_annotation",
    "Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "ref_key",
    "reduced_mass", "force_constant", "irrep",
    # d_CA (2026-07-08, data_score.csv retirement): central/hub-atom
    # displacement amplitude for ONE internal mode -- ||mode_vector[central
    # atom]||. Populated only for internal rows of molecules tagged
    # mol_type=='ideal'/'non-ideal' in mol_list_method.csv (a genuine
    # single-hub AB_n topology) AND whose bond graph resolves to exactly one
    # atom bonded to every other atom (see _central_atom_index()); null for
    # every external row, every multi-centre molecule (e.g. C6H6 -- no
    # single hub exists), and any ideal/non-ideal molecule where the degree
    # check itself fails to find a unique hub. This replaces
    # data_score.csv's old, no-longer-read `|d_CA|` column with a real,
    # engine-derived quantity computed directly from bonds + the mode's own
    # displacement vector -- see score_geometry_molecule()/
    # _central_atom_index() below.
    "d_CA",
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


def _central_atom_index(bonds, n_atoms):
    """0-based index of the unique atom bonded to every other atom (degree
    == n_atoms - 1 -- the AB_n hub of a genuine single-centre topology), or
    None if no such atom exists or more than one does (multi-centre / no
    well-defined hub -- e.g. benzene's six degree-2 ring carbons, none of
    which is bonded to the other five). `bonds` is a list of 0-based
    (i, j) tuples (ModeScorer.bList / GaussianParser's parsed connectivity).
    """
    if n_atoms <= 1:
        return None
    degree = [0] * n_atoms
    for i, j in bonds:
        degree[i] += 1
        degree[j] += 1
    candidates = [a for a in range(n_atoms) if degree[a] == n_atoms - 1]
    return candidates[0] if len(candidates) == 1 else None


def score_geometry_molecule(base, data_dir="data", thresholds=None, mol_type=None):
    """Run the real engine (Steps 1-4) on one on-disk molecule and return a
    list of row dicts (one per external T/R slot + one per internal 'Vib i'
    mode) in the library_scores.csv schema, EXCLUDING 'molecule' (the caller
    attaches that) and with ref_label/ideal/ref_key left None (attached
    separately by attach_labels/attach_ideal_tags, label-only).

    `mol_type` (mol_list_method.csv's 'ideal'/'non-ideal'/'multi-centre'
    column for this molecule, optional) gates the new `d_CA` column: only
    computed for 'ideal'/'non-ideal' molecules (a genuine single-hub AB_n
    topology is the whole premise of "central atom"), and even then only if
    _central_atom_index() finds exactly one atom bonded to every other atom.
    `mol_type=None` (the default, e.g. direct test calls that don't go
    through the roster) never computes d_CA -- same as an explicit
    'multi-centre' tag.

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

    central_idx = None
    if mol_type in ("ideal", "non-ideal"):
        central_idx = _central_atom_index(scorer.bList, scorer.n)

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
                # construct_T/construct_R never set these (synthetic, not
                # parsed from a Gaussian frequency block) -- None by design.
                "reduced_mass": m.get("reduced_mass"),
                "force_constant": m.get("force_constant"),
                "irrep": m.get("irrep"),
                "d_CA": None,  # central-atom amplitude is an INTERNAL-mode
                                # quantity only; construct_T/construct_R
                                # slots never set it.
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

        d_ca = float(np.linalg.norm(mode_vec[central_idx])) if central_idx is not None else None

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
            "reduced_mass": m.get("reduced_mass"),
            "force_constant": m.get("force_constant"),
            "irrep": m.get("irrep"),
            "d_CA": d_ca,
        })
    return rows


def multi_centre_molecules(data_dir="data"):
    """Molecule names tagged mol_type=='multi-centre' in
    mol_list_method.csv -- the authoritative multi-centre / no-single-hub-
    atom classification. Repoints (2026-07-08) src/calibrate.py's formerly
    hardcoded SINGLE_CENTRE_ONLY_EXCLUDE frozenset (C2H2, C2H4, C2H6, H2O2,
    C6H6, iso-C4H10, n-C4H10 -- the 2026-07-05 scope decision) to this
    single, roster-driven source. Six of those seven names are outside the
    finalized 72-molecule roster entirely (removed 2026-07-07) and were
    already no-ops for filter_single_centre_library() on the real,
    roster-driven library_scores.csv -- only C6H6 was ever actually present
    to be dropped by that filter. mol_list_method.csv's mol_type column
    already tags C6H6 'multi-centre' (verified, not assumed -- see
    IMPLEMENTATION_PLAN.md), so this function reproduces byte-identical
    FILTERING behavior to the old hardcoded set on every real/current
    dataset, while being forward-compatible with any future multi-centre
    molecule added to the roster (no code edit needed, just a CSV edit).
    """
    roster = load_mol_roster(data_dir)
    if "mol_type" not in roster.columns:
        raise ValueError("mol_list_method.csv is missing the 'mol_type' column")
    return frozenset(roster.loc[roster["mol_type"] == "multi-centre", "molecule"])


def attach_labels(df, csv_tables, label_lookup, freq_atol=0.05, freq_rtol=1e-4):
    """Join ref_label/ref_key onto `df`'s internal rows, molecule by
    molecule, gated by a whole-molecule frequency-agreement check against
    `csv_tables["characterised_modes"]` (data/characterised_modes.csv --
    2026-07-08: repointed off data/data_score.csv, which is no longer read
    anywhere in this module; see the module docstring). This gate is NOT
    vestigial: it catches genuine mode-index mismatches between the
    engine's own parsed frequency and characterised_modes.csv's
    independently-curated expectation for that (molecule, mode_index). Since
    characterised_modes.csv's own `freq` column is now itself regenerated
    directly from the same on-disk log the engine parses (see
    `regenerate_characterised_modes()`), a gate failure at this point should
    only ever mean a real mode-count/index problem, not staleness -- kept as
    a safety net against future manual edits drifting characterised_modes.csv
    out of sync, not removed. External rows are untouched (already correct,
    structural).

    Once a molecule passes the gate, the actual ref_label/ref_key VALUES
    written come from `label_lookup`
    (``src.csv_label_ingest.build_label_lookup()``), not from the
    characterised_modes.csv row directly; that row is used only for the
    frequency gate. `ideal` is NOT set here at all -- see
    `attach_ideal_tags()` below, which sources it from
    mol_list_method.csv's per-molecule `mol_type` column instead (a
    structural roster property, not a per-mode literature match, so it does
    not need this gate).

    Returns (df, skip_report) where skip_report is a list of {'molecule',
    'n_mismatched', 'example': (mode_index, engine_freq, cm_freq)} dicts, one
    per molecule whose internal-row label join was skipped (cm_freq is None
    if no characterised_modes.csv row exists at all for that mode index --
    treated identically to a numeric mismatch, not silently skipped, per the
    fail-loud guarantee established in a prior session). A molecule entirely
    absent from characterised_modes.csv (no row at all, e.g. one of the
    still-open T-shaped/see-saw families awaiting literature back-fill) is
    left untouched with no skip-report entry -- not an error, just no ground
    truth to gate against.
    """
    cm = csv_tables["characterised_modes"].copy()
    cm["freq"] = pd.to_numeric(cm["freq"], errors="coerce")
    cm["mode"] = pd.to_numeric(cm["mode"], errors="coerce")
    cm_molecules = set(cm["molecule"].dropna().unique())

    cm_by_key = {}
    for _, r in cm.iterrows():
        if pd.isna(r["mode"]):
            continue
        cm_by_key[(r["molecule"], int(r["mode"]))] = r

    df = df.copy()
    skip_report = []
    for mol in df["molecule"].unique():
        if mol not in cm_molecules:
            continue  # no characterised_modes.csv counterpart at all -- not an error

        internal_idx = df.index[(df["molecule"] == mol) & (df["kind"] == "internal")]
        matched = []
        mismatches = []
        for idx in internal_idx:
            mode_index = int(df.at[idx, "mode_index"])
            engine_freq = df.at[idx, "freq"]
            cm_row = cm_by_key.get((mol, mode_index))
            if cm_row is None:
                mismatches.append((mode_index, engine_freq, None))
                continue
            cm_freq = cm_row["freq"]
            if pd.isna(cm_freq) or not np.isclose(engine_freq, cm_freq,
                                                    atol=freq_atol, rtol=freq_rtol):
                mismatches.append((mode_index, engine_freq,
                                    None if pd.isna(cm_freq) else float(cm_freq)))
                continue
            matched.append((idx, cm_row))

        if mismatches:
            skip_report.append({
                "molecule": mol, "n_mismatched": len(mismatches),
                "example": mismatches[0],
            })
            continue

        for idx, cm_row in matched:
            mode_index = int(df.at[idx, "mode_index"])
            ref_label, ref_key = csv_label_ingest.get_label(label_lookup, mol, mode_index)
            df.at[idx, "ref_label"] = ref_label
            df.at[idx, "ref_key"] = ref_key

    return df, skip_report


_MOL_TYPE_TO_IDEAL = {"ideal": "yes", "non-ideal": "no"}  # multi-centre -> None (n/a)


def attach_ideal_tags(df, roster):
    """Populate df['ideal'] for every INTERNAL row from
    mol_list_method.csv's per-molecule `mol_type` column ('ideal' -> 'yes',
    'non-ideal' -> 'no', 'multi-centre' -> None). Unconditional -- NOT gated
    by attach_labels()'s frequency-agreement check, because mol_type is a
    structural, roster-level property of the whole molecule, not a per-mode
    literature match that could be stale/mismatched the way ref_label is.
    External rows are left untouched (always None, as
    score_geometry_molecule() sets them) -- confusion_matrix_stats()'s
    ideal_tier_mask already treats every external row as ideal-tier via
    `kind == 'external'` regardless of the `ideal` column's value, so no
    consumer ever reads it there. Multi-centre molecules (mol_type ==
    'multi-centre', e.g. C6H6) get `ideal = None` -- they are already
    excluded from every ideal/non-ideal statistic via
    `src.calibrate.filter_single_centre_library()`, so no downstream
    consumer needs a yes/no tag for them either.
    """
    mol_type = dict(zip(roster["molecule"], roster["mol_type"]))
    df = df.copy()
    internal_mask = df["kind"] == "internal"
    df.loc[internal_mask, "ideal"] = df.loc[internal_mask, "molecule"].map(
        lambda m: _MOL_TYPE_TO_IDEAL.get(mol_type.get(m)))
    return df


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
    mol_type_by_molecule = dict(zip(roster["molecule"], roster.get("mol_type", pd.Series(dtype=str))))

    all_rows = []
    load_errors = []
    for _, r in roster.iterrows():
        molecule, base = r["molecule"], r["basename"]
        mol_type = mol_type_by_molecule.get(molecule)
        try:
            rows = score_geometry_molecule(base, data_dir, thresholds, mol_type=mol_type)
        except (FileNotFoundError, ValueError) as e:
            load_errors.append((molecule, base, str(e)))
            continue
        for row in rows:
            row["molecule"] = molecule
        all_rows.extend(rows)

    df = pd.DataFrame(all_rows, columns=SCHEMA_COLUMNS)
    df, skip_report = attach_labels(df, csv_tables, label_lookup)
    df = attach_ideal_tags(df, roster)

    for molecule, base, err in load_errors:
        warnings.warn(
            f"library_ingest: skipped roster molecule '{molecule}' (basename "
            f"'{base}') -- parse/scoring failed ({err}) -- not included in "
            "library_scores.csv at all.",
            stacklevel=3)
    for entry in skip_report:
        mode_index, eng_f, cm_f = entry["example"]
        cm_f_str = f"{cm_f:.4f}" if cm_f is not None else "NO characterised_modes.csv ROW FOUND"
        warnings.warn(
            f"library_ingest: '{entry['molecule']}' internal rows NOT label-"
            f"joined ({entry['n_mismatched']} mode(s) mismatched, e.g. mode "
            f"{mode_index}: engine {eng_f:.4f} vs characterised_modes.csv {cm_f_str} cm-1) -- "
            "this molecule's characterised_modes.csv row likely came from a different "
            "calculation than data/logs/ (e.g. a fallback level of theory -- see "
            "mol_list_method.csv's current_method column), or characterised_modes.csv "
            "needs to be regenerated (regenerate_characterised_modes()). Its scores "
            "(V_Stretch, Tx..Rz, predicted_label, ...) are still the real "
            "engine's own and are NOT affected; only ref_label/ref_key "
            "are left null (ideal is unaffected -- it comes from mol_list_method.csv's "
            "mol_type, not this gate).",
            stacklevel=3)

    # d_CA reporting: warn once per ideal/non-ideal molecule with NO unique
    # central atom found (multi-centre molecules are never attempted at
    # all -- see score_geometry_molecule()'s mol_type gate -- so they are
    # deliberately excluded from this check).
    for molecule, mol_type in mol_type_by_molecule.items():
        if mol_type not in ("ideal", "non-ideal"):
            continue
        sub = df[(df["molecule"] == molecule) & (df["kind"] == "internal")]
        if len(sub) and sub["d_CA"].isna().all():
            warnings.warn(
                f"library_ingest: '{molecule}' (mol_type='{mol_type}') has no unique "
                "central/hub atom (no single atom bonded to every other atom) -- "
                "d_CA left null for every internal row.",
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


# --- Phase 2 (2026-07-08): resync characterised_modes.csv (data_score.csv,
# too, until it was deleted from disk 2026-07-09) freq/mu/k from the on-disk
# log, fixing attach_labels()'s stale-freq gate failures at the source. See
# IMPLEMENTATION_PLAN.md RESUME HERE + the plan file
# `before-that-the-program-ethereal-penguin.md`, Phase 2.
#
# `irrep` is deliberately NOT touched by this resync (see
# `resync_reference_metadata`'s docstring "irrep is out of scope" note) --
# confirmed via `plot_irrep_coupling` (src/figures.py) that
# characterised_modes.csv's irrep column is read downstream with hardcoded
# Unicode-subscript/prime strings (e.g. "A₂\"", "E'", "B₂") matched by exact
# equality; the raw Gaussian-log irrep token is a different, ASCII-only
# alphabet ("A2\"", "E", "B2") that would silently break every category
# match if written in place. Building a correct, fully-general
# ASCII->Unicode-subscript/prime/Greek irrep translator (Sigma/Pi for linear
# groups, primes for D3h, etc.) is a real but separable piece of work
# matching the already-recorded, author-approved deferral ("defer
# regenerating the irrep/bond-length columns ... to a follow-up session") --
# left untouched here rather than risk a silent mistranslation. A second,
# independent reason: at least one molecule (AlCl3-class, Gaussian's own
# near-degenerate "?A"/"?B" placeholder irreps) has an EXISTING irrep value
# that is the author's own manual resolution of an ambiguity Gaussian itself
# could not resolve -- the raw log value is strictly less informative there,
# so overwriting would be a regression, not a resync, even before the
# formatting problem.
#
# `k` (originally verified against the now-deleted data_score.csv) WAS
# independently confirmed to be a real force constant column, not an
# unrelated scoring metric: for every molecule whose freq already agreed
# with the on-disk log (i.e. no staleness), those `k` values matched
# `Force constants ---` in the log to 4 decimal places exactly (e.g. SnH4
# mode 1: log 0.3149 == 0.3149; TeH4 mode 6: log 1.6307 == 1.6307).
# `characterised_modes.csv` has an explicit, unambiguous `k` (force
# constant) + `μ` (reduced mass) column pair with the same semantics,
# safe to overwrite in place.
#
# 2026-07-09 (OH4/OF4 exclusion session): `data/data_score.csv` was deleted
# from disk entirely (an unrelated cleanup of an already-retired file, see
# the module docstring above). This function is updated to resync
# `characterised_modes.csv` only -- it no longer reads or writes
# `data_score.csv`.

def _resolve_log_path(base, data_dir="data"):
    """basename -> full path to its .log/.out file, or None if neither exists."""
    for ext in (".log", ".out"):
        candidate = os.path.join(data_dir, "logs", base + ext)
        if os.path.exists(candidate):
            return candidate
    return None


def _fmt_trim(value, dp=4):
    """Format `value` to `dp` decimal places, then strip trailing zeros (and
    a trailing bare '.') -- matches data_score.csv/characterised_modes.csv's
    own existing convention (e.g. Gaussian's '485.3180' is stored there as
    '485.318', not '485.3180'), so a resync doesn't introduce a purely
    cosmetic reformatting diff on top of the real value change."""
    s = f"{value:.{dp}f}"
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return s if s else "0"


def resync_reference_metadata(data_dir="data", write=True):
    """Resync `freq` (+ `k` force constant, + `μ` reduced mass) in
    `characterised_modes.csv` from the on-disk Gaussian `.log` the engine
    actually scores, for every `mol_list_method.csv` roster molecule that
    already has existing rows in that CSV. `irrep` is intentionally left
    untouched -- see the module comment above this function for why.

    Formerly (2026-07-08) also resynced `data/data_score.csv` in the same
    pass; that file was deleted from disk entirely 2026-07-09 (an unrelated
    cleanup of an already-retired file, see the module docstring), so this
    function now touches `characterised_modes.csv` only. The report/changes
    schema is unchanged (`changes` entries still carry a `csv` field) since
    every entry has always been, and remains, `"characterised_modes.csv"`.

    Rows are matched by 1-based internal mode index (the `mode` column),
    identical to the indexing `main.py build_scorer_and_final()` assigns as
    "Vib i" labels and `attach_labels()` gates on -- i.e. position within
    `GaussianParser.parse(parse_modes=True)['modes']` (Gaussian's own
    frequency-block order; Gaussian already projects out translation/
    rotation, so this list is exactly the 3N-6/3N-5 internal modes with no
    offset bookkeeping needed).

    A molecule's mode COUNT (or the max `mode` index) disagreeing with the
    on-disk log's parsed mode count is treated as a real structural mismatch
    (different geometry/atom count), not a staleness problem, and is
    skipped (that molecule's rows left untouched) -- reported, never
    silently guessed.

    Returns a report dict:
      ``resynced``: [{"molecule", "basename", "n_modes",
                       "changes": [{"csv", "mode", "field", "old", "new"}]}]
        -- one entry per roster molecule actually touched; `changes` lists
        only cells whose VALUE materially changed (same-value overwrites are
        not reported as changes, though the cell is still rewritten).
      ``skipped_mode_count_mismatch``: [{"molecule", "basename", "detail"}]
      ``skipped_no_log``: [{"molecule", "basename", "detail"}]
      ``skipped_no_rows``: [molecule, ...] -- roster molecules with NO row
        in `characterised_modes.csv` at all (the separately-tracked,
        still-open TODO: back-filling literature labels for any remaining
        unlabeled molecules -- resync cannot invent rows that were never
        entered).
    """
    roster = load_mol_roster(data_dir)
    cm_path = os.path.join(data_dir, "characterised_modes.csv")
    cm = pd.read_csv(cm_path, dtype=str, keep_default_na=False)

    report = {
        "resynced": [],
        "skipped_mode_count_mismatch": [],
        "skipped_no_log": [],
        "skipped_no_rows": [],
    }

    for _, r in roster.iterrows():
        molecule, base = r["molecule"], r["basename"]
        cm_idx = list(cm.index[cm["molecule"] == molecule])
        if not cm_idx:
            report["skipped_no_rows"].append(molecule)
            continue

        log_path = _resolve_log_path(base, data_dir)
        if log_path is None:
            report["skipped_no_log"].append(
                {"molecule": molecule, "basename": base, "detail": "no .log/.out on disk"})
            continue
        try:
            engine_modes = GaussianParser(log_path).parse(parse_modes=True)["modes"]
        except Exception as e:
            report["skipped_no_log"].append(
                {"molecule": molecule, "basename": base, "detail": f"parse failed: {e}"})
            continue
        n_engine = len(engine_modes)

        # Validate EXACT agreement (not just count/max -- a formula-auditor
        # review, 2026-07-08, caught that count+max alone would silently
        # pass a CSV with duplicate/non-contiguous mode indices, e.g.
        # [1, 1, 3, 4, 5] for n_engine=5 -- len==5, max==5, but mode 2 is
        # missing and mode 1 is duplicated) between the CSV's `mode` values
        # and {1, ..., n_engine}. A blank "mode" cell is its OWN reported
        # problem (not silently excluded from the count, which could mask a
        # real mismatch and would otherwise crash `int(float(''))` in the
        # update loop below) -- filtered rows are only ever used past this
        # point once the whole molecule has cleanly passed.
        problems = []
        blank_rows = [i for i in cm_idx if cm.at[i, "mode"] == ""]
        valid_idx = [i for i in cm_idx if cm.at[i, "mode"] != ""]
        modes_in_csv = sorted(int(float(cm.at[i, "mode"])) for i in valid_idx)
        if blank_rows:
            problems.append(f"characterised_modes.csv: {len(blank_rows)} row(s) with a blank 'mode' field")
        if modes_in_csv != list(range(1, n_engine + 1)):
            problems.append(
                f"characterised_modes.csv: mode indices {modes_in_csv} do not exactly match "
                f"engine's 1..{n_engine} ({n_engine} parsed modes)")
        cm_idx = valid_idx
        if problems:
            report["skipped_mode_count_mismatch"].append(
                {"molecule": molecule, "basename": base, "detail": "; ".join(problems)})
            continue

        changes = []
        for i in cm_idx:
            mode_i = int(float(cm.at[i, "mode"]))
            m = engine_modes[mode_i - 1]
            new_freq = _fmt_trim(m["frequency"])
            if cm.at[i, "freq"] != new_freq:
                changes.append({"csv": "characterised_modes.csv", "mode": mode_i, "field": "freq",
                                 "old": cm.at[i, "freq"], "new": new_freq})
            cm.at[i, "freq"] = new_freq
            if m.get("force_constant") is not None:
                new_k = _fmt_trim(m["force_constant"])
                if cm.at[i, "k"] != new_k:
                    changes.append({"csv": "characterised_modes.csv", "mode": mode_i, "field": "k",
                                     "old": cm.at[i, "k"], "new": new_k})
                cm.at[i, "k"] = new_k
            if m.get("reduced_mass") is not None:
                new_mu = _fmt_trim(m["reduced_mass"])
                if cm.at[i, "μ"] != new_mu:
                    changes.append({"csv": "characterised_modes.csv", "mode": mode_i, "field": "mu",
                                     "old": cm.at[i, "μ"], "new": new_mu})
                cm.at[i, "μ"] = new_mu

        report["resynced"].append({
            "molecule": molecule, "basename": base, "n_modes": n_engine, "changes": changes,
        })

    if write:
        cm.to_csv(cm_path, index=False, encoding="utf-8-sig")

    return report


# --- Phase 3 (2026-07-08): regenerate characterised_modes.csv from a direct ---
# disk scan, and retire data_score.csv from the code path entirely. See
# IMPLEMENTATION_PLAN.md RESUME HERE for the full author directive. Distinct
# from resync_reference_metadata() above: that function only ever UPDATES
# freq/k/mu on EXISTING rows for the (now-frozen) mol_list_method.csv
# roster; this one determines the row SET itself directly from disk (any
# .log+.gjf/.com pair under data/logs+data/gjf, not restricted to the
# roster), so a molecule dropped into data/logs+data/gjf later is picked up
# automatically without a mol_list_method.csv edit first.

def _basename_to_molecule_map(data_dir="data"):
    """basename -> canonical mol_list_method.csv molecule name (roster-only
    lookup; empty dict if the roster can't be read for any reason). Used by
    regenerate_characterised_modes() to translate an on-disk basename (e.g.
    'SbH3') to the canonical short name ('SbH3') that
    characterised_modes.csv's existing rows are keyed by -- a basename with
    no roster row falls back to using the raw basename itself (no canonical
    alternative exists for it)."""
    try:
        roster = load_mol_roster(data_dir)
    except Exception:
        return {}
    return dict(zip(roster["basename"], roster["molecule"]))


# The 7 manually-curated columns preserved verbatim across a regeneration
# (never engine-derived); "shape" is included here per the author's
# explicit 2026-07-08 decision -- unlike freq/mu/k/irrep, a molecule's
# shape label is not literally reparsed from the log, it is an
# author-assigned classification. Order matches characterised_modes.csv's
# own header so a freshly-built row dict lines up 1:1 with old rows.
_CHARACTERISED_MODES_MANUAL_COLUMNS = (
    "shape", "type", "sym", "description", "νₖ", "ref", "Note",
)


def regenerate_characterised_modes(data_dir="data", write=True):
    """Regenerate data/characterised_modes.csv's engine-derivable columns
    (molecule, mode, freq, mu (`μ`), k) from a DIRECT ON-DISK SCAN of
    data/logs + data/gjf (via discover_geometry_molecules() -- every
    .log/.out + .com/.gjf basename pair actually present, NOT restricted to
    mol_list_method.csv's 72-row roster, so a file pair dropped in later is
    picked up automatically without a roster edit first).

    Preserves the author's manually-curated columns
    (_CHARACTERISED_MODES_MANUAL_COLUMNS: shape, type, sym, description,
    ref, Note, and the literature mode label νₖ) for every (molecule, mode)
    that already had a row; leaves them blank for genuinely new
    (molecule, mode) pairs that never had one.

    `irrep` is DELIBERATELY carved out of the "auto columns overwritten"
    set and treated like a manual column for any row that ALREADY has one:
    an existing row's irrep is preserved untouched; only a genuinely NEW
    row (nothing to lose) gets the raw engine-parsed token. This mirrors
    resync_reference_metadata()'s own documented reasoning (see the module
    comment above that function) -- the on-disk log's irrep token is plain
    ASCII ("B2", "A1'"-with-ASCII-prime, "?A"/"?B" near-degeneracy
    placeholders), while existing characterised_modes.csv values are
    hand-verified Unicode-subscript/prime strings ("B₂", "A₁′") that
    src/figures.py::plot_irrep_coupling matches by exact string equality;
    blindly overwriting every row's irrep with the raw ASCII token would
    silently break that matching (and in at least one known case, regress
    an author-resolved Gaussian placeholder to strictly less information).
    Building a general ASCII->Unicode irrep translator remains the
    already-documented, deliberately deferred follow-up (IMPLEMENTATION_
    PLAN.md) -- not attempted here.

    Molecule naming: see _basename_to_molecule_map()'s docstring.

    Returns a report dict:
      'n_disk_basenames': int -- basenames found on disk (data/logs+data/gjf
        intersection).
      'n_old_molecules': int -- distinct molecules in the OLD
        characterised_modes.csv (before this call).
      'n_new_molecules': int -- distinct molecules in the regenerated file.
      'dropped_molecules': sorted [str, ...] -- molecules present in the OLD
        file with NO on-disk log+gjf match anymore (surfaced, never silent
        -- git history preserves the row regardless).
      'added_molecules': sorted [str, ...] -- molecules newly present on
        disk with no OLD row at all (get blank manual columns).
      'parse_failures': [{'molecule', 'basename', 'detail'}, ...] -- on-disk
        pairs that failed to parse (excluded from the output; reported, not
        silently dropped).
      'n_rows_written': int.
      'n_rows_with_preserved_manual_labels': int -- rows whose manual
        columns were carried over from an existing (molecule, mode) match
        (a strict subset of n_rows_written; the rest are genuinely new).
    """
    old_path = os.path.join(data_dir, "characterised_modes.csv")
    old = pd.read_csv(old_path, dtype=str, keep_default_na=False)
    old_by_key = {}
    for _, r in old.iterrows():
        if r["mode"] == "":
            continue
        old_by_key[(r["molecule"], int(float(r["mode"])))] = r
    old_molecules = set(old["molecule"].unique())

    base_to_mol = _basename_to_molecule_map(data_dir)
    bases = discover_geometry_molecules(data_dir)

    new_rows = []
    parse_failures = []
    new_molecules = set()
    n_preserved = 0
    for base in bases:
        molecule = base_to_mol.get(base, base)
        log_path = _resolve_log_path(base, data_dir)
        if log_path is None:
            parse_failures.append({"molecule": molecule, "basename": base,
                                    "detail": "no .log/.out on disk"})
            continue
        try:
            modes = GaussianParser(log_path).parse(parse_modes=True)["modes"]
        except Exception as e:
            parse_failures.append({"molecule": molecule, "basename": base,
                                    "detail": f"parse failed: {e}"})
            continue

        new_molecules.add(molecule)
        for i, m in enumerate(modes, start=1):
            key = (molecule, i)
            old_row = old_by_key.get(key)
            has_prior = old_row is not None
            if has_prior:
                n_preserved += 1
            row = {col: (old_row[col] if has_prior else "")
                   for col in _CHARACTERISED_MODES_MANUAL_COLUMNS}
            row["molecule"] = molecule
            row["mode"] = str(i)
            row["freq"] = _fmt_trim(m["frequency"])
            row["μ"] = _fmt_trim(m["reduced_mass"]) if m.get("reduced_mass") is not None else ""
            row["k"] = _fmt_trim(m["force_constant"]) if m.get("force_constant") is not None else ""
            row["irrep"] = (old_row["irrep"] if has_prior
                             else (m.get("irrep") or ""))
            new_rows.append(row)

    new_df = pd.DataFrame(new_rows, columns=old.columns.tolist())

    report = {
        "n_disk_basenames": len(bases),
        "n_old_molecules": len(old_molecules),
        "n_new_molecules": len(new_molecules),
        "dropped_molecules": sorted(old_molecules - new_molecules),
        "added_molecules": sorted(new_molecules - old_molecules),
        "parse_failures": parse_failures,
        "n_rows_written": len(new_df),
        "n_rows_with_preserved_manual_labels": n_preserved,
    }

    if write:
        new_df.to_csv(old_path, index=False, encoding="utf-8-sig")

    return report
