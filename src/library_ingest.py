"""Library ingest: builds data/results/library_scores.csv for every molecule
in the JCC paper's roster, ``data/mol_list_method.csv`` (77 rows: 10 "ideal"
single-centre AB_n shapes, 57 "non-ideal" substituted variants, 1
"multi-centre" = benzene, 9 "test" = a held-out transferability-test set).
``data/mol_list_method.csv`` is the single source of truth; since 2026-08-11
its ``molecule`` column doubles as the on-disk basename (the old, redundant
``basename`` column was dropped -- every on-disk file was already renamed to
match ``molecule`` in an earlier commit). This module never reads a
precomputed spreadsheet -- every score is a fresh real-engine recompute
(Steps 1-4: ``main.load_inputs`` -> ``main.build_scorer_and_final`` ->
``src.classifier.classify_all_modes``, ``ModeScorer.score_bonds()``).

Pipeline
--------
1. ``load_mol_roster()`` reads the roster.
2. ``check_roster_disk_consistency()`` cross-checks every roster molecule
   name (used directly as the on-disk basename) against what's on disk. A
   roster row with no matching ``.log``+``.gjf`` pair is a fatal
   ``FileNotFoundError`` (100% coverage is expected). An on-disk basename
   with no roster row (e.g. the gramicidin ``1grm_MM_UFF`` companion-paper
   inputs, out of scope here) is an expected orphan: non-fatal warning,
   excluded from the output.
3. ``score_geometry_molecule()`` runs the real engine per molecule. A
   molecule whose files exist but fail to parse/score is warned and skipped
   (excluded from the CSV), not raised -- distinct from "missing from disk".
4. ``attach_labels()`` joins ``ref_label``/``ref_key`` from
   ``src/csv_label_ingest.py``'s CSVs, gated by a per-molecule frequency
   check (see that function's docstring). ``attach_ideal_tags()`` separately
   populates ``ideal`` from the roster's ``mol_type`` column.

``data/data_score.csv`` is retired: every quantity it used to supply now
comes from the engine (``d_CA``) or from ``characterised_modes.csv``/
``mol_list_method.csv``. ``regenerate_characterised_modes()`` keeps
``characterised_modes.csv`` in sync with on-disk logs via a direct disk scan
(not roster-restricted), preserving manually-curated literature columns.

Output schema (``SCHEMA_COLUMNS``, unchanged names so ``src/calibrate.py``/
``src/figures.py`` keep working). ``has_geometry`` is unconditionally
``True`` (every roster molecule has on-disk geometry) but kept in the schema
so downstream readers (e.g. ``src/calibrate.py``'s ``_load_geometry_pool``)
don't need to change.
"""
import os
import warnings

import numpy as np
import pandas as pd

from src import csv_label_ingest
from src.parser import GaussianParser
from src.scoring import format_bond_map, get_v_weighting
from src.utils import find_file

_EXTERNAL_SLOTS = ("Tx", "Ty", "Tz", "Rx", "Ry", "Rz")

# Locked output schema (src/calibrate.py and src/figures.py read these exact
# column names, by name not position). `reduced_mass`/`force_constant` are
# named that way (not bare `k`) to avoid colliding with this CSV's existing
# `k` column, a differently-scoped scoring metric adjacent to `vib_scr`.
SCHEMA_COLUMNS = [
    "molecule", "mode_index", "kind", "freq", "ref_label", "ideal",
    "V_Stretch", "delta_b_mean", "s_AB", "rel_db", "has_geometry",
    "predicted_label", "predicted_annotation",
    "Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "ref_key",
    "reduced_mass", "force_constant", "irrep",
    # d_CA: central/hub-atom displacement amplitude, ||mode_vector[central
    # atom]||, for internal rows of ideal/non-ideal molecules with a unique
    # hub atom (_central_atom_index()); null otherwise (external rows,
    # multi-centre molecules, no-unique-hub cases).
    "d_CA",
    # Which eq:vscore bond weighting produced V_Stretch/s_AB on this row, so
    # the CSV is self-describing about its own scoring definition.
    "v_weighting",
]


def load_mol_roster(data_dir="data"):
    """Read data/mol_list_method.csv: the authoritative molecule roster
    (canonical `molecule` name, which also doubles as the on-disk basename --
    the `basename` column was dropped 2026-08-11 as redundant). Fail loud if
    the file or the required column is missing."""
    path = os.path.join(data_dir, "mol_list_method.csv")
    df = pd.read_csv(path)
    missing_cols = {"molecule"} - set(df.columns)
    if missing_cols:
        raise ValueError(f"mol_list_method.csv missing required column(s): {missing_cols}")
    return df


def load_library_scores(data_dir="data", lib_df=None):
    """Pass through `lib_df` if already loaded, else read
    data/results/library_scores.csv fresh."""
    if lib_df is not None:
        return lib_df
    return pd.read_csv(os.path.join(data_dir, "results", "library_scores.csv"))


def check_roster_disk_consistency(roster, data_dir="data"):
    """Cross-check roster['molecule'] (== on-disk basename) against
    discover_geometry_molecules()'s directory intersection. Returns
    (missing, orphaned):
      missing  -- [(molecule, basename), ...] roster rows whose .log+.gjf pair
                  is NOT present on disk (basename == molecule here).
      orphaned -- sorted list of on-disk basenames not referenced by any
                  roster row's molecule column.
    """
    disk_bases = set(discover_geometry_molecules(data_dir))
    roster_bases = set(roster["molecule"])
    missing = [(m, m) for m in roster["molecule"] if m not in disk_bases]
    orphaned = sorted(disk_bases - roster_bases)
    return missing, orphaned


def resolve_log_basename(molecule, data_dir="data"):
    """molecule (canonical mol_list_method.csv name) -> on-disk basename, or
    None if not in the roster. Kept as its own function (name/signature
    unchanged) since src/calibrate.py calls it directly. Since the
    `basename` column was dropped 2026-08-11, this is now effectively a
    roster-membership check that returns the molecule name itself."""
    roster = load_mol_roster(data_dir)
    matches = roster.loc[roster["molecule"] == molecule, "molecule"]
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
    mode) in the library_scores.csv schema, excluding 'molecule' (caller
    attaches it) with ref_label/ideal/ref_key left None (attached separately).

    `mol_type` gates `d_CA`: computed only for 'ideal'/'non-ideal' molecules
    with a unique hub atom (_central_atom_index()); None/'multi-centre' never
    computes it.

    Raises ValueError if no bond connectivity is available (propagated from
    build_scorer_and_final) -- a real data problem, not silently skipped
    here; the caller decides whether to skip-and-warn.
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
                "v_weighting": scorer.v_weighting,
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
        s_ab_str = format_bond_map(bonds, "s_AB")
        rel_db_str = format_bond_map(bonds, "rel_db")
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
            "v_weighting": scorer.v_weighting,
        })
    return rows


def multi_centre_molecules(data_dir="data"):
    """Molecule names tagged mol_type=='multi-centre' in mol_list_method.csv
    -- the authoritative multi-centre / no-single-hub-atom classification.
    Of the historical 7-name exclusion set, only C6H6 is actually in the
    current 77-molecule roster, so this is the single roster-driven source
    now (forward compatible with any future multi-centre addition, no code
    edit needed). Still used for figure filtering; src/calibrate.py's
    SINGLE_CENTRE_ONLY_EXCLUDE is now driven by the broader
    out_of_calibration_scope_molecules() instead (see that function).
    """
    roster = load_mol_roster(data_dir)
    if "mol_type" not in roster.columns:
        raise ValueError("mol_list_method.csv is missing the 'mol_type' column")
    return frozenset(roster.loc[roster["mol_type"] == "multi-centre", "molecule"])


def out_of_calibration_scope_molecules(data_dir="data"):
    """Molecule names OUTSIDE the stretch/bend calibration's inclusion scope
    (mol_type not in {'ideal', 'non-ideal'}) -- currently the multi-centre
    molecule (C6H6) plus the 9 held-out 'test' transferability molecules.
    Drives src/calibrate.py's SINGLE_CENTRE_ONLY_EXCLUDE: an INCLUSION
    filter (ideal/non-ideal only) expressed as its complement, so a future
    mol_type category is automatically excluded from calibration with no
    hardcoded special-case (2026-08-11 decision, replacing the old
    multi-centre-only exclusion list now that mol_type=='test' exists)."""
    roster = load_mol_roster(data_dir)
    return frozenset(roster.loc[~roster["mol_type"].isin(("ideal", "non-ideal")), "molecule"])


def attach_labels(df, csv_tables, label_lookup, freq_atol=0.05, freq_rtol=1e-4):
    """Join ref_label/ref_key onto `df`'s internal rows, molecule by
    molecule, gated by a whole-molecule frequency-agreement check against
    `csv_tables["characterised_modes"]`. This gate is NOT dead code: a
    single mismatched or missing mode disqualifies the WHOLE molecule's
    label join (no half-merge), catching real mode-index mismatches between
    the engine's parsed frequency and characterised_modes.csv's independent
    expectation. Scores themselves are never affected by a failed join --
    only ref_label/ref_key are left null. `ideal` is set separately by
    `attach_ideal_tags()` (a structural roster property, ungated).

    Returns (df, skip_report): skip_report is a list of {'molecule',
    'n_mismatched', 'example': (mode_index, engine_freq, cm_freq)} dicts, one
    per molecule whose join was skipped. A molecule entirely absent from
    characterised_modes.csv is left untouched with no skip-report entry --
    not an error, just no ground truth to gate against.
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
        molecule = r["molecule"]
        base = molecule  # 2026-08-11: `basename` column dropped; molecule == on-disk basename
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

    # Warn once per ideal/non-ideal molecule with no unique central atom.
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


def _thresholds_for_current_weighting():
    """Calibrated thresholds if they match the active V-score weighting, else
    provisional bootstrap ones with a warning.

    Breaks the chicken-and-egg after a variant switch: this pass has to score
    the library to produce the distribution the new tau_S/tau_B are derived
    FROM, but the thresholds.json on disk is still stamped for the old
    definition and classify_all_modes() would refuse it. Only the
    predicted_label column depends on thresholds -- V_Stretch, which drives the
    derivation, does not -- so a provisional pass is sound, and --library is
    re-run after --calibrate to make the labels consistent.
    """
    from src.classifier import Thresholds

    thresholds = Thresholds.calibrated()
    active = get_v_weighting()
    if thresholds.v_weighting not in ("*", active):
        warnings.warn(
            f"thresholds.json is calibrated for v_weighting="
            f"{thresholds.v_weighting!r} but scoring is running under "
            f"{active!r}; using provisional bootstrap thresholds for this pass. "
            f"Run --calibrate next, then --library again so predicted_label "
            f"reflects the recalibrated cut points.")
        return Thresholds.bootstrap()
    return thresholds


def build_library_scores(data_dir="data", thresholds=None, return_skip_report=False):
    """Build the library_scores DataFrame: real-engine recompute for every
    one of the 77 mol_list_method.csv roster molecules. Returns a DataFrame
    with exactly SCHEMA_COLUMNS.

    Raises FileNotFoundError if any roster row's .log/.gjf pair is missing
    from disk (see _build_library_scores / check_roster_disk_consistency).
    """
    if thresholds is None:
        thresholds = _thresholds_for_current_weighting()
    return _build_library_scores(data_dir, thresholds, return_skip_report)


def run_ingest_pipeline(data_dir="data", thresholds=None, write=True):
    """Headless entry point: build library_scores and optionally write the
    CSV."""
    df, skip_report = build_library_scores(data_dir, thresholds, return_skip_report=True)
    out_path = os.path.join(data_dir, "results", "library_scores.csv")
    if write:
        df.to_csv(out_path, index=False)
    return df, out_path, skip_report


# resync_reference_metadata() below resyncs characterised_modes.csv's
# freq/mu/k from the on-disk log, fixing attach_labels()'s stale-freq gate
# failures at the source. `irrep` is deliberately excluded from this resync:
# see regenerate_characterised_modes()'s docstring for the full ASCII-vs-
# Unicode irrep carve-out reasoning (authoritative site).
#
# `k` was independently cross-checked as a genuine force constant, not an
# unrelated scoring metric: for molecules whose freq already agreed with the
# log, `k` matched the log's "Force constants ---" line to 4 dp exactly
# (e.g. SnH4 mode 1: 0.3149 == 0.3149; TeH4 mode 6: 1.6307 == 1.6307).

def _resolve_log_path(base, data_dir="data"):
    """basename -> full path to its .log/.out file, or None if neither exists."""
    return find_file(os.path.join(data_dir, "logs"), base, (".log", ".out"))


def _fmt_trim(value, dp=4):
    """Format `value` to `dp` dp, stripping trailing zeros/'.' -- matches
    characterised_modes.csv's existing convention (e.g. '485.3180' stored as
    '485.318'), avoiding a purely cosmetic reformatting diff."""
    s = f"{value:.{dp}f}"
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return s if s else "0"


def resync_reference_metadata(data_dir="data", write=True):
    """Resync `freq`/`k`/`mu` in `characterised_modes.csv` from the on-disk
    Gaussian `.log` for every roster molecule with existing rows there.
    `irrep` is left untouched (see regenerate_characterised_modes()).

    Rows are matched by 1-based internal mode index (`mode` column),
    matching `main.py build_scorer_and_final()`'s "Vib i" numbering --
    position within `GaussianParser.parse(parse_modes=True)['modes']`.

    A molecule's mode count disagreeing with the on-disk log's parsed count
    is a real structural mismatch (not staleness); skipped and reported,
    never silently guessed.

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
        molecule = r["molecule"]
        base = molecule  # 2026-08-11: `basename` column dropped; molecule == on-disk basename
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

        # Validate EXACT agreement (not just count/max): count+max alone
        # would silently pass duplicate/non-contiguous mode indices, e.g.
        # [1, 1, 3, 4, 5] for n_engine=5 (mode 2 missing, mode 1 duplicated).
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


# regenerate_characterised_modes() below determines the row SET directly
# from disk (unlike resync_reference_metadata(), which only updates existing
# rows for the fixed roster), so a molecule dropped into data/logs+data/gjf
# later is picked up automatically without a mol_list_method.csv edit.

def _basename_to_molecule_map(data_dir="data"):
    """basename -> canonical roster molecule name (empty dict if roster
    unreadable); falls back to the raw basename if no roster row matches.
    Since the `basename` column was dropped 2026-08-11 (molecule == on-disk
    basename now), this is an identity map over the roster's molecule
    names."""
    try:
        roster = load_mol_roster(data_dir)
    except Exception:
        return {}
    return dict(zip(roster["molecule"], roster["molecule"]))


# The 7 manually-curated columns preserved verbatim across a regeneration
# (never engine-derived; "shape" is an author-assigned classification, not
# reparsed from the log). Order matches characterised_modes.csv's header.
_CHARACTERISED_MODES_MANUAL_COLUMNS = (
    "shape", "type", "sym", "description", "νₖ", "ref", "Note",
)


def regenerate_characterised_modes(data_dir="data", write=True, force_refresh_irrep=False):
    """Regenerate data/characterised_modes.csv's engine-derivable columns
    (molecule, mode, freq, mu, k) from a direct on-disk scan of data/logs +
    data/gjf (not restricted to the roster, so new file pairs are picked up
    automatically). Preserves manually-curated columns
    (_CHARACTERISED_MODES_MANUAL_COLUMNS) for existing (molecule, mode) rows;
    blank for genuinely new ones.

    `irrep` carve-out (authoritative reasoning; referenced elsewhere in this
    module and in src/figures.py): by default (`force_refresh_irrep=False`),
    an existing row's irrep is preserved untouched, never overwritten with
    the raw engine-parsed token. The on-disk log's irrep is plain ASCII
    ("B2", "A1'"), while existing characterised_modes.csv values are
    hand-verified Unicode-subscript/prime strings ("B₂", "A₁′") that
    src/figures.py::plot_irrep_coupling matches by exact string equality --
    e.g. BBr3/OCl2 are known cases with this ASCII-vs-Unicode mismatch.
    Overwriting would silently break that matching (and regress at least one
    author-resolved Gaussian near-degeneracy placeholder to less
    information). A general ASCII->Unicode translator is deliberately
    deferred, not attempted here. Only genuinely new rows (nothing to lose)
    get the raw engine-parsed irrep token by default.

    `force_refresh_irrep=True` overrides the carve-out: every row's irrep
    (existing or new) is overwritten with the fresh engine-parsed token from
    THIS call's own `data_dir` logs. This intentionally regresses the
    Unicode/hand-verified formatting described above -- it exists for
    exploratory/consistency-check runs against an alternate `data_dir` (e.g.
    a different Gaussian-version rerun mirror), where the whole point is to
    see exactly what that run's own logs say, not what a prior curation
    recorded. Never pass True when regenerating the canonical
    data/characterised_modes.csv in place.

    Returns a report dict:
      'n_disk_basenames', 'n_old_molecules', 'n_new_molecules': counts.
      'dropped_molecules' / 'added_molecules': sorted molecule-name lists
        (no on-disk match anymore / newly present on disk).
      'parse_failures': [{'molecule', 'basename', 'detail'}, ...].
      'n_rows_written', 'n_rows_with_preserved_manual_labels': int.
      'n_rows_irrep_overridden': int -- rows whose irrep was force-refreshed
        AND actually changed value (0 whenever force_refresh_irrep=False).
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
    n_irrep_overridden = 0
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
            engine_irrep = m.get("irrep") or ""
            if has_prior and not force_refresh_irrep:
                row["irrep"] = old_row["irrep"]
            else:
                row["irrep"] = engine_irrep
                if has_prior and old_row["irrep"] != engine_irrep:
                    n_irrep_overridden += 1
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
        "n_rows_irrep_overridden": n_irrep_overridden,
    }

    if write:
        new_df.to_csv(old_path, index=False, encoding="utf-8-sig")

    return report
