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
                # construct_T/construct_R never set these (synthetic, not
                # parsed from a Gaussian frequency block) -- None by design.
                "reduced_mass": m.get("reduced_mass"),
                "force_constant": m.get("force_constant"),
                "irrep": m.get("irrep"),
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
            "reduced_mass": m.get("reduced_mass"),
            "force_constant": m.get("force_constant"),
            "irrep": m.get("irrep"),
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


# --- Phase 2 (2026-07-08): resync data_score.csv/characterised_modes.csv ---
# freq/mu/k from the on-disk log, fixing attach_labels()'s stale-freq gate
# failures at the source. See IMPLEMENTATION_PLAN.md RESUME HERE + the plan
# file `before-that-the-program-ethereal-penguin.md`, Phase 2.
#
# `irrep` is deliberately NOT touched by this resync (see
# `resync_reference_metadata`'s docstring "irrep is out of scope" note) --
# confirmed via `plot_irrep_coupling` (src/figures.py) that data_score.csv's
# irrep column is read downstream with hardcoded Unicode-subscript/prime
# strings (e.g. "A₂\"", "E'", "B₂") matched by exact equality; the raw
# Gaussian-log irrep token is a different, ASCII-only alphabet ("A2\"", "E",
# "B2") that would silently break every category match if written in place.
# Building a correct, fully-general ASCII->Unicode-subscript/prime/Greek
# irrep translator (Sigma/Pi for linear groups, primes for D3h, etc.) is a
# real but separable piece of work matching the already-recorded,
# author-approved deferral ("defer regenerating data/data_score.csv's
# irrep/bond-length columns ... to a follow-up session") -- left untouched
# here rather than risk a silent mistranslation. A second, independent
# reason: at least one molecule (AlCl3-class, Gaussian's own near-degenerate
# "?A"/"?B" placeholder irreps) has an EXISTING data_score.csv irrep value
# that is the author's own manual resolution of an ambiguity Gaussian itself
# could not resolve -- the raw log value is strictly less informative there,
# so overwriting would be a regression, not a resync, even before the
# formatting problem.
#
# `k` (data_score.csv) WAS independently verified this session (not assumed)
# to be a real force constant column, not an unrelated scoring metric: for
# every molecule whose freq already agrees with the on-disk log (i.e. no
# staleness), data_score.csv's `k` values match `Force constants ---` in the
# log to 4 decimal places exactly (e.g. SnH4 mode 1: log 0.3149 == ds 0.3149;
# TeH4 mode 6: log 1.6307 == ds 1.6307). `characterised_modes.csv` already
# has an explicit, unambiguous `k` (force constant) + `μ` (reduced mass)
# column pair with the same semantics. Both are safe to overwrite in place.

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
    """Resync `freq` (+ `k` force constant, + `characterised_modes.csv`'s
    `μ` reduced mass) in `data_score.csv`/`characterised_modes.csv` from the
    on-disk Gaussian `.log` the engine actually scores, for every
    `mol_list_method.csv` roster molecule that already has existing rows in
    either CSV. `irrep` is intentionally left untouched -- see the module
    comment above this function for why.

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
    skipped (both CSVs' rows for that molecule left untouched) -- reported,
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
        in EITHER csv at all (the separately-tracked, still-open TODO:
        back-filling literature labels for the 23 T-shaped/see-saw/BBr3/OCl2
        molecules -- resync cannot invent rows that were never entered).
    """
    roster = load_mol_roster(data_dir)
    ds_path = os.path.join(data_dir, "data_score.csv")
    cm_path = os.path.join(data_dir, "characterised_modes.csv")
    ds = pd.read_csv(ds_path, dtype=str, keep_default_na=False)
    cm = pd.read_csv(cm_path, dtype=str, keep_default_na=False)

    report = {
        "resynced": [],
        "skipped_mode_count_mismatch": [],
        "skipped_no_log": [],
        "skipped_no_rows": [],
    }

    for _, r in roster.iterrows():
        molecule, base = r["molecule"], r["basename"]
        ds_idx = list(ds.index[ds["molecule"] == molecule])
        cm_idx = list(cm.index[cm["molecule"] == molecule])
        if not ds_idx and not cm_idx:
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
        filtered_idx = {}
        for label, idx, table in (("data_score.csv", ds_idx, ds),
                                   ("characterised_modes.csv", cm_idx, cm)):
            if not idx:
                filtered_idx[label] = []
                continue
            blank_rows = [i for i in idx if table.at[i, "mode"] == ""]
            valid_idx = [i for i in idx if table.at[i, "mode"] != ""]
            modes_in_csv = sorted(int(float(table.at[i, "mode"])) for i in valid_idx)
            if blank_rows:
                problems.append(f"{label}: {len(blank_rows)} row(s) with a blank 'mode' field")
            if modes_in_csv != list(range(1, n_engine + 1)):
                problems.append(
                    f"{label}: mode indices {modes_in_csv} do not exactly match "
                    f"engine's 1..{n_engine} ({n_engine} parsed modes)")
            filtered_idx[label] = valid_idx
        if problems:
            report["skipped_mode_count_mismatch"].append(
                {"molecule": molecule, "basename": base, "detail": "; ".join(problems)})
            continue
        ds_idx = filtered_idx["data_score.csv"]
        cm_idx = filtered_idx["characterised_modes.csv"]

        changes = []
        for i in ds_idx:
            mode_i = int(float(ds.at[i, "mode"]))
            m = engine_modes[mode_i - 1]
            new_freq = _fmt_trim(m["frequency"])
            if ds.at[i, "freq"] != new_freq:
                changes.append({"csv": "data_score.csv", "mode": mode_i, "field": "freq",
                                 "old": ds.at[i, "freq"], "new": new_freq})
            ds.at[i, "freq"] = new_freq
            if m.get("force_constant") is not None:
                new_k = _fmt_trim(m["force_constant"])
                if ds.at[i, "k"] != new_k:
                    changes.append({"csv": "data_score.csv", "mode": mode_i, "field": "k",
                                     "old": ds.at[i, "k"], "new": new_k})
                ds.at[i, "k"] = new_k

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
        ds.to_csv(ds_path, index=False, encoding="utf-8-sig")
        cm.to_csv(cm_path, index=False, encoding="utf-8-sig")

    return report
