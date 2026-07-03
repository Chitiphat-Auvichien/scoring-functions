"""Library ingest: data/logs/+data/gjf/ (real engine) -> data/results/library_scores.csv,
with data/vibrational-scoring-functions.xlsx used ONLY as a ground-truth label
lookup (``ref_label``/``ideal``).

**Architecture (rearchitected 2026-07-03, author directive; supersedes the
earlier "ingest, don't recompute" locked decision that this module's docstring
described until this session -- see IMPLEMENTATION_PLAN.md's RESUME HERE for
the full record).** The author is supplying real Gaussian ``.log``/``.gjf``
pairs for the hydride-library molecule set incrementally over time (25 of
~70 Excel molecules have a pair on disk as of this session). Every score
column is now recomputed from the raw files via the real engine
(``main.load_inputs`` -> ``build_scorer_and_final`` -> ``classify_all_modes``,
``ModeScorer.score_bonds``) for every molecule physically present in
``data/logs/`` ∩ ``data/gjf/`` -- the Excel workbook no longer drives the loop
or supplies any score. It supplies exactly two things, joined on
(molecule, mode index) with a frequency sanity check gating the join (never
gating score computation): the literature ``ref_label`` ('bend'/'stretch',
from ``data_score``'s ``type`` column) and the ``ideal`` tag. The subset of
on-disk molecules is expected to grow toward the full ~70 as the author drops
in more files; **no code changes are needed for that** -- the ingest loop is
keyed off ``discover_geometry_molecules()``, i.e. what is physically present
in ``data/logs/``+``data/gjf/``, not off Excel rows. Excel-only molecules
(no on-disk geometry) no longer appear in ``library_scores.csv`` at all (a
real, disclosed shrinkage of the calibration population from ~70 to 25
molecules until more files arrive -- not a bug; see IMPLEMENTATION_PLAN.md's
RESUME HERE for the resulting N and the validation performed against it).

Per-row provenance (locked schema, unchanged column names so
``src/calibrate.py``/``src/figures.py`` keep working):
  - ``freq``, ``V_Stretch``, ``Tx..Rz``, ``predicted_label``,
    ``predicted_annotation`` -- ALWAYS the real engine's own numbers
    (``classify_all_modes``'s Step 1-4 output). ``has_geometry`` is therefore
    ``True`` for every row now (only geometry-backed molecules are ever
    ingested); the column is kept for schema/backward-compatibility rather
    than dropped.
  - ``s_AB`` -- per-bond eq:bondscore contributions (``ModeScorer.score_bonds()``),
    computed from THIS mode's own displacement, for every internal row
    (not just STRETCHING/MIXED_STRETCH_BEND rows -- ``classify_all_modes()``
    itself only keeps bonds for those two labels, so this module calls
    ``score_bonds()`` again directly rather than reusing its filtered output).
    Formatted identically to before: semicolon-joined ``"i-j:value"`` with
    1-based atom indices.
  - ``rel_db`` -- per-bond SIGNED relative bond-length change,
    ``rel_db_AB = (|b_AB + Δd_AB| - |b_AB|) / |b_AB|``, where ``b_AB`` is the
    equilibrium bond vector and ``Δd_AB = d_B - d_A`` is the mode-displacement
    difference (the SAME convention ``s_AB``/``Vscore`` already use -- see
    ``ModeScorer._bond_contributions()`` in src/scoring.py). This diagnostic
    quantity (used as fig:bondscores' x-axis) is not a named JCC equation; its
    definition was originally reverse-engineered THIS way from the Excel
    workbook's own ``"(|b'|-|b|)/|b|"`` per-bond column (verified there, in an
    earlier session, to 1e-16 self-consistency against Excel's own mode-
    averaged ``Delta|b|`` companion column) -- now computed from geometry
    instead of read from the spreadsheet, same formula.
  - ``delta_b_mean`` -- mean of ``|rel_db_AB|`` over a mode's bonds (the
    mode-level companion to the per-bond ``rel_db`` string, matching the
    Excel ``Delta|b|`` column's own definition).
  - ``ref_label``/``ideal`` -- ONLY from Excel's ``data_score`` sheet, ONLY
    for internal rows, ONLY if the engine's own computed frequency for that
    mode agrees with Excel's ``freq`` for the SAME (molecule, mode index)
    within tolerance (``atol=0.05, rtol=1e-4``) for EVERY internal mode of
    that molecule (checked before writing anything -- one mismatched mode
    disqualifies the whole molecule's internal-row label join, same
    all-or-nothing guarantee as before). External (T/R) rows get
    ``ref_label`` = 'translation'/'rotation' unconditionally -- that is
    structural truth (Eckart-Sayvetz completeness), not Excel-derived.
    A molecule with NO Excel counterpart at all is not an error: its internal
    rows are still fully scored, just with ``ref_label``/``ideal`` left null.

Excel sheet layout used (unchanged from the prior session's direct
inspection): ``data_score`` has one row per (molecule, 1-based internal
vibrational-mode index) with columns ``molecule``, ``mode``, ``freq``,
``type`` ('bend'/'stretch'), ``ideal`` ('yes'/'no'). This module reads ONLY
those columns now (``sum_|d2-d1|^2/...`` V-score/bond columns are no longer
read at all -- they're recomputed, not ingested). ``data_mode&bond`` is no
longer read by this module either (its per-bond columns are superseded by
``ModeScorer.score_bonds()``'s own geometry-based computation).

**Discovered in the prior (Excel-driven) session, still true and still
handled the same way under this architecture:** water's Excel 'H2O' row does
NOT match this repo's ``water.log`` (Excel freqs 1722.454/3501.541/3660.830
cm-1 vs. the engine's 1628.029/3887.192/4005.506 cm-1 -- 94-386 cm-1 off, not
rounding) -- likely a different calculation/basis than the log used
elsewhere in this repo. Water's 3 internal rows are correctly left with null
``ref_label``/``ideal`` (frequency mismatch -> skipped label join, warned),
while water's 6 (engine-computed, always-correct) external T/R rows are
unaffected. Same story for OF2/Cl2O/Br2O (also historically mismatched).

EXCEL_TO_LOG maps an Excel ``data_score`` molecule name to its on-disk
``data/logs/``+``data/gjf/`` basename for the subset that ships under a
different filename/casing convention; every other on-disk basename is tried
as a direct-name match against ``data_score``'s own molecule column (see
``resolve_excel_molecule_name``). ``resolve_log_basename`` (the reverse
direction: Excel/library molecule name -> on-disk basename) is kept because
``src/calibrate.py`` still needs it to re-load geometry for a given
``library_scores.csv`` molecule name.

Excluded: 'Gly5' (a gramicidin fragment; Decision 5, deferred to the
companion paper) -- has no on-disk .log/.gjf pair in this repo anyway, so it
is naturally absent from the disk-driven loop; still filtered out of the
Excel side too (belt-and-suspenders) in case that ever changes.
"""
import os
import warnings

import numpy as np
import openpyxl
import pandas as pd

EXCLUDED_MOLECULES = {"Gly5"}  # Decision 5: gramicidin fragment, out of scope this paper.

# Excel molecule name -> data/logs & data/gjf basename, for the subset that
# ships with real Gaussian geometry in this repo under a different filename
# (case or naming convention). Everything else is tried as a direct/case
# match against the molecule name itself (see resolve_log_basename /
# resolve_excel_molecule_name).
EXCEL_TO_LOG = {
    "H2O": "water",
    "C6H6": "benzene",
    "CO2": "co2_mp2_3-21g",
    "Cl2O": "ocl2",
    "OF2": "of2",
    "Br2O": "br2o",
    "SeBr2": "SeBr2-cc",
}
_LOG_TO_EXCEL = {v: k for k, v in EXCEL_TO_LOG.items()}

_EXTERNAL_SLOTS = ("Tx", "Ty", "Tz", "Rx", "Ry", "Rz")

# Locked output schema (src/calibrate.py and src/figures.py read these exact
# column names).
SCHEMA_COLUMNS = [
    "molecule", "mode_index", "kind", "freq", "ref_label", "ideal",
    "V_Stretch", "delta_b_mean", "s_AB", "rel_db", "has_geometry",
    "predicted_label", "predicted_annotation",
    "Tx", "Ty", "Tz", "Rx", "Ry", "Rz",
]


def _read_sheet(xlsx_path, sheet_name):
    wb = openpyxl.load_workbook(xlsx_path, data_only=True)
    ws = wb[sheet_name]
    rows = list(ws.iter_rows(values_only=True))
    header, body = rows[0], rows[1:]
    return pd.DataFrame(body, columns=header)


def load_excel_tables(xlsx_path):
    """Read the data_score sheet this module needs (label lookup only)."""
    return {"data_score": _read_sheet(xlsx_path, "data_score")}


def resolve_log_basename(molecule, data_dir="data"):
    """Return the data/logs/<basename> stem for `molecule` if a matching .log
    (or .out) AND a connectivity .com/.gjf pair both exist in this repo;
    otherwise None (Excel-only molecule -- no geometry available here).

    Kept for callers that go Excel-name -> on-disk basename (src/calibrate.py
    still needs this to re-load geometry for a given library_scores.csv
    molecule name).
    """
    logs_dir = os.path.join(data_dir, "logs")
    gjf_dir = os.path.join(data_dir, "gjf")
    candidates = [EXCEL_TO_LOG.get(molecule, molecule), molecule]
    seen = set()
    for base in candidates:
        if base in seen:
            continue
        seen.add(base)
        has_log = any(os.path.exists(os.path.join(logs_dir, base + ext)) for ext in (".log", ".out"))
        has_gjf = any(os.path.exists(os.path.join(gjf_dir, base + ext)) for ext in (".com", ".gjf"))
        if has_log and has_gjf:
            return base
    return None


def discover_geometry_molecules(data_dir="data"):
    """Sorted list of basenames present in BOTH data/logs/ (.log or .out) and
    data/gjf/ (.com or .gjf) -- the disk-driven source of truth for this
    module's molecule loop. Scales automatically as more files are dropped
    in; no code changes needed."""
    logs_dir = os.path.join(data_dir, "logs")
    gjf_dir = os.path.join(data_dir, "gjf")
    log_bases = {os.path.splitext(f)[0] for f in os.listdir(logs_dir)
                 if f.lower().endswith((".log", ".out"))}
    gjf_bases = {os.path.splitext(f)[0] for f in os.listdir(gjf_dir)
                 if f.lower().endswith((".com", ".gjf"))}
    return sorted(log_bases & gjf_bases)


def resolve_excel_molecule_name(base, excel_molecules):
    """On-disk basename -> Excel data_score molecule name, or None if `base`
    has no Excel counterpart at all (not an error -- see module docstring
    point 5)."""
    if base in _LOG_TO_EXCEL:
        return _LOG_TO_EXCEL[base]
    if base in excel_molecules:
        return base
    return None


def score_geometry_molecule(base, data_dir="data", thresholds=None):
    """Run the real engine (Steps 1-4) on one on-disk molecule and return a
    list of row dicts (one per external T/R slot + one per internal 'Vib i'
    mode) in the library_scores.csv schema, EXCLUDING 'molecule' (the caller
    attaches that) and with ref_label/ideal left None (attached separately by
    attach_excel_labels, label-only, per this module's architecture).

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

        def _label(idx):
            # Atom-symbol + 1-based-index label (e.g. "C1", "H7"), matching
            # the convention the old Excel-sourced s_AB strings used (see
            # src/benzene_validation.py's C-C/C-H bond-type parsing, which
            # relies on this exact "<symbol><1-based index>" format).
            return f"{scorer.atoms[idx].symbol}{idx + 1}"

        s_ab_str = ";".join(f"{_label(b['i'])}-{_label(b['j'])}:{b['s_AB']:.4f}" for b in bonds)
        rel_db_str = ";".join(f"{_label(b['i'])}-{_label(b['j'])}:{b['rel_db']:.4f}" for b in bonds)
        delta_b_mean = float(np.mean([abs(b["rel_db"]) for b in bonds])) if bonds else None

        rows.append({
            "mode_index": mode_index, "kind": "internal",
            "freq": m["frequency"],
            "ref_label": None, "ideal": None,  # attached later from Excel, label-only
            "V_Stretch": m["V"], "delta_b_mean": delta_b_mean,
            "s_AB": s_ab_str, "rel_db": rel_db_str,
            "has_geometry": True,
            "predicted_label": m["classification"],
            "predicted_annotation": m["annotation"],
            "Tx": m["T"]["x"], "Ty": m["T"]["y"], "Tz": m["T"]["z"],
            "Rx": m["R"]["x"], "Ry": m["R"]["y"], "Rz": m["R"]["z"],
        })
    return rows


def attach_excel_labels(df, tables, freq_atol=0.05, freq_rtol=1e-4):
    """Join ref_label/ideal from Excel's data_score sheet onto `df`'s internal
    rows, molecule by molecule, gated by a whole-molecule frequency-agreement
    check (see module docstring). External rows are untouched (already
    correct, structural). Returns (df, skip_report) where skip_report is a
    list of {'molecule', 'n_mismatched', 'example': (mode_index, engine_freq,
    excel_freq)} dicts, one per molecule whose internal-row label join was
    skipped (excel_freq is None if no Excel row exists at all for that mode
    index -- treated identically to a numeric mismatch, not silently
    skipped, per the fail-loud guarantee established in the prior session).
    """
    ds = tables["data_score"].copy()
    ds = ds[~ds["molecule"].isin(EXCLUDED_MOLECULES)].copy()
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
            continue  # no Excel counterpart at all -- not an error (schema doc point 5)

        internal_idx = df.index[(df["molecule"] == mol) & (df["kind"] == "internal")]
        matched = []
        mismatches = []
        for idx in internal_idx:
            mode_index = int(df.at[idx, "mode_index"])
            engine_freq = df.at[idx, "freq"]
            ex_row = ds_by_key.get((mol, mode_index))
            if ex_row is None:
                mismatches.append((mode_index, engine_freq, None))
                continue
            excel_freq = ex_row["freq"]
            if pd.isna(excel_freq) or not np.isclose(engine_freq, excel_freq,
                                                       atol=freq_atol, rtol=freq_rtol):
                mismatches.append((mode_index, engine_freq,
                                    None if pd.isna(excel_freq) else float(excel_freq)))
                continue
            matched.append((idx, ex_row))

        if mismatches:
            skip_report.append({
                "molecule": mol, "n_mismatched": len(mismatches),
                "example": mismatches[0],
            })
            continue

        for idx, ex_row in matched:
            ref_label = ex_row["type"] if ex_row["type"] in ("bend", "stretch") else None
            df.at[idx, "ref_label"] = ref_label
            df.at[idx, "ideal"] = ex_row["ideal"]

    return df, skip_report


def build_library_scores(xlsx_path="data/vibrational-scoring-functions.xlsx",
                          data_dir="data", thresholds=None, return_skip_report=False):
    """Disk-driven library build: score EVERY molecule in
    discover_geometry_molecules() with the real engine, then join
    ref_label/ideal from Excel where a matching, frequency-consistent row
    exists. See module docstring for the full contract.
    """
    tables = load_excel_tables(xlsx_path)
    ds_molecules = set(tables["data_score"]["molecule"].dropna().unique()) - EXCLUDED_MOLECULES

    bases = discover_geometry_molecules(data_dir)
    all_rows = []
    load_errors = []
    for base in bases:
        molecule_label = resolve_excel_molecule_name(base, ds_molecules) or base
        try:
            rows = score_geometry_molecule(base, data_dir, thresholds)
        except (FileNotFoundError, ValueError) as e:
            load_errors.append((base, str(e)))
            continue
        for r in rows:
            r["molecule"] = molecule_label
        all_rows.extend(rows)

    df = pd.DataFrame(all_rows, columns=SCHEMA_COLUMNS)
    df, skip_report = attach_excel_labels(df, tables)

    for base, err in load_errors:
        warnings.warn(
            f"excel_ingest: skipped on-disk molecule '{base}' -- parse/scoring "
            f"failed ({err}) -- not included in library_scores.csv at all.",
            stacklevel=2)
    for entry in skip_report:
        mode_index, eng_f, exc_f = entry["example"]
        exc_f_str = f"{exc_f:.4f}" if exc_f is not None else "NO EXCEL ROW FOUND"
        warnings.warn(
            f"excel_ingest: '{entry['molecule']}' internal rows NOT label-"
            f"joined ({entry['n_mismatched']} mode(s) mismatched, e.g. mode "
            f"{mode_index}: engine {eng_f:.4f} vs Excel {exc_f_str} cm-1) -- "
            "this molecule's Excel data_score row likely came from a "
            "different calculation than data/logs/. Its scores (V_Stretch, "
            "Tx..Rz, predicted_label, ...) are still the real engine's own "
            "and are NOT affected; only ref_label/ideal are left null.",
            stacklevel=2)

    if return_skip_report:
        return df, skip_report
    return df


def run_ingest_pipeline(xlsx_path="data/vibrational-scoring-functions.xlsx",
                         data_dir="data", thresholds=None, write=True):
    """Headless entry point: build library_scores and optionally write the CSV."""
    df, skip_report = build_library_scores(xlsx_path, data_dir, thresholds,
                                            return_skip_report=True)
    out_path = os.path.join(data_dir, "results", "library_scores.csv")
    if write:
        df.to_csv(out_path, index=False)
    return df, out_path, skip_report
