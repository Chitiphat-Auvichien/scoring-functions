"""Library ingest: builds data/results/library_scores.csv from TWO
switchable data sources, selected by the ``source`` parameter on
``build_library_scores()``/``run_ingest_pipeline()`` (``"excel"`` or
``"gaussian"``). Both produce the exact same locked ``SCHEMA_COLUMNS`` output
schema, so ``src/calibrate.py``/``src/figures.py`` never need to know or care
which source produced a given ``library_scores.csv``.

**Why two sources (2026-07-04 author decision -- refinement/interim mode of
the 2026-07-03 "recompute from Gaussian" directive below, NOT a reversal of
it; see IMPLEMENTATION_PLAN.md's RESUME HERE for the full record).** The
manuscript's already-typeset figures (fig:confusion, fig:bondscores,
fig:boxplots, fig:modemixing) were built from the OLD Excel-scored ~70-
molecule ``library_scores.csv``, before the 2026-07-03 disk-driven
rearchitecture shrank the population to ~25 molecules (only those with
on-disk Gaussian ``.log``+``.gjf`` pairs). The author is adding real Gaussian
files for the rest of the hydride library incrementally but is having
Gaussian trouble right now, so the interim default is ``source="excel"``
(reproducing the pre-rearchitecture, manuscript-matching population and
values) while the disk-driven ``source="gaussian"`` path is kept fully
working, ready to become the default again once enough ``.log``/``.gjf``
pairs exist on disk.

``source="excel"`` (default; commit-149fc62 logic, restored verbatim):
  Reads the hydride-library's PRECOMPUTED scores directly out of
  ``data/vibrational-scoring-functions.xlsx`` -- never re-derives s[V_S] from
  the Gaussian logs (Phase 0 already spot-verified the eq:vscore column
  identity on H2S/SF2 to ~1e-5, well within 3 dp; see IMPLEMENTATION_PLAN.md
  Changelog). ``data_score`` supplies one row per (molecule, 1-based internal
  vibrational-mode index) directly: ``V_Stretch``, ``delta_b_mean``,
  ``ref_label`` ('type'), ``ideal``, ``freq``. ``data_mode&bond`` supplies the
  per-bond detail (``s_AB``, ``rel_db``) that ``data_score`` does not carry.
  For the subset of molecules that ALSO have a real ``.log``/``.gjf`` pair on
  disk, ``attach_geometry_classification`` runs the real engine
  (``main.build_scorer_and_final`` + ``classify_all_modes``) and OVERLAYS
  ``predicted_label``/``predicted_annotation``/``Tx..Rz``/``has_geometry`` onto
  the matching internal rows (gated by a whole-molecule engine-vs-Excel
  frequency-agreement check -- one mismatched mode skips the WHOLE molecule's
  overlay, so a bad molecule cannot half-merge) and ALWAYS appends the n_T+n_R
  ideal external T/R rows (construct_T()/construct_R() need only geometry, not
  a frequency match). ``V_Stretch``/``delta_b_mean``/``s_AB``/``rel_db``/
  ``ref_label``/``ideal`` on internal rows are ALWAYS Excel's own values in
  this mode, even for geometry-backed molecules -- only the engine-derived
  overlay columns (``predicted_label``, ``predicted_annotation``, ``Tx..Rz``,
  ``has_geometry``) and the always-geometry-only external rows come from the
  engine. Water's Excel 'H2O' row does NOT match this repo's ``water.log``
  (Excel freqs 1722.454/3501.541/3660.830 cm-1 vs. the engine's parsed
  1628.029/3887.192/4005.506 cm-1) -- so water's internal rows correctly stay
  ``has_geometry=False`` (Excel-only scores, not force-merged), while its 6
  ideal T/R rows are attached normally. Same story for OF2/Cl2O/Br2O.
  Excel-only molecules (no on-disk ``.log``/``.gjf`` pair at all -- ~45 of the
  ~70) get Step-4-relevant columns only (``V_Stretch``, ``freq``,
  ``delta_b_mean``, ``ideal``, ``ref_label``, per-bond detail);
  ``predicted_label``/``Tx..Rz`` are left blank (``has_geometry=False``) -- a
  genuine data-availability limit, not a design shortcut.

``source="gaussian"`` (2026-07-03 disk-driven rearchitecture, unchanged):
  ``discover_geometry_molecules()`` lists every basename present in BOTH
  ``data/logs/`` (``.log``/``.out``) AND ``data/gjf/`` (``.com``/``.gjf``);
  ``score_geometry_molecule()`` runs the real engine
  (``main.load_inputs`` -> ``build_scorer_and_final`` ->
  ``src.classifier.classify_all_modes``, ``ModeScorer.score_bonds()``) on
  each one to build EVERY score column (``V_Stretch``/``Tx..Rz``/
  ``predicted_label``/``s_AB``/``rel_db``/``delta_b_mean``) from scratch, for
  every molecule physically present on disk -- Excel supplies ONLY
  ``ref_label``/``ideal`` (``attach_excel_labels()``, joined by
  (molecule, mode index), gated by the same whole-molecule frequency check).
  Molecules with no on-disk geometry are absent from the output entirely.
  This scales with zero code changes as more ``.log``/``.gjf`` pairs are
  dropped in -- the loop is keyed off ``os.listdir(data/logs)`` ∩
  ``os.listdir(data/gjf)``, not off Excel rows. This is the intended
  EVENTUAL default once enough Gaussian files exist for the full library;
  it is being kept alive and ready, just not the default output right now.

Both sources ultimately write the exact same locked schema (unchanged column
names so downstream consumers keep working):
  ``molecule``, ``mode_index``, ``kind``, ``freq``, ``ref_label``, ``ideal``,
  ``V_Stretch``, ``delta_b_mean``, ``s_AB``, ``rel_db``, ``has_geometry``,
  ``predicted_label``, ``predicted_annotation``, ``Tx``, ``Ty``, ``Tz``,
  ``Rx``, ``Ry``, ``Rz``.

Excluded (both sources): 'Gly5' (a gramicidin fragment; Decision 5, deferred
to the companion paper).

``EXCEL_TO_LOG``/``resolve_log_basename``/``resolve_excel_molecule_name`` are
shared by both sources (name-mapping between the Excel molecule name and the
on-disk ``data/logs``+``data/gjf`` basename, for the subset that ships under a
different filename/casing convention).
"""
import os
import warnings

import numpy as np
import openpyxl
import pandas as pd

# --- Column identities in the Excel workbook (verified 2026-07 session;
# used only by source="excel") ---
_VS_COL = "sum_|d₂-d₁|²/sum(|d₂-d₁|²)*cosθ"
_BOND_S_AB_COL = "|d₂-d₁|²/sum(|d₂-d₁|²)*|cosθ|"
_BOND_RELDB_COL = "(|b'|-|b|)/|b|"
_DELTA_B_COL = "Δ|b|"

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
# column names) -- identical regardless of source.
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


def load_excel_tables(xlsx_path, include_bonds=False):
    """Read the Excel sheet(s) needed by the current source. ``data_score``
    (label lookup, both sources) is always read; ``data_mode&bond`` (per-bond
    precomputed scores, source="excel" only) is read only if
    ``include_bonds=True`` to avoid the extra I/O cost for source="gaussian",
    which never needs it."""
    tables = {"data_score": _read_sheet(xlsx_path, "data_score")}
    if include_bonds:
        tables["data_mode_bond"] = _read_sheet(xlsx_path, "data_mode&bond")
    return tables


def resolve_log_basename(molecule, data_dir="data"):
    """Return the data/logs/<basename> stem for `molecule` if a matching .log
    (or .out) AND a connectivity .com/.gjf pair both exist in this repo;
    otherwise None (Excel-only molecule -- no geometry available here).
    Shared by both sources.
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
    data/gjf/ (.com or .gjf) -- the disk-driven source of truth for
    source="gaussian"'s molecule loop. Scales automatically as more files are
    dropped in; no code changes needed."""
    logs_dir = os.path.join(data_dir, "logs")
    gjf_dir = os.path.join(data_dir, "gjf")
    log_bases = {os.path.splitext(f)[0] for f in os.listdir(logs_dir)
                 if f.lower().endswith((".log", ".out"))}
    gjf_bases = {os.path.splitext(f)[0] for f in os.listdir(gjf_dir)
                 if f.lower().endswith((".com", ".gjf"))}
    return sorted(log_bases & gjf_bases)


def resolve_excel_molecule_name(base, excel_molecules):
    """On-disk basename -> Excel data_score molecule name, or None if `base`
    has no Excel counterpart at all (not an error -- used by source="gaussian")."""
    if base in _LOG_TO_EXCEL:
        return _LOG_TO_EXCEL[base]
    if base in excel_molecules:
        return base
    return None


def _bond_string(rows, col):
    """Format a group of per-bond rows as 'atom1-atom2:value;atom1-atom2:value'
    (source="excel" helper, reading the precomputed per-bond Excel columns)."""
    if rows is None or len(rows) == 0:
        return ""
    parts = []
    for _, row in rows.iterrows():
        val = row[col]
        if pd.isna(val):
            continue
        parts.append(f"{row['atom1']}-{row['atom2']}:{float(val):.4f}")
    return ";".join(parts)


# ---------------------------------------------------------------------------
# source="gaussian": disk-driven, real-engine recompute (2026-07-03).
# ---------------------------------------------------------------------------

def score_geometry_molecule(base, data_dir="data", thresholds=None):
    """Run the real engine (Steps 1-4) on one on-disk molecule and return a
    list of row dicts (one per external T/R slot + one per internal 'Vib i'
    mode) in the library_scores.csv schema, EXCLUDING 'molecule' (the caller
    attaches that) and with ref_label/ideal left None (attached separately by
    attach_excel_labels, label-only, per the source="gaussian" architecture).

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

        # i_label/j_label are ModeScorer.score_bonds()'s own atom-symbol +
        # 1-based-index labels (e.g. "C1", "H7"), matching the convention
        # src/benzene_validation.py's C-C/C-H bond-type parsing relies on.
        s_ab_str = ";".join(f"{b['i_label']}-{b['j_label']}:{b['s_AB']:.4f}" for b in bonds)
        rel_db_str = ";".join(f"{b['i_label']}-{b['j_label']}:{b['rel_db']:.4f}" for b in bonds)
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
    """(source="gaussian") Join ref_label/ideal from Excel's data_score sheet
    onto `df`'s internal rows, molecule by molecule, gated by a whole-molecule
    frequency-agreement check (see module docstring). External rows are
    untouched (already correct, structural). Returns (df, skip_report) where
    skip_report is a list of {'molecule', 'n_mismatched', 'example':
    (mode_index, engine_freq, excel_freq)} dicts, one per molecule whose
    internal-row label join was skipped (excel_freq is None if no Excel row
    exists at all for that mode index -- treated identically to a numeric
    mismatch, not silently skipped, per the fail-loud guarantee established in
    a prior session).
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
            continue  # no Excel counterpart at all -- not an error

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


def _build_library_scores_gaussian(xlsx_path, data_dir, thresholds, return_skip_report):
    """source="gaussian": disk-driven, real-engine recompute for every
    molecule with a .log+.gjf pair on disk; Excel supplies ref_label/ideal
    only. See module docstring."""
    tables = load_excel_tables(xlsx_path, include_bonds=False)
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
            f"excel_ingest[gaussian]: skipped on-disk molecule '{base}' -- parse/scoring "
            f"failed ({err}) -- not included in library_scores.csv at all.",
            stacklevel=3)
    for entry in skip_report:
        mode_index, eng_f, exc_f = entry["example"]
        exc_f_str = f"{exc_f:.4f}" if exc_f is not None else "NO EXCEL ROW FOUND"
        warnings.warn(
            f"excel_ingest[gaussian]: '{entry['molecule']}' internal rows NOT label-"
            f"joined ({entry['n_mismatched']} mode(s) mismatched, e.g. mode "
            f"{mode_index}: engine {eng_f:.4f} vs Excel {exc_f_str} cm-1) -- "
            "this molecule's Excel data_score row likely came from a "
            "different calculation than data/logs/. Its scores (V_Stretch, "
            "Tx..Rz, predicted_label, ...) are still the real engine's own "
            "and are NOT affected; only ref_label/ideal are left null.",
            stacklevel=3)

    if return_skip_report:
        return df, skip_report
    return df


# ---------------------------------------------------------------------------
# source="excel": precomputed Excel scores + geometry overlay where available
# (commit-149fc62 logic, restored 2026-07-04).
# ---------------------------------------------------------------------------

def ingest_internal_rows(tables):
    """(source="excel") Build one row per (molecule, mode) internal vibration
    from data_score + data_mode&bond. No geometry/classification columns yet
    (see attach_geometry_classification) -- this function only ever reads
    precomputed Excel values, never recomputes a score.
    """
    ds = tables["data_score"].copy()
    dmb = tables["data_mode_bond"].copy()
    ds = ds[~ds["molecule"].isin(EXCLUDED_MOLECULES)].copy()
    dmb = dmb[~dmb["molecule"].isin(EXCLUDED_MOLECULES)].copy()

    for col in ("freq", _VS_COL, _DELTA_B_COL):
        ds[col] = pd.to_numeric(ds[col], errors="coerce")
    for col in (_BOND_S_AB_COL, _BOND_RELDB_COL):
        dmb[col] = pd.to_numeric(dmb[col], errors="coerce")

    bond_groups = {key: grp for key, grp in dmb.groupby(["molecule", "mode"])}

    rows = []
    for _, r in ds.iterrows():
        key = (r["molecule"], r["mode"])
        bonds = bond_groups.get(key)
        ref_label = r["type"] if pd.notna(r["type"]) and r["type"] in ("bend", "stretch") else None
        rows.append({
            "molecule": r["molecule"],
            "mode_index": int(r["mode"]),
            "kind": "internal",
            "freq": r["freq"],
            "ref_label": ref_label,
            "ideal": r["ideal"],
            "V_Stretch": r[_VS_COL],
            "delta_b_mean": r[_DELTA_B_COL],
            "s_AB": _bond_string(bonds, _BOND_S_AB_COL),
            "rel_db": _bond_string(bonds, _BOND_RELDB_COL),
            "has_geometry": False,
            "predicted_label": None,
            "predicted_annotation": None,
            "Tx": None, "Ty": None, "Tz": None,
            "Rx": None, "Ry": None, "Rz": None,
        })
    return pd.DataFrame(rows, columns=SCHEMA_COLUMNS)


def attach_geometry_classification(df, data_dir="data", thresholds=None,
                                    freq_atol=0.05, freq_rtol=1e-4):
    """(source="excel") For every molecule in `df` that has a real .log/.gjf
    pair on disk, run the actual classifier (Algorithm 1) on its real
    geometry.

    Two independent attachments, handled separately because they have
    different validity conditions:
      (a) internal-row overlay (predicted_label/annotation/T,R scores onto the
          matching "Vib i" rows) -- requires the WHOLE molecule's engine vs.
          Excel frequencies to agree within tolerance (checked before writing
          anything, so one bad molecule cannot half-merge); a molecule that
          fails this check is skipped (has_geometry stays False for its
          internal rows) and recorded in the returned skip report rather than
          raising -- this is a real, isolated data-provenance mismatch (see
          the water/H2O case in the module docstring), not necessarily a code
          bug, so it must not crash the whole ingest. V_Stretch/delta_b_mean/
          s_AB/rel_db/ref_label/ideal are NEVER touched here -- those stay
          Excel's own values from ingest_internal_rows regardless.
      (b) external-row append (the n_T+n_R ideal T/R rows) -- always done
          when the geometry parses at all, since construct_T()/construct_R()
          build these directly from geometry and do not depend on the
          vibrational-mode frequency match at all.

    Returns (df, skip_report) where skip_report is a list of
    {'molecule', 'n_mismatched', 'example': (mode_index, engine_freq, excel_freq)}
    dicts, one per molecule whose internal-row overlay was skipped.
    """
    from main import load_inputs, build_scorer_and_final
    from src.classifier import classify_all_modes

    df = df.copy()
    extra_rows = []
    skip_report = []
    for mol in df["molecule"].unique():
        base = resolve_log_basename(mol, data_dir)
        if base is None:
            continue
        try:
            raw, _ = load_inputs(base, "normal", data_dir)
            scorer, final = build_scorer_and_final(raw, "normal")
        except (FileNotFoundError, ValueError):
            # Missing bonds / parse trouble -- skip rather than fabricate geometry.
            continue
        scored = classify_all_modes(scorer, final, thresholds)
        by_name = {m["name"]: m for m in scored}

        mol_mask = df["molecule"] == mol
        mol_rows = df[mol_mask]

        matched = []          # (idx, engine_mode_dict) pairs that align
        mismatches = []        # (mode_index, engine_freq, excel_freq)
        for idx, row in mol_rows.iterrows():
            m = by_name.get(f"Vib {row['mode_index']}")
            if m is None:
                continue
            if np.isclose(m["frequency"], row["freq"], atol=freq_atol, rtol=freq_rtol):
                matched.append((idx, m))
            else:
                mismatches.append((row["mode_index"], m["frequency"], row["freq"]))

        if mismatches:
            skip_report.append({
                "molecule": mol, "n_mismatched": len(mismatches),
                "example": mismatches[0],
            })
        else:
            for idx, m in matched:
                df.at[idx, "has_geometry"] = True
                df.at[idx, "predicted_label"] = m["classification"]
                df.at[idx, "predicted_annotation"] = m["annotation"]
                for ax in "xyz":
                    df.at[idx, f"T{ax}"] = m["T"][ax]
                    df.at[idx, f"R{ax}"] = m["R"][ax]

        # External (T/R) rows are independent of the internal-row frequency
        # match above -- attach them whenever geometry parsed successfully.
        for m in scored:
            if m["name"] in _EXTERNAL_SLOTS:
                extra_rows.append({
                    "molecule": mol, "mode_index": m["name"], "kind": "external",
                    "freq": m["frequency"],
                    "ref_label": "translation" if m["name"][0] == "T" else "rotation",
                    "ideal": None, "V_Stretch": m["V"], "delta_b_mean": None,
                    "s_AB": "", "rel_db": "", "has_geometry": True,
                    "predicted_label": m["classification"],
                    "predicted_annotation": m["annotation"],
                    "Tx": m["T"]["x"], "Ty": m["T"]["y"], "Tz": m["T"]["z"],
                    "Rx": m["R"]["x"], "Ry": m["R"]["y"], "Rz": m["R"]["z"],
                })
    if extra_rows:
        df = pd.concat([df, pd.DataFrame(extra_rows, columns=SCHEMA_COLUMNS)], ignore_index=True)
    return df, skip_report


def _build_library_scores_excel(xlsx_path, data_dir, thresholds, return_skip_report):
    """source="excel": precomputed Excel scores for the full ~70-molecule
    library, with a real-engine geometry overlay (predicted_label/Tx..Rz +
    appended external T/R rows) for the subset that also has on-disk
    geometry. See module docstring."""
    tables = load_excel_tables(xlsx_path, include_bonds=True)
    df = ingest_internal_rows(tables)
    df, skip_report = attach_geometry_classification(df, data_dir, thresholds)
    df = df.reindex(columns=SCHEMA_COLUMNS)
    for entry in skip_report:
        mode_index, eng_f, exc_f = entry["example"]
        warnings.warn(
            f"excel_ingest[excel]: '{entry['molecule']}' internal rows NOT geometry-"
            f"overlaid ({entry['n_mismatched']} mode(s) mismatched, e.g. mode "
            f"{mode_index}: engine {eng_f:.4f} vs Excel {exc_f:.4f} cm-1) -- "
            "this molecule's Excel row likely came from a different "
            "calculation than data/logs/. Its ideal T/R rows are still "
            "attached (geometry-only, unaffected); its internal-row scores "
            "remain Excel's own (V_Stretch/s_AB/rel_db/ref_label/ideal are "
            "never touched by this overlay).", stacklevel=3)
    if return_skip_report:
        return df, skip_report
    return df


# ---------------------------------------------------------------------------
# Public, source-switchable entry points.
# ---------------------------------------------------------------------------

def build_library_scores(xlsx_path="data/vibrational-scoring-functions.xlsx",
                          data_dir="data", thresholds=None, return_skip_report=False,
                          source="excel"):
    """Build the library_scores DataFrame from either data source.

    source="excel" (default, interim per the 2026-07-04 author decision --
        see module docstring): precomputed Excel scores for the full
        hydride library (~70 molecules), geometry-overlaid where an on-disk
        .log/.gjf pair also exists. Matches the population/values the
        manuscript's currently-typeset figures were built from.
    source="gaussian": disk-driven real-engine recompute for every molecule
        with a .log+.gjf pair on disk (~25 today, growing as more are added).
        The intended eventual default once the full library has on-disk
        geometry.

    Either way the returned DataFrame has exactly SCHEMA_COLUMNS.
    """
    if source == "excel":
        return _build_library_scores_excel(xlsx_path, data_dir, thresholds, return_skip_report)
    elif source == "gaussian":
        return _build_library_scores_gaussian(xlsx_path, data_dir, thresholds, return_skip_report)
    else:
        raise ValueError(f"build_library_scores: source must be 'excel' or 'gaussian', got {source!r}")


def run_ingest_pipeline(xlsx_path="data/vibrational-scoring-functions.xlsx",
                         data_dir="data", thresholds=None, write=True, source="excel"):
    """Headless entry point: build library_scores (from `source`) and
    optionally write the CSV. See build_library_scores()'s docstring for the
    source="excel"/"gaussian" contract."""
    df, skip_report = build_library_scores(xlsx_path, data_dir, thresholds,
                                            return_skip_report=True, source=source)
    out_path = os.path.join(data_dir, "results", "library_scores.csv")
    if write:
        df.to_csv(out_path, index=False)
    return df, out_path, skip_report
