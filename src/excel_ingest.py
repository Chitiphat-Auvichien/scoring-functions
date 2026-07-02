"""Phase-3 library ingest: data/vibrational-scoring-functions.xlsx -> data/results/library_scores.csv.

Per the locked decision in IMPLEMENTATION_PLAN.md ("ingest library scores rather
than recompute"), this module reads the hydride/AB_n library's PRECOMPUTED
scores directly out of the Excel workbook -- it never re-derives s[V_S] from
the Gaussian logs (Phase 0 already spot-verified the eq:vscore column identity
on H2S/SF2 to ~1e-5, well within 3 dp; see IMPLEMENTATION_PLAN.md Changelog).

Excel sheet layout (verified this session by direct inspection, not assumed
from sheet names -- the sibling ``Eckart``/``Eckart vs score`` sheets turned
out NOT to hold what their names implied, so nothing here is taken on faith):

``data_score`` -- one row per (molecule, mode), where "mode" is a 1-based
INTERNAL vibrational-mode index (3N-6 rows per molecule; there are NO
translation/rotation rows in this sheet at all -- confirmed for every
molecule, e.g. SnH4 has exactly 9 rows, C6H6 exactly 30). Columns used:
  - ``sum_|d2-d1|^2/sum(|d2-d1|^2)*cos(theta)``  -> s[V_S] (eq:vscore; Phase-0 verified)
  - ``ideal``   ('yes'/'no')     -> the ideal/non-ideal tag, straight from the sheet
                                    (no need to re-derive from tab:ideal/tab:nonideal)
  - ``type``    ('bend'/'stretch') -> the literature reference label. Cross-checked
                                    this session: byte-identical to the independent
                                    ``characterised modes`` sheet's own 'type' column
                                    for every one of the 600 rows that has an entry
                                    in both (the only rows without a ``characterised
                                    modes`` counterpart are the 87 'Gly5' rows, which
                                    are excluded anyway -- see EXCLUDED_MOLECULES).
  - ``freq``, ``Delta|b|`` (mean |relative bond-length change| over that mode's
    bonds -- verified this session against the independent per-bond average in
    ``data_mode&bond`` to 1e-16, i.e. it IS that average, not a different quantity).

``data_mode&bond`` -- one row per (molecule, mode, bond); supplies the per-bond
detail data_score does not carry:
  - ``|d2-d1|^2/sum(|d2-d1|^2)*|cos(theta)|`` -- per-bond s_AB (eq:bondscore).
    Verified this session: summed within each (molecule, mode) group, this
    reproduces data_score's s[V_S] column to ~5e-6 (600/600 rows checked) --
    it is not a different candidate formula, it is s_AB itself.
  - ``(|b'|-|b|)/|b|`` -- signed per-bond relative bond-length change (the
    per-bond counterpart of data_score's mode-averaged ``Delta|b|``).

T/R reference labels: data_score/data_mode&bond carry NO external-mode rows
at all, and constructing an ideal T/R reference requires real atomic geometry
that this repo simply does not have for most of the ~70 Excel molecules --
only 22 hydride-library molecules (+ water/benzene/CO2) ship with a matching
.log/.gjf pair in data/logs//data/gjf/ (see EXCEL_TO_LOG for the name-mapping
between the Excel molecule name and the on-disk log basename, plus a direct
case-match fallback for everything else; confirmed by cross-checking the
engine's parsed frequency against the Excel 'freq' column for 5 hydride
spot-check molecules this session: max abs deviation 0.0042 cm-1, consistent
with rounding, not misalignment -- CO2 and benzene likewise match to <=0.015
cm-1 across all their modes). For that geometry-backed subset,
``attach_geometry_classification`` runs the REAL classifier
(main.build_scorer_and_final + classify_all_modes) and (a) attaches predicted
labels/annotations/T,R scores onto the matching internal rows -- ONLY if every
one of that molecule's Excel-vs-engine frequencies agrees within tolerance
first (checked for the whole molecule before writing anything, so a bad
molecule cannot half-merge) -- and (b) always appends new rows (kind=
'external') for the n_T+n_R ideal T/R modes, which do not exist in the Excel
sheet at all and whose construction (Eckart-Sayvetz, from geometry alone) does
not depend on the vibrational-mode frequency match at all. **Discovered this
session:** water's Excel 'H2O' row does NOT match this repo's water.log --
Excel frequencies 1722.454/3501.541/3660.830 cm-1 vs. the engine's parsed
1628.029/3887.192/4005.506 cm-1 (differences of 94-386 cm-1, not rounding) --
so water's 3 internal rows are correctly left has_geometry=False (skipped, not
force-merged); water's 6 ideal T/R rows are still attached normally, since
those do not depend on this mismatch. This looks like the Excel 'H2O' row was
populated from a different water calculation (e.g. a different basis/method,
or literature/experimental frequencies) than the water.log used elsewhere in
this repo for tab:water -- flagged for the lead-author, not silently patched
over. CO2 and benzene both match cleanly (see above) and merge normally. For
the remaining ~45 Excel-only molecules (no .log/.gjf pair in this repo at
all), the output carries Step-4-relevant columns only (V_Stretch, freq,
Delta|b|, ideal tag, ref_label, per-bond detail); T/R/predicted-label columns
are left blank -- a genuine data-availability limit (no coordinates to build a
reference from), not a design shortcut, and is documented in the output
schema (``has_geometry`` column) so downstream consumers cannot miss it.

Excluded: 'Gly5' (87 rows, ``type=='0'`` in data_score -- a gramicidin
fragment; Decision 5, deferred to the companion paper). 'Gramicidin_A' never
appears in data_score at all (only in 'characterised modes', with no computed
scores), so it is naturally absent already.

Library-vs-manuscript-table coverage (checked exhaustively this session
against tab:ideal / tab:nonideal in JCC_temp_LaTeXtemplate.tex): ALL 11
tab:ideal entries (SnO2, TeH2, InH3, SbH3, IH3, SnH4, XeH4, TeH4, SbH5, XeOH4,
TeH6) and the full tab:nonideal roster -- bent AB2 (O/S/Se x H/F/Cl/Br) +
H2O2, trig-planar AB3 (B/Al/Ga x H/F/Cl/Br, minus BBr3/GaBr3) + C2H4,
trig-pyramidal AB3 (N/P/As x H/F/Cl/Br, complete) + tetrahedral AB4 (C/Si/Ge x
H/F/Cl/Br, minus CBr4/SiBr4/GeBr4) + C2H6 -- are present in ``data_score``.
The only real gap: 5 specific bromides (BBr3, GaBr3, CBr4, SiBr4, GeBr4) are
absent from the workbook entirely. Flagged for the lead-author/tex-data-sync;
not a large enough gap to force a manuscript table edit (5 of ~75 named
entries), but should be noted (e.g. a footnote) rather than silently ignored.
"""
import os
import warnings

import numpy as np
import openpyxl
import pandas as pd

# --- Column identities in the Excel workbook (verified 2026-07 session) ---
_VS_COL = "sum_|d₂-d₁|²/sum(|d₂-d₁|²)*cosθ"
_BOND_S_AB_COL = "|d₂-d₁|²/sum(|d₂-d₁|²)*|cosθ|"
_BOND_RELDB_COL = "(|b'|-|b|)/|b|"
_DELTA_B_COL = "Δ|b|"

EXCLUDED_MOLECULES = {"Gly5"}  # Decision 5: gramicidin fragment, out of scope this paper.

# Excel molecule name -> data/logs & data/gjf basename, for the subset that
# ships with real Gaussian geometry in this repo under a different filename
# (case or naming convention). Everything else is tried as a direct/case
# match against the molecule name itself (see resolve_log_basename).
EXCEL_TO_LOG = {
    "H2O": "water",
    "C6H6": "benzene",
    "CO2": "co2_mp2_3-21g",
    "Cl2O": "ocl2",
    "OF2": "of2",
    "Br2O": "br2o",
    "SeBr2": "SeBr2-cc",
}

_EXTERNAL_SLOTS = ("Tx", "Ty", "Tz", "Rx", "Ry", "Rz")


def _read_sheet(xlsx_path, sheet_name):
    wb = openpyxl.load_workbook(xlsx_path, data_only=True)
    ws = wb[sheet_name]
    rows = list(ws.iter_rows(values_only=True))
    header, body = rows[0], rows[1:]
    return pd.DataFrame(body, columns=header)


def load_excel_tables(xlsx_path):
    """Read the two sheets this module needs into DataFrames."""
    return {
        "data_score": _read_sheet(xlsx_path, "data_score"),
        "data_mode_bond": _read_sheet(xlsx_path, "data_mode&bond"),
    }


def resolve_log_basename(molecule, data_dir="data"):
    """Return the data/logs/<basename> stem for `molecule` if a matching .log
    (or .out) AND a connectivity .com/.gjf pair both exist in this repo;
    otherwise None (Excel-only molecule -- no geometry available here).
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


def _bond_string(rows, col):
    """Format a group of per-bond rows as 'atom1-atom2:value;atom1-atom2:value'."""
    if rows is None or len(rows) == 0:
        return ""
    parts = []
    for _, row in rows.iterrows():
        val = row[col]
        if pd.isna(val):
            continue
        parts.append(f"{row['atom1']}-{row['atom2']}:{float(val):.4f}")
    return ";".join(parts)


def ingest_internal_rows(tables):
    """Build one row per (molecule, mode) internal vibration from data_score
    + data_mode&bond. No geometry/classification columns yet (see
    attach_geometry_classification) -- this function only ever reads
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
    return pd.DataFrame(rows)


def attach_geometry_classification(df, data_dir="data", thresholds=None,
                                    freq_atol=0.05, freq_rtol=1e-4):
    """For every molecule in `df` that has a real .log/.gjf pair on disk, run
    the actual classifier (Algorithm 1) on its real geometry.

    Two independent attachments, handled separately because they have
    different validity conditions:
      (a) internal-row merge (predicted_label/annotation/T,R scores onto the
          matching "Vib i" rows) -- requires the WHOLE molecule's engine vs.
          Excel frequencies to agree within tolerance (checked before writing
          anything, so one bad molecule cannot half-merge); a molecule that
          fails this check is skipped (has_geometry stays False for its
          internal rows) and recorded in the returned skip report rather than
          raising -- this is a real, isolated data-provenance mismatch (see
          the water/H2O case in the module docstring), not necessarily a code
          bug, so it must not crash the whole ingest. A row whose "Vib i"
          name resolves to NO engine mode at all (m is None -- not currently
          known to happen for any molecule, but not proven impossible) is
          treated identically to a frequency mismatch, not silently skipped:
          it is appended to `mismatches` with engine_freq=None so it still
          trips the same all-or-nothing gate and shows up in the skip report
          (fixed 2026-07-02, formula-auditor finding -- the previous `continue`
          could in principle let a molecule half-merge with no report at all).
      (b) external-row append (the n_T+n_R ideal T/R rows) -- always done
          when the geometry parses at all, since construct_T()/construct_R()
          build these directly from geometry and do not depend on the
          vibrational-mode frequency match at all.

    Returns (df, skip_report) where skip_report is a list of
    {'molecule', 'n_mismatched', 'example': (mode_index, engine_freq, excel_freq)}
    dicts, one per molecule whose internal-row merge was skipped.
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
                # No engine mode resolves for this Excel row at all (not
                # currently known to happen for any molecule -- see the
                # module docstring -- but not silently skipped either: an
                # unresolved row is exactly as disqualifying as a frequency
                # mismatch, so it counts toward the mismatch gate below and
                # the whole molecule's internal-row merge is excluded and
                # reported, matching the documented all-or-nothing guarantee.
                mismatches.append((row["mode_index"], None, row["freq"]))
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
        df = pd.concat([df, pd.DataFrame(extra_rows)], ignore_index=True)
    return df, skip_report


def build_library_scores(xlsx_path="data/vibrational-scoring-functions.xlsx",
                          data_dir="data", thresholds=None, return_skip_report=False):
    tables = load_excel_tables(xlsx_path)
    df = ingest_internal_rows(tables)
    df, skip_report = attach_geometry_classification(df, data_dir, thresholds)
    for entry in skip_report:
        mode_index, eng_f, exc_f = entry["example"]
        eng_f_str = f"{eng_f:.4f}" if eng_f is not None else "NO ENGINE MODE RESOLVED"
        warnings.warn(
            f"excel_ingest: '{entry['molecule']}' internal rows NOT geometry-"
            f"merged ({entry['n_mismatched']} mode(s) mismatched, e.g. mode "
            f"{mode_index}: engine {eng_f_str} vs Excel {exc_f:.4f} cm-1) -- "
            "this molecule's Excel row likely came from a different "
            "calculation than data/logs/. Its ideal T/R rows are still "
            "attached (geometry-only, unaffected).", stacklevel=2)
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
