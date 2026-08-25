import argparse
import os
import numpy as np
import pandas as pd

from src.parser import GaussianParser, EMITParser, IntermediateIO
from src.scoring import ModeScorer, V_WEIGHTINGS, DEFAULT_V_WEIGHTING, set_v_weighting
from src.classifier import classify_all_modes, classify_to_rows, is_linear, Thresholds
from src.projection import (build_reference_basis, project_emit,
                             build_reference_basis_cartesian, project_emit_cartesian)
from src.utils import find_file


def resolve_dirs(data_dir="data"):
    """Return the standard data subdirectories, creating them if needed."""
    dirs = {k: os.path.join(data_dir, k)
            for k in ("logs", "EMIT", "gjf", "intermediate", "results")}
    for d in dirs.values():
        os.makedirs(d, exist_ok=True)
    return dirs


def _find_log(logs_dir, mol_name):
    p = find_file(logs_dir, mol_name, (".log", ".out"))
    if p is None:
        raise FileNotFoundError(f"No log file for '{mol_name}' in {logs_dir} (expected .log or .out)")
    return p


def _find_emit(emit_dir, mol_name):
    for name in (f"{mol_name}_EMIT.txt", f"{mol_name}.txt"):
        p = os.path.join(emit_dir, name)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        f"No EMIT file for '{mol_name}' in {emit_dir} (expected {mol_name}_EMIT.txt)")


def _find_gjf(gjf_dir, mol_name):
    return find_file(gjf_dir, mol_name, (".com", ".gjf"))


def intermediate_path(dirs, mol_name, mode_type):
    """data/intermediate/<mol>_<normal|emit>_data.txt -- mode-type-suffixed so
    a molecule with both a normal-mode and an EMIT intermediate doesn't have
    the two overwrite each other (the older, unsuffixed '<mol>_data.txt'
    naming had exactly that collision)."""
    return os.path.join(dirs["intermediate"], f"{mol_name}_{mode_type}_data.txt")


def _cache_is_fresh(inter_path, source_paths):
    """True iff `inter_path` exists and is at least as new as every file in
    `source_paths` (mtime-based invalidation), so a stale cache (source
    regenerated/edited) is regenerated. Manually-added bonds persist across
    runs because saving them bumps the intermediate's mtime past its sources."""
    if not os.path.exists(inter_path):
        return False
    inter_mtime = os.path.getmtime(inter_path)
    return all(inter_mtime >= os.path.getmtime(p) for p in source_paths)


def _normal_cache_is_current_format(inter_path):
    """True iff a 'normal' intermediate is in the new Gaussian-direct format,
    not the old 'MOLECULE_DATA' format it used to share with EMIT
    intermediates. A pre-rework cache is mtime-fresh but structurally lacks
    reduced_mass/force_constant/irrep, and mtime alone can't detect that --
    this forces one self-healing reparse-and-rewrite in the new format."""
    with open(inter_path, 'r') as f:
        first_line = f.readline().strip()
    return first_line != "MOLECULE_DATA"


def load_inputs(mol_name, mode_type, data_dir="data", use_cache=True):
    """Parse geometry, modes, and connectivity for a molecule, via the
    mtime-invalidated intermediate cache (see _cache_is_fresh). mode_type is
    'normal' or 'emit'. Returns (raw_data, dirs); raises on missing files or
    wrong mode counts. Does NOT prompt; suitable for headless/reproduce use.
    Pass `use_cache=False` to force a fresh parse (still refreshes the cache).

    If `raw['bonds']` comes back empty, main()'s interactive fallback lets the
    user hand-edit bonds into the written intermediate; that edit's mtime
    keeps the cache fresh afterward so the connectivity persists.
    """
    if mode_type not in ("normal", "emit"):
        raise ValueError(f"mode_type must be 'normal' or 'emit', got {mode_type!r}")
    dirs = resolve_dirs(data_dir)
    inter_path = intermediate_path(dirs, mol_name, mode_type)

    log_path = _find_log(dirs["logs"], mol_name)
    source_paths = [log_path]
    gjf_path = _find_gjf(dirs["gjf"], mol_name)
    if gjf_path:
        source_paths.append(gjf_path)
    if mode_type == "emit":
        emit_path = _find_emit(dirs["EMIT"], mol_name)
        source_paths.append(emit_path)

    cache_fresh = use_cache and _cache_is_fresh(inter_path, source_paths)
    if cache_fresh and mode_type == "normal" and os.path.exists(inter_path):
        cache_fresh = _normal_cache_is_current_format(inter_path)
    if cache_fresh:
        return IntermediateIO.load(inter_path), dirs

    gp = GaussianParser(log_path)
    if mode_type == "normal":
        raw = gp.parse(parse_modes=True)
    else:
        raw = gp.parse(parse_modes=False)
        raw["modes"] = EMITParser(emit_path, len(raw["atoms"])).parse()

    IntermediateIO.save(raw, inter_path)
    return raw, dirs


def build_scorer_and_final(raw, mode_type, v_weighting=None):
    """Align to principal axes and build the candidate mode pool ('final'),
    shared by run_scoring_pipeline() and the classifier so both score the
    identical mode list.
      - 'normal': MIT rotates the molecule AND the mode vectors
        (rotate_modes=True); the 3 ideal T + 3 ideal R references are
        prepended to the real vibrational modes.
      - 'emit'  : EMIT modes already live in principal axes, so MIT rotates
        only the molecule (rotate_modes=False); no ideal references are
        added -- all 3N raw EMIT eigenvectors are the candidate pool.

    v_weighting selects the eq:vscore bond weighting ('mu' or 'none'); None
    means the process-active default set by --v-weighting. This function is the
    only place the app constructs a ModeScorer, so every caller that routes
    through it (library_ingest, calibrate, merge_ped_scores, compare_rerun)
    inherits the variant automatically.

    Raises ValueError if no bond connectivity is available (fail loud rather
    than silently corrupt the V-score). Returns (scorer, final).

    Linear-molecule invariant (authoritative site): a linear molecule has
    n_R=2, but MIT() always places the linear (smallest-moment) axis on the
    new X axis, so construct_R()'s "Rx" reference is an all-zero placeholder,
    not a genuine external mode. It is dropped here -- left in, it would fall
    through Step 2 (external_slots() excludes Rx for linear molecules, so no
    slot claims it) into Step 4 and be mislabeled BENDING (V=0 <= tau_B).
    """
    if not raw["bonds"]:
        raise ValueError(
            "No bond connectivity available; provide a .com/.gjf file in data/gjf/ "
            "or add bonds to the intermediate file before scoring.")
    rotate = (mode_type == "normal")
    scorer = ModeScorer(raw["atoms"], raw["coords"], raw["bonds"],
                        v_weighting=v_weighting)
    rotated = scorer.MIT(raw["modes"], rotate_modes=rotate)
    if mode_type == "normal":
        for i, m in enumerate(rotated):
            m["label"] = f"Vib {i+1}"
        ideal_R = scorer.construct_R()
        if is_linear(scorer):
            ideal_R = [m for m in ideal_R if m["label"] != "Rx"]
        final = scorer.construct_T() + ideal_R + rotated
    else:
        final = rotated
    return scorer, final


def run_scoring_pipeline(mol_name, mode_type, data_dir="data", thresholds=None, write=True,
                          scheme="binary", identify_tr=True):
    """Headless pipeline: load inputs -> classify_all_modes (Algorithm 1) -> CSV.

    This is the single, complete per-molecule result: scores (Tx..Rz,
    V_Stretch), Mu/K/Irrep, and the Step 2/3 classification (vib_label,
    tr_label, tr_score, s_AB) all in one row per mode -- there is no separate
    scores-only output. Writes data/results/<mol>_{normal,EMIT}.csv.

    `scheme` ("binary" default, paper-standard, or "threeway" explicit
    opt-in) is passed straight through to classify_all_modes() -- see that
    function's docstring. Default tracks classify_all_modes()'s own default
    so callers that don't pass `scheme` (e.g. reproduce.py's per-molecule
    stages) automatically follow the global scheme choice. `identify_tr`
    (default True, matching classify_all_modes()'s own default) toggles the
    optional Step 3 T/R identification; every reproduce.py/headless default
    run keeps it on.
    """
    raw, dirs = load_inputs(mol_name, mode_type, data_dir)
    scorer, final = build_scorer_and_final(raw, mode_type)
    scored = classify_all_modes(scorer, final, thresholds, scheme=scheme, identify_tr=identify_tr)
    df = pd.DataFrame(classify_to_rows(scored))
    suffix = "normal" if mode_type == "normal" else "EMIT"
    output_file = os.path.join(dirs["results"], f"{mol_name}_{suffix}.csv")
    if write:
        df.to_csv(output_file, index=False, float_format="%.4f")
    return df, output_file


def run_projection_pipeline(mol_name, data_dir="data", thresholds=None, write=True):
    """Headless EMIT -> normal-mode projection pipeline (Phase 2, eq:emitproj).

    Builds the mass-weighted normal-mode reference basis from data/logs/<mol>.log,
    projects the raw EMIT eigenvectors from data/EMIT/<mol>_EMIT.txt onto it, and
    merges the grouped T/R/V fractions (C2_Tx..C2_VMix) into
    data/results/<mol>_EMIT.csv in place (must already exist -- run
    `python main.py -m <mol> --mode emit` first), plus a separate
    <mol>_EMIT_full.csv with the per-reference-mode detail.

    Also runs the explicitly non-orthonormal plain-Cartesian-overlap
    comparison pathway (src.projection.build_reference_basis_cartesian /
    project_emit_cartesian -- see that module's docstring for why it exists
    only as a point of comparison, and why it is NOT squared like the C2_*
    columns above): merges "Ocart_Tx".."Ocart_Sum" (raw signed overlaps) into
    the same <mol>_EMIT.csv, and writes the per-reference-mode detail to a
    separate <mol>_EMIT_full_cartesian.csv.

    Returns (df_emit, df_full, (path_emit, path_full, path_full_cartesian)).
    """
    raw_n, dirs = load_inputs(mol_name, "normal", data_dir)
    scorer_n, final_n = build_scorer_and_final(raw_n, "normal")
    raw_e, _ = load_inputs(mol_name, "emit", data_dir)
    scorer_e, final_e = build_scorer_and_final(raw_e, "emit")

    # Q and Theta are only frame-consistent because both parses read the SAME
    # log geometry and MIT()'s eigendecomposition is deterministic -- assert
    # it rather than relying on that silently.
    coords_n = scorer_n.coords
    coords_e = scorer_e.coords
    if not np.allclose(coords_n, coords_e, atol=1e-6):
        raise ValueError(
            f"'{mol_name}': the 'normal' and 'emit' parses of the log geometry "
            "rotated into different principal-axis frames (max coord diff "
            f"{np.max(np.abs(coords_n - coords_e)):.3g}) -- Q and Theta would "
            "not be projection-comparable. Check for a second/different "
            "'Standard orientation' block or non-deterministic MIT() sign fix."
        )

    emit_path = os.path.join(dirs["results"], f"{mol_name}_EMIT.csv")
    if not os.path.exists(emit_path):
        raise FileNotFoundError(
            f"{emit_path} not found -- run `python main.py -m {mol_name} --mode emit` first.")
    df_emit = pd.read_csv(emit_path)
    # Re-running --emit-projection (e.g. after a fresh --mode emit) must not
    # duplicate C2_*/Ocart_* columns via merge's _x/_y suffixing -- drop any
    # already merged in from a prior run first, so this is idempotent.
    # ("C2cart_" is dropped too: an earlier revision of this pathway squared
    # the overlap under that prefix before Ocart_ replaced it -- see
    # src/projection.py's project_emit_cartesian docstring.)
    df_emit = df_emit.drop(columns=[c for c in df_emit.columns
                                     if c.startswith("C2_") or c.startswith("C2cart_")
                                     or c.startswith("Ocart_")])

    ref = build_reference_basis(scorer_n, final_n, thresholds)
    rows, full_rows = project_emit(ref, final_e)

    contrib = pd.DataFrame(rows).drop(columns=["Eigenvalue"])
    df_emit = df_emit.merge(contrib, on="Mode", how="left")
    df_full = pd.DataFrame(full_rows)
    out_full = os.path.join(dirs["results"], f"{mol_name}_EMIT_full.csv")

    ref_cart = build_reference_basis_cartesian(scorer_n, final_n, thresholds)
    rows_cart, full_rows_cart = project_emit_cartesian(ref_cart, final_e)

    contrib_cart = pd.DataFrame(rows_cart).drop(columns=["Eigenvalue"])
    df_emit = df_emit.merge(contrib_cart, on="Mode", how="left")
    df_full_cart = pd.DataFrame(full_rows_cart)
    out_full_cart = os.path.join(dirs["results"], f"{mol_name}_EMIT_full_cartesian.csv")

    if write:
        df_emit.to_csv(emit_path, index=False, float_format="%.4f")
        df_full.to_csv(out_full, index=False, float_format="%.6f")
        df_full_cart.to_csv(out_full_cart, index=False, float_format="%.6f")
    return df_emit, df_full, (emit_path, out_full, out_full_cart)


def _print_table(df, title):
    """Echo a results DataFrame via pandas to_string() -- borderless and paste-friendly into Excel."""
    print(f"\n--- {title} ---")
    print(df.to_string(index=False, float_format="%.4f"))


def _run_flag_pipelines(args, thresholds=None):
    """Handle the developer/maintainer CLI flags (--emit-projection,
    --ped-merge, --ped-merge-all, --library, --calibrate, --figures) by
    wiring up the existing headless pipeline functions. The plain user
    workflow (-m/--molecule [--mode]) is NOT handled here -- see main().

      - --library/--calibrate/--figures/--ped-merge-all are GLOBAL and ignore
        -m/--molecule.
      - --emit-projection/--ped-merge are PER-MOLECULE and require -m.
        --emit-projection requires <mol>_EMIT.csv to already exist (run
        `-m <mol> --mode emit` first); --ped-merge requires <mol>_normal.csv
        to already exist (run `-m <mol> --mode normal` first).
      - Combined flags run in fixed order: --library, --calibrate,
        --emit-projection, --ped-merge, --ped-merge-all, --figures. Fails
        loud on a bad combination rather than silently doing nothing.

    `thresholds` (the effective, possibly CLI-overridden Thresholds built in
    main()) is passed through to --emit-projection's run_projection_pipeline
    call; defaults to Thresholds.calibrated() if not given (e.g. a direct
    caller/test that doesn't build one itself).

    Returns True if at least one flag was handled (caller should stop),
    False otherwise (caller falls through to the plain user workflow).
    """
    thresholds = thresholds if thresholds is not None else Thresholds.calibrated()
    any_flag = (args.library or args.calibrate or args.emit_projection
                or args.figures or args.ped_merge or args.ped_merge_all)
    if not any_flag:
        return False

    if args.library:
        if args.molecule:
            print("Note: --library is global and ignores -m/--molecule.")
        print("Running library ingest (src.library_ingest.run_ingest_pipeline) "
              "-- this is known-slow. Please wait...")
        from src.library_ingest import run_ingest_pipeline
        df_lib, path, skip_report = run_ingest_pipeline()
        print(f"Wrote {len(df_lib)} rows -> {path}")
        if skip_report:
            print(f"  {len(skip_report)} molecule(s) had their internal-row ref_label/ref_key "
                  "join skipped (frequency mismatch vs. characterised_modes.csv) -- see warnings "
                  "above. Their scores (V_Stretch, Tx..Rz, predicted_label, ...) and their "
                  "'ideal' tag (sourced separately from mol_list_method.csv) are unaffected.")

    if args.calibrate:
        if args.molecule and not args.library:
            print("Note: --calibrate is global and ignores -m/--molecule.")
        print("Running threshold calibration (src.calibrate.run_calibration_pipeline)...")
        from src.calibrate import run_calibration_pipeline
        # Named distinctly from this function's own `thresholds` param (the
        # effective, possibly CLI-overridden Thresholds passed in from
        # main()) -- this is the freshly-recalibrated result, not to be
        # confused with or silently substituted for it.
        calibrated_thresholds, result, sweep_df, (path_json, path_sweep) = run_calibration_pipeline()
        print(f"Frozen thresholds tau_TR={calibrated_thresholds.tau_TR}, "
              f"tau_S={calibrated_thresholds.tau_S}, tau_B={calibrated_thresholds.tau_B} "
              f"-> {path_json}")
        print(f"Wrote {len(sweep_df)}-row sensitivity sweep -> {path_sweep}")

        print("Running tau_SB error-vs-threshold sweep "
              "(src.calibrate.run_tau_sb_error_analysis)...")
        from src.calibrate import run_tau_sb_error_analysis
        _sb_sweep_df, plateau_all, plateau_test, path_sb_sweep, plateau_sc = run_tau_sb_error_analysis()
        _lo_all, _hi_all, mid_all, err_all = plateau_all
        _lo_test, _hi_test, mid_test, err_test = plateau_test
        _lo_sc, _hi_sc, mid_sc, err_sc = plateau_sc
        print(f"Wrote {len(_sb_sweep_df)}-row tau_SB sensitivity sweep -> {path_sb_sweep}")
        print(f"Optimal tau_SB (min error, ADVISORY ONLY): "
              f"all={mid_all} (err={err_all:.4f}), test={mid_test} (err={err_test:.4f}), "
              f"single-centre={mid_sc} (err={err_sc:.4f}) -- "
              f"the active default stays tau_SB={calibrated_thresholds.tau_SB}; adopt via "
              f"--tau-sb <value> or by editing thresholds.json.")

    if args.emit_projection or args.ped_merge:
        if not args.molecule:
            print("Error: --emit-projection/--ped-merge require -m/--molecule.")
            return True

    if args.emit_projection:
        try:
            df, df_full, (path_emit, path_full, path_full_cart) = run_projection_pipeline(
                args.molecule, thresholds=thresholds)
        except FileNotFoundError as e:
            print(f"Error: --emit-projection for '{args.molecule}' needs data/EMIT/, the "
                  f"normal-mode log, AND an existing <mol>_EMIT.csv ({e})")
            return True
        except ValueError as e:
            print(f"Error running --emit-projection for '{args.molecule}': {e}")
            return True
        _print_table(df, f"{args.molecule} EMIT->normal-mode contributions (merged)")
        print(f"Merged {len(df)}-mode grouped contributions into -> {path_emit}")
        print("(per-reference-mode full detail is wide -- not echoed here; "
              f"see the CSV) Wrote {len(df_full)}-row full projection detail -> {path_full}")
        print(f"Wrote {len(df_full)}-row Cartesian-overlap comparison detail -> {path_full_cart}")

    if args.ped_merge:
        from ped.merge_ped_scores import merge_molecule_ped
        try:
            _df, vib_rows, path = merge_molecule_ped(args.molecule)
        except (FileNotFoundError, ValueError) as e:
            print(f"Error running --ped-merge for '{args.molecule}': {e}")
            return True
        print(f"Wrote {len(vib_rows)}-mode PED merge -> {path}")

    if args.ped_merge_all:
        if args.molecule and not args.ped_merge:
            print("Note: --ped-merge-all is global and ignores -m/--molecule.")
        from ped.merge_ped_scores import build_combined_table, _read_roster, _REPO_ROOT
        molecules = _read_roster(_REPO_ROOT)
        combined_df, per_molecule_paths = build_combined_table(molecules, _REPO_ROOT, skip_missing=True)
        out_path = args.combined_output or os.path.join("data", "results", "combined_ped_vs_scores.csv")
        combined_df.to_csv(out_path, index=False, float_format="%.4f")
        print(f"Merged PED for {len(per_molecule_paths)} of {len(molecules)} roster molecule(s) "
              f"(rest skipped -- see messages above for missing data/ved/<mol>.ved/.vdf).")
        print(f"Wrote {len(combined_df)}-row combined table -> {out_path}")

    if args.figures:
        print("Regenerating all manuscript figures (src.figures.regenerate_all)...")
        from src.figures import regenerate_all
        regenerate_all()

    return True


def main():
    ap = argparse.ArgumentParser(description="Calculate Molecular Mode Scores")

    user_group = ap.add_argument_group(
        "User workflow",
        "Score + classify one molecule's modes. This is the whole job for most users: "
        "pick a molecule and a mode type, get one CSV back.")
    user_group.add_argument("-m", "--molecule", help="Molecule name (without extension)")
    user_group.add_argument("--mode", choices=["normal", "emit"],
                             help="Run non-interactively with this mode type (skips the prompt). "
                                  "Writes data/results/<mol>_{normal,EMIT}.csv (scores + Mu/K/Irrep "
                                  "+ Steps 2-4 classification, all in one file).")
    user_group.add_argument("--v-weighting", choices=list(V_WEIGHTINGS),
                             default=DEFAULT_V_WEIGHTING, dest="v_weighting",
                             help="Bond weighting in the V-score (eq:vscore). 'mu' (default) "
                                  "weights each bond by its reduced mass mu_AB = m_A*m_B/(m_A+m_B), "
                                  "so a bond counts for as much as the kinetic energy its "
                                  "stretching motion carries. 'none' is the original unweighted "
                                  "definition. The two agree exactly for homoleptic AB_n molecules "
                                  "(mu cancels); they differ only where bond types are mixed. "
                                  "Recorded in data/results/thresholds.json -- scoring and "
                                  "thresholds must be produced under the same setting.")
    user_group.add_argument("--tau-tr", type=float, default=None, dest="tau_tr",
                             help="Override tau_TR for this run only -- starts from "
                                  "Thresholds.calibrated() and is never written back to "
                                  "data/results/thresholds.json. DIAGNOSTIC-ONLY as of the "
                                  "2026-08-25 restructuring: Step 3 (T/R identification) no "
                                  "longer gates on tau_TR at all, so this value no longer changes "
                                  "any output label -- it is kept only for descriptive/sensitivity "
                                  "reporting (src.calibrate's tau_TR sweep). Default: the "
                                  "calibrated value.")
    user_group.add_argument("--tau-s", type=float, default=None, dest="tau_s",
                             help="Override tau_S (three-way stretching bar) for this run only -- "
                                  "see --tau-tr. Only affects scheme=threeway (no longer the "
                                  "default; pass --classify-scheme threeway to use it).")
    user_group.add_argument("--tau-b", type=float, default=None, dest="tau_b",
                             help="Override tau_B (three-way bending bar) for this run only -- "
                                  "see --tau-tr. Only affects scheme=threeway's Step-2 vibrational "
                                  "vocabulary (no longer the default; pass --classify-scheme "
                                  "threeway to use it) -- as of the 2026-08-25 restructuring this "
                                  "is purely a Step-2 boundary, not a purity gate of any kind.")
    user_group.add_argument("--tau-sb", type=float, default=None, dest="tau_sb",
                             help="Override tau_SB (single-cutoff binary S/B split) for this run "
                                  "only -- see --tau-tr. Only affects scheme=binary.")
    user_group.add_argument("--classify-scheme", choices=["threeway", "binary"],
                             default="binary", dest="classify_scheme",
                             help="Step-2 vibrational-label vocabulary: 'binary' (default, "
                                  "paper-standard) forces every mode to 'S' or 'B' via the single "
                                  "tau_SB cutoff, never 'SB'; 'threeway' (explicit opt-in) may "
                                  "label a mode 'SB' (mixed stretch/bend) via the tau_S/tau_B "
                                  "split. Does not affect Step 3 (the optional T/R identification), "
                                  "which uses no scheme or threshold at all.")
    user_group.add_argument("--no-identify-tr", action="store_false", dest="identify_tr",
                             help="Skip Step 3 (the optional T/R identification) -- every mode's "
                                  "tr_label/tr_score are left blank, only the Step-2 vib_label is "
                                  "reported. Default: Step 3 runs (identify_tr=True), matching the "
                                  "flowchart's default path; every reproduce.py/headless run keeps "
                                  "this on.")

    dev_group = ap.add_argument_group(
        "Developer / maintainer workflow",
        "Diagnostics, cross-checks, and manuscript-support pipelines. Not needed for ordinary use.")
    dev_group.add_argument("--emit-projection", action="store_true", dest="emit_projection",
                            help="Project raw EMIT eigenvectors onto the normal-mode reference "
                                 "basis for -m <molecule>. Requires <mol>_EMIT.csv to already "
                                 "exist (run -m <mol> --mode emit first); merges C2_Tx..C2_VMix "
                                 "(mass-weighted overlap) into that file in place and writes "
                                 "<mol>_EMIT_full.csv (per-reference-mode detail). Also runs the "
                                 "plain-Cartesian-overlap comparison pathway (raw signed overlap, "
                                 "NOT squared), merging Ocart_Tx..Ocart_Sum and writing "
                                 "<mol>_EMIT_full_cartesian.csv "
                                 "-- see src/projection.py's docstring for why it is kept only as "
                                 "a comparison, not a replacement for the mass-weighted result.")
    dev_group.add_argument("--ped-merge", action="store_true", dest="ped_merge",
                            help="Merge real VEDA4 PED (data/ved/<mol>.ved+.vdf) into -m "
                                 "<molecule>'s <mol>_normal.csv (must already exist -- run "
                                 "-m <mol> --mode normal first). Appends "
                                 "PED_Stretch_pct/PED_Bend_pct/per-bond-type columns in place. "
                                 "Per-molecule; requires -m.")
    dev_group.add_argument("--ped-merge-all", action="store_true", dest="ped_merge_all",
                            help="Run --ped-merge for every molecule in data/mol_list_method.csv's "
                                 "roster, skipping (with a message) any missing its <mol>_normal.csv "
                                 "or data/ved/<mol>.ved/.vdf pair, and write one combined table "
                                 "(default data/results/combined_ped_vs_scores.csv, override with "
                                 "--combined-output). Global (ignores -m).")
    dev_group.add_argument("--combined-output", dest="combined_output", default=None,
                            help="Output path for the --ped-merge-all combined table "
                                 "(default data/results/combined_ped_vs_scores.csv).")
    dev_group.add_argument("--library", action="store_true",
                            help="Ingest the full mol_list_method.csv roster (real-engine recompute "
                                 "for every molecule's on-disk .log/.gjf pair) -> "
                                 "data/results/library_scores.csv. Global (ignores -m); slow.")
    dev_group.add_argument("--calibrate", action="store_true",
                            help="Calibrate tau_TR/tau_S/tau_B against the ingested library -> "
                                 "data/results/thresholds.json + tau_sensitivity_sweep.csv. Also "
                                 "runs the advisory tau_SB error-vs-threshold sweep -> "
                                 "tau_sb_sensitivity_sweep.csv (prints a suggested tau_SB; never "
                                 "overwrites the frozen tau_SB default -- see --tau-sb). "
                                 "Global (ignores -m).")
    dev_group.add_argument("--figures", action="store_true",
                            help="Regenerate all manuscript figures from data/results/*.csv -> "
                                 "data/figures/*.{pdf,png}. Global (ignores -m).")
    args = ap.parse_args()

    # Before anything scores: every downstream pipeline reads this default
    # rather than taking the variant as an argument.
    set_v_weighting(args.v_weighting)

    # Effective thresholds: start from the calibrated (frozen) values, apply
    # any per-run CLI overrides on top. Never written back to
    # data/results/thresholds.json -- these are per-run only.
    base = Thresholds.calibrated()
    thresholds = Thresholds(
        tau_TR=args.tau_tr if args.tau_tr is not None else base.tau_TR,
        tau_S=args.tau_s if args.tau_s is not None else base.tau_S,
        tau_B=args.tau_b if args.tau_b is not None else base.tau_B,
        tau_SB=args.tau_sb if args.tau_sb is not None else base.tau_SB,
        v_weighting=base.v_weighting,
    )
    if any(v is not None for v in (args.tau_tr, args.tau_s, args.tau_b, args.tau_sb)) \
            or args.classify_scheme != "binary" or not args.identify_tr:
        print(f"Note: threshold/scheme override active for this run only (not written to "
              f"thresholds.json): tau_TR={thresholds.tau_TR} tau_S={thresholds.tau_S} "
              f"tau_B={thresholds.tau_B} tau_SB={thresholds.tau_SB} "
              f"scheme={args.classify_scheme} identify_tr={args.identify_tr}")

    if _run_flag_pipelines(args, thresholds):
        return

    if not args.molecule:
        print("Error: -m/--molecule is required (unless using --library/--calibrate/--figures alone).")
        return
    mol_name = args.molecule

    if args.mode:
        mode_type = args.mode
    else:
        print("-------------------------------------------------------")
        print(f" Processing Molecule: {mol_name}")
        print("-------------------------------------------------------")
        print("Choose calculation type:")
        print("  1. Normal Modes (reads from data/logs/)")
        print("  2. EMIT Modes (reads from data/EMIT/)")
        mode_type = {"1": "normal", "2": "emit"}.get(input("Enter 1 or 2: ").strip())
        if mode_type is None:
            print("Invalid choice.")
            return

    try:
        raw, dirs = load_inputs(mol_name, mode_type, "data")
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")
        return

    # Interactive fallback: let the user add bonds to the intermediate file
    # then reload (headless run_scoring_pipeline raises instead); see _cache_is_fresh.
    if not raw["bonds"]:
        inter = intermediate_path(dirs, mol_name, mode_type)
        print("\n" + "!" * 70)
        print(" ATTENTION: No bonding information found.")
        print(f" Open {inter} and add bonds (e.g. '1 2' per line), then save.")
        print("!" * 70 + "\n")
        input(">> Press ENTER after saving the file...")
        raw = IntermediateIO.load(inter)

    try:
        df, output_file = run_scoring_pipeline(mol_name, mode_type, thresholds=thresholds,
                                                scheme=args.classify_scheme,
                                                identify_tr=args.identify_tr)
    except ValueError as e:
        print(f"Error: {e}")
        return

    _print_table(df, "Scoring + Classification Results")
    print(f"\nResults saved to: {output_file}")


if __name__ == "__main__":
    main()
