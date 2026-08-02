import argparse
import os
import numpy as np
import pandas as pd

from src.parser import GaussianParser, EMITParser, IntermediateIO
from src.scoring import ModeScorer
from src.classifier import classify_all_modes, classify_to_rows, is_linear
from src.projection import build_reference_basis, project_emit
from src.utils import find_file

# Column order for the results table / CSV.
_SCORE_COLS = ["Tx", "Ty", "Tz", "Rx", "Ry", "Rz", "V_Stretch"]


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


def build_scorer_and_final(raw, mode_type):
    """Align to principal axes and build the candidate mode pool ('final'),
    shared by score_modes() and the classifier so both score the identical
    mode list.
      - 'normal': MIT rotates the molecule AND the mode vectors
        (rotate_modes=True); the 3 ideal T + 3 ideal R references are
        prepended to the real vibrational modes.
      - 'emit'  : EMIT modes already live in principal axes, so MIT rotates
        only the molecule (rotate_modes=False); no ideal references are
        added -- all 3N raw EMIT eigenvectors are the candidate pool.

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
    scorer = ModeScorer(raw["atoms"], raw["coords"], raw["bonds"])
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


def score_modes(raw, mode_type):
    """Score every mode in build_scorer_and_final's candidate pool; returns a list of result-row dicts."""
    scorer, final = build_scorer_and_final(raw, mode_type)

    rows = []
    for i, mode in enumerate(final):
        sc = scorer.calculate_scores(mode["vector"])
        is_emit = mode.get("is_emit", False)
        rows.append({
            "Mode": mode.get("label", f"Mode {i+1}"),
            ("Eigenvalue" if is_emit else "Freq"): mode["frequency"],
            "Tx": sc["T"]["x"], "Ty": sc["T"]["y"], "Tz": sc["T"]["z"],
            "Rx": sc["R"]["x"], "Ry": sc["R"]["y"], "Rz": sc["R"]["z"],
            "V_Stretch": sc["V"],
            # Only real Gaussian normal modes carry mu/k/irrep; .get() yields
            # None for EMIT modes and the synthetic ideal T/R references.
            "Mu": mode.get("reduced_mass"),
            "K": mode.get("force_constant"),
            "Irrep": mode.get("irrep"),
        })
    return rows


def run_classify_pipeline(mol_name, mode_type, data_dir="data", thresholds=None, write=True):
    """Headless pipeline: load inputs -> classify_all_modes (Algorithm 1) -> CSV."""
    raw, dirs = load_inputs(mol_name, mode_type, data_dir)
    scorer, final = build_scorer_and_final(raw, mode_type)
    scored = classify_all_modes(scorer, final, thresholds)
    df = pd.DataFrame(classify_to_rows(scored))
    suffix = "normal" if mode_type == "normal" else "EMIT"
    output_file = os.path.join(dirs["results"], f"{mol_name}_{suffix}_classified.csv")
    if write:
        df.to_csv(output_file, index=False, float_format="%.4f")
    return df, output_file


def run_projection_pipeline(mol_name, data_dir="data", thresholds=None, write=True):
    """Headless EMIT -> normal-mode projection pipeline (Phase 2, eq:emitproj).

    Builds the mass-weighted normal-mode reference basis from data/logs/<mol>.log,
    projects the raw EMIT eigenvectors from data/EMIT/<mol>_EMIT.txt onto it, and
    writes two files: <mol>_EMIT_contributions.csv (grouped T/R/V fractions per
    EMIT mode) and <mol>_EMIT_projection_full.csv (per-reference-mode detail).
    Returns (df_grouped, df_full, (path_grouped, path_full)).
    """
    raw_n, dirs = load_inputs(mol_name, "normal", data_dir)
    scorer_n, final_n = build_scorer_and_final(raw_n, "normal")
    raw_e, _ = load_inputs(mol_name, "emit", data_dir)
    scorer_e, final_e = build_scorer_and_final(raw_e, "emit")

    # Q and Theta are only frame-consistent because both parses read the SAME
    # log geometry and MIT()'s eigendecomposition is deterministic -- assert
    # it rather than relying on that silently.
    coords_n = np.array([[a.x(), a.y(), a.z()] for a in scorer_n.atoms])
    coords_e = np.array([[a.x(), a.y(), a.z()] for a in scorer_e.atoms])
    if not np.allclose(coords_n, coords_e, atol=1e-6):
        raise ValueError(
            f"'{mol_name}': the 'normal' and 'emit' parses of the log geometry "
            "rotated into different principal-axis frames (max coord diff "
            f"{np.max(np.abs(coords_n - coords_e)):.3g}) -- Q and Theta would "
            "not be projection-comparable. Check for a second/different "
            "'Standard orientation' block or non-deterministic MIT() sign fix."
        )

    ref = build_reference_basis(scorer_n, final_n, thresholds)
    rows, full_rows = project_emit(ref, final_e)

    df = pd.DataFrame(rows)
    df_full = pd.DataFrame(full_rows)
    out_grouped = os.path.join(dirs["results"], f"{mol_name}_EMIT_contributions.csv")
    out_full = os.path.join(dirs["results"], f"{mol_name}_EMIT_projection_full.csv")
    if write:
        df.to_csv(out_grouped, index=False, float_format="%.6f")
        df_full.to_csv(out_full, index=False, float_format="%.6f")
    return df, df_full, (out_grouped, out_full)


def _print_table(df, title):
    """Echo a results DataFrame via pandas to_string() -- borderless and paste-friendly into Excel."""
    print(f"\n--- {title} ---")
    print(df.to_string(index=False, float_format="%.4f"))


def _run_flag_pipelines(args):
    """Handle the --classify / --emit-projection / --library / --calibrate /
    --figures CLI flags by wiring up the existing headless pipeline functions.

      - --library/--calibrate/--figures are GLOBAL and ignore -m/--molecule.
      - --classify/--emit-projection are PER-MOLECULE and require -m;
        --classify additionally requires --mode (normal|emit).
      - Combined flags run in fixed order: --library, --calibrate, --classify,
        --emit-projection, --figures. Fails loud on a bad combination rather
        than silently doing nothing.

    Returns True if at least one flag was handled (caller should stop),
    False otherwise (caller falls through to the interactive/--mode path).
    """
    any_flag = args.library or args.calibrate or args.classify or args.emit_projection or args.figures
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
        thresholds, result, sweep_df, (path_json, path_sweep) = run_calibration_pipeline()
        print(f"Frozen thresholds tau_TR={thresholds.tau_TR}, tau_S={thresholds.tau_S}, "
              f"tau_B={thresholds.tau_B} -> {path_json}")
        print(f"Wrote {len(sweep_df)}-row sensitivity sweep -> {path_sweep}")

    if args.classify or args.emit_projection:
        if not args.molecule:
            print("Error: --classify/--emit-projection require -m/--molecule.")
            return True

    if args.classify:
        if not args.mode:
            print("Error: --classify requires --mode {normal,emit} "
                  "(it determines which classified CSV gets written).")
            return True
        try:
            df, path = run_classify_pipeline(args.molecule, args.mode)
        except (FileNotFoundError, ValueError) as e:
            print(f"Error running --classify for '{args.molecule}' ({args.mode}): {e}")
            return True
        _print_table(df, f"{args.molecule} classification ({args.mode})")
        print(f"Wrote {len(df)}-row classification -> {path}")

    if args.emit_projection:
        try:
            df, df_full, (path_grouped, path_full) = run_projection_pipeline(args.molecule)
        except FileNotFoundError as e:
            print(f"Error: --emit-projection for '{args.molecule}' needs both normal-mode "
                  f"AND EMIT input data present ({e})")
            return True
        except ValueError as e:
            print(f"Error running --emit-projection for '{args.molecule}': {e}")
            return True
        _print_table(df, f"{args.molecule} EMIT->normal-mode contributions")
        print(f"Wrote {len(df)}-row grouped contributions -> {path_grouped}")
        print("(per-reference-mode full detail is wide -- not echoed here; "
              f"see the CSV) Wrote {len(df_full)}-row full projection detail -> {path_full}")

    if args.figures:
        print("Regenerating all manuscript figures (src.figures.regenerate_all)...")
        from src.figures import regenerate_all
        regenerate_all()

    return True


def main():
    ap = argparse.ArgumentParser(description="Calculate Molecular Mode Scores")
    ap.add_argument("-m", "--molecule", help="Molecule name (without extension)")
    ap.add_argument("--mode", choices=["normal", "emit"],
                    help="Run non-interactively with this mode type (skips the prompt).")
    ap.add_argument("--classify", action="store_true",
                    help="Run Steps 2-4 classification for -m <molecule> "
                         "(requires --mode). Writes <mol>_{normal,emit}_classified.csv.")
    ap.add_argument("--emit-projection", action="store_true", dest="emit_projection",
                    help="Project raw EMIT eigenvectors onto the normal-mode reference "
                         "basis for -m <molecule>. Writes <mol>_EMIT_contributions.csv "
                         "and <mol>_EMIT_projection_full.csv.")
    ap.add_argument("--library", action="store_true",
                    help="Ingest the full mol_list_method.csv roster (real-engine recompute "
                         "for every molecule's on-disk .log/.gjf pair) -> "
                         "data/results/library_scores.csv. Global (ignores -m); slow.")
    ap.add_argument("--calibrate", action="store_true",
                    help="Calibrate tau_TR/tau_S/tau_B against the ingested library -> "
                         "data/results/thresholds.json + tau_sensitivity_sweep.csv. Global (ignores -m).")
    ap.add_argument("--figures", action="store_true",
                    help="Regenerate all manuscript figures from data/results/*.csv -> "
                         "data/figures/*.{pdf,png}. Global (ignores -m).")
    args = ap.parse_args()

    if _run_flag_pipelines(args):
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
    # then reload (headless run_pipeline/score_modes raise instead); see _cache_is_fresh.
    if not raw["bonds"]:
        inter = intermediate_path(dirs, mol_name, mode_type)
        print("\n" + "!" * 70)
        print(" ATTENTION: No bonding information found.")
        print(f" Open {inter} and add bonds (e.g. '1 2' per line), then save.")
        print("!" * 70 + "\n")
        input(">> Press ENTER after saving the file...")
        raw = IntermediateIO.load(inter)

    try:
        rows = score_modes(raw, mode_type)
    except ValueError as e:
        print(f"Error: {e}")
        return

    df = pd.DataFrame(rows)
    suffix = "normal" if mode_type == "normal" else "EMIT"
    output_file = os.path.join(dirs["results"], f"{mol_name}_{suffix}_scores.csv")
    _print_table(df, "Scoring Results")
    df.to_csv(output_file, index=False, float_format="%.4f")
    print(f"\nResults saved to: {output_file}")


if __name__ == "__main__":
    main()
