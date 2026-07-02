import argparse
import os
import numpy as np
import pandas as pd

from src.parser import GaussianParser, EMITParser, IntermediateIO
from src.scoring import ModeScorer
from src.classifier import classify_all_modes, classify_to_rows, is_linear
from src.projection import build_reference_basis, project_emit

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
    for ext in (".log", ".out"):
        p = os.path.join(logs_dir, f"{mol_name}{ext}")
        if os.path.exists(p):
            return p
    raise FileNotFoundError(f"No log file for '{mol_name}' in {logs_dir} (expected .log or .out)")


def load_inputs(mol_name, mode_type, data_dir="data"):
    """Parse geometry, modes, and connectivity for a molecule.

    mode_type is 'normal' (Gaussian vibrational modes) or 'emit' (EMIT modes).
    Returns (raw_data, dirs). Raises on missing files or wrong mode counts (the
    parsers fail loud). Does NOT prompt; suitable for headless / reproduce use.
    """
    if mode_type not in ("normal", "emit"):
        raise ValueError(f"mode_type must be 'normal' or 'emit', got {mode_type!r}")
    dirs = resolve_dirs(data_dir)
    gp = GaussianParser(_find_log(dirs["logs"], mol_name))
    if mode_type == "normal":
        raw = gp.parse(parse_modes=True)
    else:
        raw = gp.parse(parse_modes=False)
        emit_path = os.path.join(dirs["EMIT"], f"{mol_name}_EMIT.txt")
        if not os.path.exists(emit_path):
            emit_path = os.path.join(dirs["EMIT"], f"{mol_name}.txt")
        if not os.path.exists(emit_path):
            raise FileNotFoundError(
                f"No EMIT file for '{mol_name}' in {dirs['EMIT']} (expected {mol_name}_EMIT.txt)")
        raw["modes"] = EMITParser(emit_path, len(raw["atoms"])).parse()
    return raw, dirs


def build_scorer_and_final(raw, mode_type):
    """Align to principal axes and build the candidate mode pool ('final').

    Shared by score_modes() and classifier.run_classification() so both stages
    construct the exact same mode list from the same raw parse:
      - 'normal': ModeScorer.MIT rotates the molecule AND the mode vectors, and
        the 3 ideal translations + 3 ideal rotations are prepended (they are
        literally in the candidate pool, alongside the real vibrational modes).
      - 'emit'  : MIT rotates the molecule only (EMIT modes already live in
        principal axes); no ideal references are added -- all 3N raw EMIT
        eigenvectors are the candidate pool.

    Raises ValueError if no bond connectivity is available (a missing bond list
    silently corrupts the V-score, so we fail loud rather than score garbage).
    Returns (scorer, final) where final is a list of mode dicts.
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
            # n_R = 2 for a linear molecule (spec). construct_R() always builds
            # 3 ideal references, but MIT() places the linear (smallest-moment)
            # axis on the new X axis, so the "Rx" reference is an all-zero
            # vector -- an ill-defined placeholder, not a genuine external
            # mode. Drop it before it enters the candidate pool: left in, it
            # would fall through classify_all_modes' Step 2 (no slot claims
            # it, since external_slots() correctly excludes Rx for linear
            # molecules) into Step 4 and be mislabeled BENDING (V=0 <= tau_B).
            ideal_R = [m for m in ideal_R if m["label"] != "Rx"]
        final = scorer.construct_T() + ideal_R + rotated
    else:
        final = rotated
    return scorer, final


def score_modes(raw, mode_type):
    """Pure scoring core: align to principal axes, build ideal T/R modes (normal
    only), and score every mode. Returns a list of result-row dicts.

    Raises ValueError if no bond connectivity is available (see
    build_scorer_and_final).
    """
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
        })
    return rows


def run_pipeline(mol_name, mode_type, data_dir="data", write=True):
    """Headless end-to-end pipeline: load inputs -> score -> optionally write CSV.

    Returns (DataFrame, output_path). Raises on any missing input, missing bonds,
    or wrong mode count. No prompts, so reproduce.py and tests can call it directly.
    """
    raw, dirs = load_inputs(mol_name, mode_type, data_dir)
    df = pd.DataFrame(score_modes(raw, mode_type))
    suffix = "normal" if mode_type == "normal" else "EMIT"
    output_file = os.path.join(dirs["results"], f"{mol_name}_{suffix}_scores.csv")
    if write:
        df.to_csv(output_file, index=False, float_format="%.4f")
    return df, output_file


def run_classify_pipeline(mol_name, mode_type, data_dir="data", thresholds=None, write=True):
    """Headless classify pipeline: load inputs -> classify_all_modes -> CSV.

    Mirrors run_pipeline() but runs the full Algorithm 1 classifier
    (src/classifier.py) instead of stopping at Step-1 scores. Returns
    (DataFrame, output_path); raises on any missing input, missing bonds, or
    wrong mode count (same fail-loud behaviour as load_inputs()/
    build_scorer_and_final()).
    """
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

    Builds the mass-weighted normal-mode reference basis (ideal T/R + real
    vibrational modes, src/projection.py's locked convention) from
    data/logs/<mol>.log, projects the raw EMIT eigenvectors from
    data/EMIT/<mol>_EMIT.txt onto it, and writes two files:
      - data/results/<mol>_EMIT_contributions.csv -- grouped fractions
        (C2_Tx..C2_Rz, C2_VS/VB/VMix) per EMIT mode; matches the columns and
        semantics of the pre-existing hand-derived ground-truth file of the
        same name (validated to ~1e-4 absolute agreement on EMIT 2/9/34-36).
      - data/results/<mol>_EMIT_projection_full.csv -- the finer-grained
        per-individual-reference-mode Theta_tilde**2 detail (one column per
        ideal T/R slot and per real normal mode), the richer
        "projection-coefficients data file" the roadmap calls for.

    Raises on any missing input / missing bonds / wrong mode count (same
    fail-loud behaviour as load_inputs()/build_scorer_and_final()).
    Returns (df_grouped, df_full, (path_grouped, path_full)).
    """
    raw_n, dirs = load_inputs(mol_name, "normal", data_dir)
    scorer_n, final_n = build_scorer_and_final(raw_n, "normal")
    raw_e, _ = load_inputs(mol_name, "emit", data_dir)
    scorer_e, final_e = build_scorer_and_final(raw_e, "emit")

    # Q (built from scorer_n's rotated geometry) and Theta (the raw, unrotated
    # EMIT eigenvectors) are only frame-consistent because both parses read
    # the SAME log geometry, and MIT()'s eigendecomposition is deterministic
    # on identical input, so scorer_n and scorer_e land in the identical
    # principal-axis frame. That coupling was implicit; assert it rather than
    # relying on parsing determinism silently (formula-auditor finding,
    # 2026-07-01).
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


def _run_flag_pipelines(args):
    """Handle the --classify / --emit-projection / --library / --calibrate /
    --figures flags (Phase-5 CLI subcommands). These wire up the EXISTING
    headless pipeline functions (run_classify_pipeline, run_projection_pipeline,
    src.excel_ingest.run_ingest_pipeline, src.calibrate.run_calibration_pipeline,
    src.figures.regenerate_all) -- no scoring/classification logic lives here.

    Design (documented in README.md "How to Use" / IMPLEMENTATION_PLAN.md):
      - --library/--calibrate/--figures are GLOBAL (molecule-independent) and
        ignore -m/--molecule if it happens to be supplied alongside them.
      - --classify/--emit-projection are PER-MOLECULE and require -m.
        --classify additionally requires --mode (normal|emit), since that is
        exactly what determines which classified CSV gets written.
      - All flags given together run in a fixed, sensible order: --library,
        --calibrate, --classify, --emit-projection, --figures (library/
        calibrate feed figures; classify/emit-projection are independent of
        each other and of library/calibrate). Any subset can be combined.
      - Fails loud (prints a clear message, returns without a stack trace) on
        a bad combination (e.g. --classify with no --mode, --emit-projection
        with no EMIT data) rather than silently doing nothing.

    Returns True if at least one flag was handled (caller should stop --
    the plain Step-1 scoring path is skipped), False if none of these flags
    were passed (caller falls through to the original interactive/--mode path).
    """
    any_flag = args.library or args.calibrate or args.classify or args.emit_projection or args.figures
    if not any_flag:
        return False

    if args.library:
        if args.molecule:
            print("Note: --library is global and ignores -m/--molecule.")
        print("Running library ingest (src.excel_ingest.run_ingest_pipeline) -- "
              "this is known-slow (~630s to parse the workbook via openpyxl). Please wait...")
        from src.excel_ingest import run_ingest_pipeline
        df_lib, path, skip_report = run_ingest_pipeline()
        print(f"Wrote {len(df_lib)} rows -> {path}")
        if skip_report:
            print(f"  {len(skip_report)} molecule(s) had their internal-row geometry "
                  "merge skipped (frequency mismatch) -- see warnings above.")

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
        print(f"Wrote {len(df)}-row grouped contributions -> {path_grouped}")
        print(f"Wrote {len(df_full)}-row full projection detail -> {path_full}")

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
                    help="Ingest the hydride-library spreadsheet -> "
                         "data/results/library_scores.csv. Global (ignores -m); slow (~630s).")
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

    # Interactive fallback: if connectivity is missing, let the user add bonds to
    # the intermediate file, then reload. (run_pipeline/score_modes raise instead.)
    if not raw["bonds"]:
        inter = os.path.join(dirs["intermediate"], f"{mol_name}_data.txt")
        IntermediateIO.save(raw, inter)
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
    print("\n--- Scoring Results ---")
    print(df.to_string(index=False, float_format="%.3f"))
    df.to_csv(output_file, index=False, float_format="%.4f")
    print(f"\nResults saved to: {output_file}")


if __name__ == "__main__":
    main()
