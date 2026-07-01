import argparse
import os
import pandas as pd

from src.parser import GaussianParser, EMITParser, IntermediateIO
from src.scoring import ModeScorer

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


def score_modes(raw, mode_type):
    """Pure scoring core: align to principal axes, build ideal T/R modes (normal
    only), and score every mode. Returns a list of result-row dicts.

    Raises ValueError if no bond connectivity is available (a missing bond list
    silently corrupts the V-score, so we fail loud rather than score garbage).
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
        final = scorer.construct_T() + scorer.construct_R() + rotated
    else:
        final = rotated

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


def main():
    ap = argparse.ArgumentParser(description="Calculate Molecular Mode Scores")
    ap.add_argument("-m", "--molecule", required=True, help="Molecule name (without extension)")
    ap.add_argument("--mode", choices=["normal", "emit"],
                    help="Run non-interactively with this mode type (skips the prompt).")
    args = ap.parse_args()
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
