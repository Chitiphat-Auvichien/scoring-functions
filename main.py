import argparse
import pandas as pd
import os
import sys
from src.parser import GaussianParser, EMITParser, IntermediateIO
from src.scoring import ModeScorer

def main():
    parser = argparse.ArgumentParser(description="Calculate Molecular Mode Scores")
    parser.add_argument("-m", "--molecule", help="Molecule name (without extension)", required=True)
    args = parser.parse_args()
    
    mol_name = args.molecule
    
    data_dir = "data"
    logs_dir = os.path.join(data_dir, "logs")
    emit_dir = os.path.join(data_dir, "EMIT")
    gjf_dir = os.path.join(data_dir, "gjf")
    inter_dir = os.path.join(data_dir, "intermediate")
    results_dir = os.path.join(data_dir, "results")
    
    for d in [logs_dir, emit_dir, gjf_dir, inter_dir, results_dir]:
        os.makedirs(d, exist_ok=True)

    print("-------------------------------------------------------")
    print(f" Processing Molecule: {mol_name}")
    print("-------------------------------------------------------")
    print("Choose calculation type:")
    print("  1. Normal Modes (reads from data/logs/)")
    print("  2. EMIT Modes (reads from data/EMIT/)")
    
    choice = input("Enter 1 or 2: ").strip()
    if choice not in ['1', '2']:
        print("Invalid choice.")
        return

    # Find Log File
    log_path = None
    for ext in ['.log', '.out']:
        p = os.path.join(logs_dir, f"{mol_name}{ext}")
        if os.path.exists(p):
            log_path = p
            break
            
    if not log_path:
        print(f"Error: Could not find log file for '{mol_name}' in {logs_dir}")
        return

    intermediate_file = os.path.join(inter_dir, f"{mol_name}_data.txt")
    
    # Common variables to be set by logic blocks
    raw_data = None
    should_rotate_modes = True # Default for Normal Modes

    if choice == '1':
        # --- NORMAL MODES ---
        output_file = os.path.join(results_dir, f"{mol_name}_normal_scores.csv")
        print(f"Reading Gaussian log: {log_path}")
        gp = GaussianParser(log_path)
        raw_data = gp.parse(parse_modes=True)
        should_rotate_modes = True

    elif choice == '2':
        # --- EMIT MODES ---
        emit_path = os.path.join(emit_dir, f"{mol_name}_EMIT.txt")
        if not os.path.exists(emit_path):
            emit_path = os.path.join(emit_dir, f"{mol_name}.txt")
            
        if not os.path.exists(emit_path):
            print(f"Error: Could not find EMIT file in {emit_dir} (expected {mol_name}_EMIT.txt)")
            return

        output_file = os.path.join(results_dir, f"{mol_name}_EMIT_scores.csv")

        # 1. Parse Geometry only from log
        print(f"Reading Geometry from: {log_path}")
        gp = GaussianParser(log_path)
        raw_data = gp.parse(parse_modes=False)
        
        # 2. Parse EMIT modes
        print(f"Reading EMIT modes from: {emit_path}")
        ep = EMITParser(emit_path, len(raw_data['atoms']))
        emit_modes = ep.parse()
        raw_data['modes'] = emit_modes
        
        # EMIT modes are already in Principal Axes, so do NOT rotate them
        should_rotate_modes = False

    # --- COMMON WORKFLOW ---
    try:
        print(f"Saving extracted data to: {intermediate_file}")
        IntermediateIO.save(raw_data, intermediate_file)
        
        if not raw_data['bonds']:
            print("\n" + "!"*70)
            print(" ATTENTION: No bonding information found.")
            print(f" Please open {intermediate_file} and add bonds manually.")
            print(" Format: '1 2' (for bond between Atom 1 and Atom 2)")
            print("!"*70 + "\n")
            input(">> Press ENTER after saving the file...")
        
        print(f"Loading data from: {intermediate_file}...")
        data = IntermediateIO.load(intermediate_file)
        print(f"Loaded {len(data['bonds'])} bonds.")

        scorer = ModeScorer(data['atoms'], data['coords'], data['bonds'])
        
        # A. Rotate Molecule (Always) and Modes (Conditional)
        print("Aligning molecule to Principal Axes (MIT)...")
        if not should_rotate_modes:
            print("  -> EMIT modes detected: Rotating MOLECULE ONLY, preserving mode vectors.")
        else:
            print("  -> Normal modes detected: Rotating BOTH molecule and mode vectors.")
            
        rotated_modes = scorer.MIT(data['modes'], rotate_modes=should_rotate_modes)
        
        final_modes = []
        
        if choice == '1':
            print("Constructing ideal T/R modes...")
            t_modes = scorer.construct_T()
            r_modes = scorer.construct_R()
            for i, m in enumerate(rotated_modes):
                m['label'] = f"Vib {i+1}"
            final_modes = t_modes + r_modes + rotated_modes
        else:
            # For EMIT Modes: Use 3N modes directly
            final_modes = rotated_modes

        results = []
        print(f"Calculating scores for {len(final_modes)} modes...")
        
        for i, mode in enumerate(final_modes):
            scores = scorer.calculate_scores(mode['vector'])
            mode_name = mode.get('label', f"Mode {i+1}")
            
            is_emit = mode.get("is_emit", False)
            val_header = "Eigenvalue" if is_emit else "Freq"
            val_data = mode['frequency']
            
            row = {
                "Mode": mode_name,
                val_header: val_data,
                "Tx": scores['T']['x'], "Ty": scores['T']['y'], "Tz": scores['T']['z'],
                "Rx": scores['R']['x'], "Ry": scores['R']['y'], "Rz": scores['R']['z'],
                "V_Stretch": scores['V']
            }
            results.append(row)

        df = pd.DataFrame(results)
        print("\n--- Scoring Results ---")
        print(df.to_string(index=False, float_format="%.3f"))
        
        df.to_csv(output_file, index=False, float_format="%.4f")
        print(f"\nResults saved to: {output_file}")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()