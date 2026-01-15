import argparse
import pandas as pd
import os
import sys
from src.parser import GaussianParser, IntermediateIO
from src.scoring import ModeScorer

def main():
    parser = argparse.ArgumentParser(description="Calculate Molecular Mode Scores (J. Chem. Educ.)")
    parser.add_argument("input_file", help="Path to Gaussian .log file (in data/logs/)")
    parser.add_argument("--out", help="Optional explicit output path (overrides auto-naming)")
    
    args = parser.parse_args()
    
    # 1. Setup Systematic Paths
    base_name = os.path.splitext(os.path.basename(args.input_file))[0]
    
    # Define folder structure
    data_dir = "data"
    inter_dir = os.path.join(data_dir, "intermediate")
    results_dir = os.path.join(data_dir, "results")
    
    # Ensure folders exist
    os.makedirs(inter_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    
    # Define filenames
    intermediate_file = os.path.join(inter_dir, f"{base_name}.txt")
    
    if args.out:
        output_file = args.out
    else:
        output_file = os.path.join(results_dir, f"{base_name}_scores.csv")
    
    try:
        # Step 2: Parse the Log File
        print(f"Reading Gaussian log: {args.input_file}")
        gp = GaussianParser(args.input_file)
        raw_data = gp.parse()
        
        # Step 3: Save to Intermediate File
        print(f"Saving extracted data to: {intermediate_file}")
        IntermediateIO.save(raw_data, intermediate_file)
        
        # Step 4: Check if bonds exist. If not, PAUSE for user input.
        if not raw_data['bonds']:
            print("\n" + "!"*70)
            print(" ATTENTION: No bonding information found in the log file.")
            print(f" Please open the file below and manually add bonds under the 'BONDS' section:")
            print(f" FILE: {os.path.abspath(intermediate_file)}")
            print("\n Format: 'AtomIndex1 AtomIndex2' (e.g., '1 2' for bond between atom 1 and 2).")
            print("!"*70 + "\n")
            
            input(">> Once you have saved the file with bonds, press ENTER to continue...")
        
        # Step 5: Reload validated data from the intermediate file
        print(f"Loading data from: {intermediate_file}...")
        data = IntermediateIO.load(intermediate_file)
        
        if not data['bonds']:
            print("Warning: Still no bonds found. Vibrational scores (V) will be 0.")
        else:
            print(f"Loaded {len(data['bonds'])} bonds.")

        # Step 6: Calculate Scores
        scorer = ModeScorer(data['atoms'], data['coords'], data['bonds'])
        results = []
        
        print(f"Calculating scores for {len(data['modes'])} modes...")
        for i, mode in enumerate(data['modes']):
            scores = scorer.calculate_scores(mode['vector'])
            results.append({
                "Mode": i + 1,
                "Freq": mode['frequency'],
                "Tx": scores['T']['x'], "Ty": scores['T']['y'], "Tz": scores['T']['z'],
                "Rx": scores['R']['x'], "Ry": scores['R']['y'], "Rz": scores['R']['z'],
                "V_Stretch": scores['V']
            })

        # Step 7: Save Results
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