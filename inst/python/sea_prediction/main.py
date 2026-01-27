"""
SEA Target Prediction Main Script

This script integrates the data download and SEA prediction functionalities.
It first ensures all necessary data files are downloaded, then performs
target prediction for the specified compound.

Usage:
    python main.py
"""

import os
from dataloader import download_data
from runsea import SEAPredictor


def main():
    ## Configuration

    # Google Drive File IDs
    DB_ID = "1SV7-vrG3BFL51Ie07jtEBaOlFv8pkK3D"
    FIT_ID = "1jviqHWp0UaXj2bvhO5WY_bnOUyqLs_mu"
    
    # Local file paths (the path you want to store the downloaded files)
    DB_PATH = "target_fp_db.pkl"
    FIT_PATH = "chembl27_rdkit_ecfp4.fit"
    
    # ========== User Parameters ==========
    OUTPUT_DIR = "./results"      # Output directory for result files
    TOP_N_PRINT = 10              # Number of top results to print in terminal
    
    # Download Data Files
    
    print("=" * 50)
    print("Step 1: Checking and downloading data files...")
    print("=" * 50)
    
    download_data(DB_PATH, DB_ID)
    download_data(FIT_PATH, FIT_ID)
    
    # Initialize SEA Predictor
    
    print("\n" + "=" * 50)
    print("Step 2: Initializing SEA Predictor...")
    print("=" * 50)
    
    if not os.path.exists(DB_PATH) or not os.path.exists(FIT_PATH):
        print("Error: Required data files are missing. Please check download.")
        return
    
    predictor = SEAPredictor(DB_PATH, FIT_PATH)
    
    # Run Prediction
    
    print("\n" + "=" * 50)
    print("Step 3: Running target prediction...")
    print("=" * 50)
    
    # Define your query compound here
    compound_name = "Urolithin A"
    test_smiles = "C1=CC2=C(C=C1O)C(=O)OC3=C2C=CC(=C3)O"
    
    print(f"\nCompound: {compound_name}")
    print(f"SMILES: {test_smiles}")
    
    results_df = predictor.predict(test_smiles)
    
    # Display and Save Results
    
    print("\n" + "=" * 50)
    print("Step 4: Results")
    print("=" * 50)
    
    if results_df.empty:
        print("No significant targets found.")
    else:
        # Display top N results
        print(f"\nTop {TOP_N_PRINT} predicted targets for {compound_name}:")
        print(results_df.head(TOP_N_PRINT).to_string(index=False))
        
        # Create output directory if not exists
        if not os.path.exists(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR)
            print(f"\nCreated output directory: {OUTPUT_DIR}")
        
        # Save full results to CSV
        safe_name = "".join([c if c.isalnum() else "_" for c in compound_name])
        output_file = os.path.join(OUTPUT_DIR, f"{safe_name[:50]}_sea_results.csv")
        results_df.to_csv(output_file, index=False)
        print(f"Results saved to: {output_file}")
        
        print(f"\nTotal targets found: {len(results_df)}")
    
    print("\n" + "=" * 50)
    print("Prediction complete!")
    print("=" * 50)


if __name__ == "__main__":
    main()
