"""
Target Prediction based on the Similarity Ensemble Approach (SEA)

Purpose:
This script predicts potential protein targets for a given molecule (SMILES). 
It calculates the similarity between the query molecule and sets of known 
ligands for various targets. It then uses a statistical model (Extreme Value 
Distribution) to determine the significance of the observed similarity.

Reference:
Keiser MJ, Roth BL, Armbruster BN, Ernsberger P, Irwin JJ, Shoichet BK. 
Relating protein pharmacology by ligand chemistry. Nat Biotech 25 (2), 197-206 (2007).
"""

import pickle
import pandas as pd
import numpy as np
import math
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import os

class SEAPredictor:
    def __init__(self, db_path, fit_path):
        """
        Constructor to initialize the SEA Predictor.
        
        Args:
            db_path (str): Path to the pickled target database (dictionary mapping 
                           target names to lists of RDKit bit-vector fingerprints).
            fit_path (str): Path to the '.fit' file containing background model parameters.
        """
        print(f"loading database: {db_path} ...")
        with open(db_path, 'rb') as f:
            raw_db = pickle.load(f)
            self.target_db = raw_db["fingerprints"]
            self.metadata = raw_db["metadata"]
            
        self.n_targets = len(self.target_db)
        self.morgan_gen = AllChem.GetMorganGenerator(radius=2, fpSize=2048)

        print(f"loading background params: {fit_path} ...")
        self.fit_params, self.default_threshold = self._load_fit_params(fit_path)
        print(f"SEA set successfully. Default Tanimoto threshold from fit file: {self.default_threshold}")

    def _load_fit_params(self, fit_path):
        """
        Parses statistical parameters from a .fit file to calibrate the background model.
        
        The background model accounts for random chemical similarity. In SEA, the 
        distribution of raw similarity scores follows an Extreme Value Distribution (EVD).
        
        Args:
            fit_path (str): Path to the text file containing 'MU' (mean) and 'SIGMA' 
                            (standard deviation) fitting coefficients.
        
        Returns:
            dict: A dictionary containing lists of coefficients for MU and SIGMA.
        """
        params = {'mu': [], 'sigma': []}

        with open(fit_path, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if not line or line.startswith("#"): 
                    continue
                
                if line.startswith("TANI"):
                    tani_threshold = float(line.split('\t')[1])
                
                elif line.startswith("MU"):
                    params['mu'] = [float(x) for x in line.split('\t')[1:]]
                elif line.startswith("SIGMA"):
                    params['sigma'] = [float(x) for x in line.split('\t')[1:]]
        
        return params, tani_threshold

    def predict(self, smiles, threshold=None, n_bits=2048):
        """
        Predicts protein targets for a query molecule by comparing it against target sets.
        
        Statistical Principle:
        1. Calculate the Tanimoto coefficient (Tc) between the query and all ligands 
           in a target set.
        2. Sum all Tc values above the 'threshold' to get a 'Raw Score'.
        3. Estimate the background Mean (MU) and StdDev (SIGMA) based on the size 
           of the target set (z_prod) using the pre-loaded fit parameters.
        4. Calculate a Z-score: (RawScore - MU) / SIGMA.
        5. Convert the Z-score to a P-value using the Gumbel distribution (EVD) 
           formula: P = 1 - exp(-exp(-z*)), where z* is the normalized Z-score.
        
        Args:
            smiles (str): SMILES string of the query molecule.
            threshold (float): Minimum Tanimoto similarity to contribute to the Raw Score.
                               Must match the threshold used during background fitting.
            n_bits (int): Length of the Morgan Fingerprint bit vector.
            
        Returns:
            pd.DataFrame: Sorted results containing Target names, Raw Scores, Z-scores, 
                          P-values, and E-values.
        """
        # set threshold (use params in .fit as default, or you can set it manually)
        use_threshold = threshold if threshold is not None else self.default_threshold

        # calculate the fingerprint of query smiles
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES"}
        
        query_fp = self.morgan_gen.GetFingerprint(mol)

        results = []
        
        for target_name, db_fps in self.target_db.items():
            sims = DataStructs.BulkTanimotoSimilarity(query_fp, db_fps)
            sims_arr = np.array(sims)
            valid_sims = sims_arr[sims_arr >= use_threshold]
            
            raw_score = np.sum(valid_sims)
            
            if raw_score == 0:
                continue
                
            max_tc = np.max(sims_arr)
            n_ref = len(db_fps)
            z_prod = float(n_ref)
            
            mu_p = self.fit_params['mu']
            sigma_p = self.fit_params['sigma']
            
            # calculation of EVD params
            bg_mean = mu_p[0] * (z_prod ** mu_p[1]) + mu_p[2]
            bg_sigma = sigma_p[0] * (z_prod ** sigma_p[1]) + sigma_p[2]
            
            if bg_sigma < 1e-6: bg_sigma = 1e-6
            
            z_score = (raw_score - bg_mean) / bg_sigma
            
            # P = 1 - exp(-exp(-z*)) = -expm1(-exp(-z*))
            gamma_val = np.euler_gamma
            exponent_term = -z_score * math.pi / math.sqrt(6) - gamma_val
            
            try:
                p_value = -math.expm1(-math.exp(exponent_term))
            except OverflowError:
                p_value = 0.0 
            
            info = self.metadata.get(target_name, {"name": "N/A", "organism": "N/A", "description": "N/A"})
            
            if p_value < 1e-300:
                p_str = "< 1e-300"
            else:
                p_str = "{:.3e}".format(p_value)
            
            e_value = p_value * self.n_targets
            
            # E-value = P-value * target_nums
            if e_value < 1e-300:
                e_str = "< 1e-300"
            else:
                e_str = "{:.3e}".format(e_value)

            results.append({
                "Target_Key": target_name,
                "Symbol": info["name"],
                "Organism": info["organism"],
                "Description": info["description"],
                "Ligand_Count": n_ref,
                "MaxTC": round(max_tc, 3),
                "RawScore": round(raw_score, 3),
                "Z_score": round(z_score, 2),
                "P_value": p_str,
                "E_value": e_str,
                "P_value_raw": p_value 
            })
            
        if not results:
            return pd.DataFrame()
        
        df = pd.DataFrame(results)
        df = df.sort_values(by="P_value_raw", ascending=True)
        df = df.drop(columns=["P_value_raw"])
        
        return df
    
    def save_results(self, df, query_name):
        """
        Exports the prediction DataFrame to a CSV file called "{query_name}_sea_results.csv".
        """
        if df.empty:
            print(f"No results to save for {query_name}.")
            return

        safe_name = "".join([c if c.isalnum() else "_" for c in query_name])
        file_name = f"{safe_name[:50]}_sea_results.csv"
        
        df.to_csv(file_name, index=False)
        print(f"Results for '{query_name}' saved to: {file_name}")


if __name__ == "__main__":
    # set working dir
    DB_PATH = "target_fp_db.pkl" 
    FIT_PATH = "chembl27_rdkit_ecfp4.fit" 
    
    if os.path.exists(DB_PATH) and os.path.exists(FIT_PATH):
        predictor = SEAPredictor(DB_PATH, FIT_PATH)
        
        ## set your target compound
        compound_name = "Urolithin A"
        test_smiles = "C1=CC2=C(C=C1O)C(=O)OC3=C2C=CC(=C3)O"

        ## prediction
        print(f"\npredicting: {compound_name} ({test_smiles})")
        res_df = predictor.predict(test_smiles)
        
        ## print top 10 results
        print("\nprediction results (Top 10):")
        print(res_df.head(10).to_string(index=False))

        ## save full results
        predictor.save_results(res_df, compound_name)
    else:
        print("Error: Missing database or .fit file.")