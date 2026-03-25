# SEA: Similarity Ensemble Approach for Target Prediction

A Python implementation of the Similarity Ensemble Approach (SEA) for predicting protein targets of small molecules based on chemical similarity.

## Background

The Similarity Ensemble Approach (SEA) is a computational method that relates proteins based on the chemical similarity of their ligands. Given a query molecule, SEA compares its structure against sets of known ligands for each protein target and uses statistical models to assess the significance of observed similarities.

The method employs the Extreme Value Distribution (EVD) to model background similarity distributions, enabling robust statistical inference of target predictions.

## Reference

Keiser MJ, Roth BL, Armbruster BN, Ernsberger P, Irwin JJ, Shoichet BK. Relating protein pharmacology by ligand chemistry. *Nature Biotechnology* 25 (2), 197-206 (2007). DOI: [10.1038/nbt1284](https://doi.org/10.1038/nbt1284)

## Repository Structure

```
sea/
├── main.py           # Main entry script
├── runsea.py         # Core SEA prediction algorithm
├── dataloader.py     # Data download utilities
├── requirements.txt  # Python dependencies
└── README.md
```

### File Descriptions

| File | Description |
|------|-------------|
| `main.py` | Integrated pipeline that handles data download, model initialization, and prediction execution |
| `runsea.py` | Contains the `SEAPredictor` class implementing the SEA algorithm with EVD-based statistical modeling |
| `dataloader.py` | Utility functions for downloading pre-computed fingerprint databases and fitting parameters from remote storage |

## Requirements

- Python >= 3.8
- numpy >= 1.20.0
- pandas >= 1.3.0
- rdkit >= 2022.03.1
- gdown >= 4.5.0

## Installation

1. Clone or download this repository

2. Install dependencies:

```bash
pip install -r requirements.txt
```

## Quick Start

### Using the Main Script

Edit the compound information in `main.py`:

```python
## use quercetin as example
compound_name = "quercetin"
test_smiles = "O=c1c(O)c(c2ccc(O)c(O)c2)oc2cc(O)cc(O)c12"
```

Run the prediction:

```bash
python main.py
```

The script will automatically download required data files on first run.


## Output Format

The prediction results contain the following columns:

| Column | Description |
|--------|-------------|
| Target_Key | Unique identifier for the target (ChEMBL ID) |
| Symbol | Gene symbol of the target protein |
| Organism | Source organism |
| Description | Functional description of the target |
| Ligand_Count | Number of known ligands in the target set |
| MaxTC | Maximum Tanimoto coefficient between query and target ligands |
| RawScore | Sum of Tanimoto coefficients above threshold |
| Z_score | Standardized score based on EVD background model |
| P_value | Statistical significance of the prediction |
| E_value | Expected number of false positives (P-value × number of targets) |

## Configuration Parameters

Parameters can be modified in `main.py`:

```python
OUTPUT_DIR = "./results"    # Output directory for CSV files
TOP_N_PRINT = 10            # Number of results to display in terminal
```

## Data Files

The following data files are required and will be downloaded automatically:

- `target_fp_db.pkl`: Pre-computed Morgan fingerprints (radius=2, 2048 bits) for known target-ligand associations derived from ChEMBL
- `chembl27_rdkit_ecfp4.fit`: Background model parameters including the Tanimoto threshold ($\tau$), mu coefficients, and sigma coefficients fitted using ChEMBL27 data

## Methodology

1. **Fingerprint Generation**: Morgan fingerprints (ECFP4-like, 2048 bits) are computed for the query molecule

2. **Similarity Calculation**: Tanimoto coefficients are computed between the query and all ligands in each target set

3. **Raw Score**: Sum of Tanimoto coefficients exceeding the threshold 

4. **Statistical Modeling**: Background distribution parameters are estimated based on target set size using pre-fitted EVD coefficients

5. **Significance Assessment**: Z-scores are converted to P-values using the Gumbel distribution; E-values provide multiple testing correction

## Statistical Framework

The statistical foundation of SEA is based on modeling chemical similarity as a stochastic process. When comparing a query molecule against a set of ligands, the distribution of similarity scores follows predictable patterns that can be modeled using Extreme Value Theory.

### Tanimoto Coefficient

The Tanimoto coefficient (Tc) measures the similarity between two molecular fingerprints:

$$T_c(A, B) = \frac{|A \cap B|}{|A \cup B|}$$

where A and B are the bit vectors of two molecular fingerprints. Values range from 0 (no similarity) to 1 (identical fingerprints).

### Raw Score Calculation

For a query molecule q and a target set T containing n ligands, the raw score is defined as:

$$S_{raw} = \sum_{i=1}^{n} T_c(q, l_i) \cdot \mathbf{1}[T_c(q, l_i) \geq \tau]$$

where $\tau$ is the similarity threshold and $\mathbf{1}[\cdot]$ is the indicator function.

### Extreme Value Distribution Model

Under the null hypothesis of random chemical similarity, raw scores follow an Extreme Value Distribution (EVD). The background distribution parameters are modeled as functions of the target set size:

$$\mu(n) = a_\mu \cdot n^{b_\mu} + c_\mu$$

$$\sigma(n) = a_\sigma \cdot n^{b_\sigma} + c_\sigma$$

where $n$ is the number of ligands in the target set, and the coefficients $(a, b, c)$ are pre-fitted from large-scale chemical databases.

### Z-score and P-value

The Z-score standardizes the raw score against the background:

$$Z = \frac{S_{raw} - \mu}{\sigma}$$

The P-value is computed using the Gumbel distribution (Type I EVD):

$$P = 1 - \exp\left(-\exp\left(-Z \cdot \frac{\pi}{\sqrt{6}} - \gamma\right)\right)$$

where $\gamma \approx 0.5772$ is the Euler-Mascheroni constant.

### E-value

The E-value provides a multiple testing correction by estimating the expected number of false positives:

$$E = P \times N_{targets}$$

where $N_{targets}$ is the total number of targets in the database. 

## Related Resources

- Official SEA Search Server: [https://sea.bkslab.org/](https://sea.bkslab.org/)