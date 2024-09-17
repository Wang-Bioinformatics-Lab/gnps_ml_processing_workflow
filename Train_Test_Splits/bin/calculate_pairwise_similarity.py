import argparse
import logging
from pathlib import Path

import pandas as pd
import numpy as np
from rdkit import Chem, DataStructs

from tqdm import tqdm


def compute_pairwise_similarities(metadata:Path, output_filename:Path):
    df = pd.read_csv(metadata)
    df['InChIKey_smiles_14'] = df['InChIKey_smiles'].str[:14]
    
    # Get unique structures by inchikey14, get most frequent corresponding smiles
    inchikey_smiles_mapping = df.groupby('InChIKey_smiles_14')['Smiles'].agg(lambda x: pd.Series.mode(x).iloc[0])
    unique_inchikeys = inchikey_smiles_mapping.index
    unique_smiles = inchikey_smiles_mapping.values

    logging.info("Computing Fingerprints")
    fps = [Chem.RDKFingerprint(Chem.MolFromSmiles(x)) for x in tqdm(unique_smiles)]

    logging.info("Computing Similarity Matrix")

    sim_matrix = np.zeros((len(fps), len(fps)))
    for i, fp1 in enumerate(tqdm(fps)):
        sim_matrix[i, i] = -1
        for j in range(i+1, len(fps)):
            fp2 = fps[j]
            sim_matrix[i, j] = DataStructs.TanimotoSimilarity(fp1, fp2)
            sim_matrix[j, i] = sim_matrix[i, j]

    # Set diagonal to 1.0
    np.fill_diagonal(sim_matrix, 1.0)
    
    sim_matrix = pd.DataFrame(sim_matrix, index=unique_inchikeys, columns=unique_inchikeys)
    sim_matrix.to_csv(output_filename)

def main():
    parser = argparse.ArgumentParser(description='Generate Structural Similarity Matrix')
    parser.add_argument('--input_metadata', type=str, help='Input CSV')
    parser.add_argument('--output_filename', type=str, help="Output CSV")
    parser.add_argument('--debug', action='store_true', help="Enable Debugging")
    args = parser.parse_args()

    logging_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(filename="calculate_pairwise_similarities.log", level=logging_level) 

    # Log all args
    for arg in vars(args):
        logging.info(f"{arg}: {getattr(args, arg)}")

    metadata = Path(args.input_metadata)
    output_filename = Path(args.output_filename)
    # Check if the input file exists
    if not metadata.exists():
        raise FileNotFoundError(f"File {metadata} does not exist")
    # Make output directories for the output file
    output_filename.parent.mkdir(parents=True, exist_ok=True)

    compute_pairwise_similarities(metadata, output_filename)

if __name__ == "__main__":
    main()