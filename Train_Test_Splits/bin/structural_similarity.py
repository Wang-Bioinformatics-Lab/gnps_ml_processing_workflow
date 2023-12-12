import argparse
import csv
import pandas as pd
import os

from joblib import Parallel, delayed
from pandarallel import pandarallel

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs.cDataStructs import BulkTanimotoSimilarity, CreateFromBitString

from tqdm import tqdm

PARALLEL_WORKERS = 32

def _generate_fingerprints_helper(smiles):
    if smiles is not None and smiles != 'nan':
        mol = Chem.MolFromSmiles(str(smiles), sanitize=True)
        if mol is not None:
            return smiles, {'Morgan_2048_2': list(AllChem.GetMorganFingerprintAsBitVect(mol,2,useChirality=False,nBits=2048)),
                            'Morgan_4096_2': list(AllChem.GetMorganFingerprintAsBitVect(mol,2,useChirality=False,nBits=4096)), 
                            'Morgan_2048_3': list(AllChem.GetMorganFingerprintAsBitVect(mol,3,useChirality=False,nBits=2048)), 
                            'Morgan_4096_3': list(AllChem.GetMorganFingerprintAsBitVect(mol,3,useChirality=False,nBits=4096))}
        else:
            raise NotImplementedError("Failed to parse smiles:", smiles)
    else:
        raise NotImplementedError("Smiles is None, or nan:", smiles)
    return smiles, {'Morgan_2048_2': None,
                    'Morgan_4096_2': None, 
                    'Morgan_2048_3': None, 
                    'Morgan_4096_3': None}

def generate_fingerprints(summary, progress=True):
    smiles_series = summary.Smiles.astype(str).str.strip()
    
    if progress:
        unique_smiles = tqdm(smiles_series.unique())
    else:
        unique_smiles = smiles_series.unique()
    mapping_dict = dict(Parallel(n_jobs=-1)(delayed(_generate_fingerprints_helper)(smiles) for smiles in unique_smiles))
    # Save mapping dict to text file
    with open('mapping_dict.txt', 'w') as f:
        f.write(str(mapping_dict))
        
    # Assert index is reset, this assignment will not work otherwise
    assert (summary.index == summary.reset_index().index).all()
    summary[['Morgan_2048_2','Morgan_4096_2','Morgan_2048_3','Morgan_4096_3']] = pd.DataFrame.from_records(smiles_series.map(mapping_dict))
    return summary

def compute_similarities_square(train_df:pd.DataFrame,
                                test_df:pd.DataFrame,
                                similarity_threshold=0.50,
                                fingerprint_column_name='Morgan_2048_3',) -> pd.DataFrame:
    if similarity_threshold < 0:
        raise ValueError("Expected arg 'similarity_threshold' to be between 0 and 1.")
    
    if not isinstance(train_df, pd.DataFrame):
        raise ValueError(f"Expected a Pandas DataFrame but got {type(train_df)}")
    if not isinstance(test_df, pd.DataFrame):
        raise ValueError(f"Expected a Pandas DataFrame but got {type(test_df)}")
       
    grouped_train_df = train_df.groupby('Smiles')
    train_df = train_df.drop_duplicates(subset='Smiles')
    
    grouped_test_df = test_df.groupby('Smiles')
    test_df = test_df.drop_duplicates(subset='Smiles')
    
    pandarallel.initialize(progress_bar=False, nb_workers=PARALLEL_WORKERS, verbose=0)
    train_fps   = train_df[fingerprint_column_name].parallel_apply(lambda x: CreateFromBitString(''.join(str(y) for y in x))).values
    test_fps    = test_df[fingerprint_column_name].parallel_apply(lambda x: CreateFromBitString(''.join(str(y) for y in x))).values
    train_fps   = list(train_fps)
    test_fps    = list(test_fps)

    # Indices are now non contiguous because the entries without structures are removed
    # This will map back to the original spectrum_id
    train_idx_mapping = {idx: group_df['spectrum_id'].values for idx, (_, group_df) in enumerate(grouped_train_df)}
    test_idx_mapping = {idx: group_df['spectrum_id'].values for idx, (_, group_df) in enumerate(grouped_test_df)}
    
    # Cleanup variables to reduce memory overhead
    del grouped_train_df
    del train_df
    del grouped_test_df
    del test_df
    
    fieldnames = ['spectrumid1', 'spectrumid2', 'Tanimoto_Similarity']
    
    output_file = "output.csv"
    with open(output_file, mode='w') as f:
        temp_writer = csv.DictWriter(f, fieldnames=fieldnames)
        temp_writer.writeheader()
    
        for train_index, train_fp in enumerate(train_fps):
            sims = BulkTanimotoSimilarity(train_fp, test_fps)
            
            # We only need the highest similarity (it's good enough to know it won't be in our train set)
            max_sim = -1.0
            max_sim_index = -1
            
            for j, this_sim in enumerate(sims): # Iterates over similarities to test structures
                if this_sim >= similarity_threshold and this_sim > max_sim:
                    max_sim = this_sim
                    max_sim_index = j
            if max_sim != -1.0:
                for edge_from in train_idx_mapping[train_index]:
                    # Note that this will output all test spectrum_ids that have an identical smiles. 
                    # This in not technically necessary but may be useful for debugging.
                    for edge_to in test_idx_mapping[max_sim_index]:
                        row = {
                            'spectrumid1': edge_from,
                            'spectrumid2': edge_to,  
                            'Tanimoto_Similarity': max_sim
                        }
                        temp_writer.writerow(row)


def generate_structural_similarity(train_csv, test_csv, similarity_threshold=0.7, fingerprint='Morgan_2048_3'):
    train_df = pd.read_csv(train_csv)
    test_df = pd.read_csv(test_csv)
    
    train_df = train_df[train_df['Smiles'].notna()].reset_index(drop=True)
    test_df = test_df[test_df['Smiles'].notna()].reset_index(drop=True)
    
    train_df = generate_fingerprints(train_df)
    test_df = generate_fingerprints(test_df)
    
    train_df.to_csv('train_df.csv')
    test_df.to_csv('test_df.csv')
    
    compute_similarities_square(train_df, test_df, similarity_threshold=similarity_threshold, fingerprint_column_name=fingerprint)

def main():
    parser = argparse.ArgumentParser(description='Generate Structural Similarity')
    parser.add_argument('--train_csv', type=str, help='Input CSV')
    parser.add_argument('--test_csv', type=str, help='Input CSV')
    parser.add_argument('--similarity_threshold', help='Similarity threshold for removing spectra from training set.', type=float)
    parser.add_argument('--fingerprint', help='Fingerprint to use for similarity calculations.', type=str, choices=['Morgan_2048_2','Morgan_4096_2','Morgan_2048_3','Morgan_4096_3'], default='Morgan_2048_3')
    args = parser.parse_args()
    
    generate_structural_similarity(args.train_csv, args.test_csv, args.similarity_threshold, fingerprint=args.fingerprint)

if __name__ == "__main__":
    main()
    
def test_structural_similarity():
    temp_file_path1 = 'test_structural_similarity_1.csv'
    temp_file_path2 = 'test_structural_similarity_2.csv'
    
    caffine_smiles  = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    asprin_smiles   = "CC(=O)OC1=CC=CC=C1C(=O)O"
    isobutane_smiles = "CC(C)C=O"
    biphenyl_phosphine_smiles = "c1ccc(cc1)c2ccccc2.P"
    
    
    df1 = pd.DataFrame({'spectrum_id':['spectrum_1','spectrum_2', 'spectrum_3', 'spectrum_4'], 'Smiles': [caffine_smiles, isobutane_smiles, isobutane_smiles,biphenyl_phosphine_smiles]})
    df2 = pd.DataFrame({'spectrum_id':['ref_spectrum_1','ref_spectrum_2', 'ref_spectrum_3','ref_spectrum_4'], 'Smiles': [caffine_smiles, asprin_smiles, isobutane_smiles, isobutane_smiles]})
    
    
    df1.to_csv(temp_file_path1)
    df2.to_csv(temp_file_path2)
    
    output_path = "output.csv"
    generate_structural_similarity(temp_file_path1,temp_file_path2,similarity_threshold=0)
    output_df = pd.read_csv(output_path)
    print(output_df)
    
    assert output_df.equals(pd.DataFrame({'spectrumid1':['spectrum_2','spectrum_3','spectrum_1','spectrum_4','spectrum_4'],
                                      'spectrumid2':['ref_spectrum_2','ref_spectrum_2','ref_spectrum_1','ref_spectrum_3','ref_spectrum_4'],
                                      'Tanimoto_Similarity':[1.0, 1.0, 1.0, 0.1538461538461538, 0.1538461538461538]}))
    
    if os.path.isfile(temp_file_path1):
        os.remove(temp_file_path1)
    if os.path.isfile(temp_file_path2):
        os.remove(temp_file_path2)
    if os.path.isfile(output_path):
        os.remove(output_path)