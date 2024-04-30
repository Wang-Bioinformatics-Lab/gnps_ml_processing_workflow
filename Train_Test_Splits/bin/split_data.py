import argparse
import pandas as pd
import numpy as np
from pyteomics.mgf import IndexedMGF
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdFMCS
from rdkit.Chem import DataStructs
from utils import calculate_murcko_histogram, calculate_murcko_histogram_distance
from tqdm import tqdm
import shutil

def synchronize_spectra(input_path, output_path, summary, progress_bar=True):
    """Reads an MGF file from input_path and generates a new_mgf file in output_path with only the spectra in spectrum_ids.

    Args:
        input_path (str): Path to the input mgf.
        output_path (str): Path to save the output mgf.
        summary (pd.Dataframe): Dataframe with spectrum ids to keep.
    """   
    for col_name in ['spectrum_id', 'scan', 'Charge']:
        if col_name not in summary.columns:
            raise ValueError("Summary must contain columns 'spectrum_id', 'scan', and 'charge'")
    
    with open(output_path, 'w') as output_mgf:
        input_mgf = IndexedMGF(input_path)
        
        if progress_bar:
            print("Syncing MGF with summary")
            mapping = tqdm(summary[['spectrum_id','scan','Charge']].itertuples())
        else:
            mapping = summary[['spectrum_id','scan','Charge']].itertuples()
        
        for _, title, scan, charge in mapping:
            spectra = input_mgf[title]
            if spectra['params']['title'] != title:
                raise ValueError("Sanity Check Failed. Expected specrum identifier did not match mgf spectrum identifier.")
            output_mgf.write("BEGIN IONS\n")
            output_mgf.write("PEPMASS={}\n".format(float(spectra['params']['precursor_mz'])))
            output_mgf.write("CHARGE={}\n".format(charge))
            output_mgf.write("TITLE={}\n".format(spectra['params']['title']))
            output_mgf.write("SCANS={}\n".format(scan))

            peaks = zip(spectra['m/z array'], spectra['intensity array'])
            for peak in peaks:
                output_mgf.write("{} {}\n".format(peak[0], peak[1]))

            output_mgf.write("END IONS\n")

def split_data_spectral(input_train_csv, input_train_mgf, input_test_csv, input_test_mgf, spectral_similarities, similarity_thresholds, progress_bar=True):
    
    shutil.copy(input_test_csv, input_test_csv[:-4] + f"_spectral.csv")
    shutil.copy(input_test_mgf, input_test_mgf[:-4] + f"_spectral.mgf")
    
    summary = pd.read_csv(input_train_csv)
    
    # Fasst search only outputs scan number in mgf, create a mapping scan -> spectrum_id
    spectra = IndexedMGF(input_train_mgf)
    scan_mapping = {int(spectrum['params']['scans']): spectrum['params']['title'] for spectrum in spectra}
    
    spectral_similarities = pd.read_table(spectral_similarities).drop(['Query Idx', 'Query Filename', 'Query Charge', 'Index UnitPM','Index IdxInUnitPM','Index Filename', 'Index Charge'], axis=1)
    spectral_similarities['Query_spectrumid'] = spectral_similarities['Query Scan'].apply(lambda x: scan_mapping[x])
    
    # Get the maximum number of data poitns to be dropped in order to hold train size fixed.
    min_threshold = min(similarity_thresholds)
    max_num_removed = len(spectral_similarities[spectral_similarities['Cosine'] >= min_threshold].Query_spectrumid.unique())
    min_num_points = len(summary) - max_num_removed
    
    for threshold in sorted(similarity_thresholds, reverse=True): # Sort in descending order to remove the most similar spectra first
        output_csv_path = input_train_csv[:-4] + f"_spectral_{threshold}.csv"
        output_mgf_path = input_train_mgf[:-4] + f"_spectral_{threshold}.mgf"
        
        to_drop = spectral_similarities.loc[spectral_similarities['Cosine'] <= threshold].Query_spectrumid.values
        summary = summary.loc[~summary['spectrum_id'].isin(to_drop)]
        
        # Subsample summary for saving
        temp = summary.sample(n=min_num_points, replace=False)
        temp.to_csv(output_csv_path)
        synchronize_spectra(input_train_mgf, output_mgf_path, temp)
        

def split_data_structural(input_train_csv, input_train_mgf, input_test_csv, input_test_mgf, structural_similarities, similarity_thresholds, progress_bar=True):
    # Column names are spectrumid1 spectrumid2 Tanimoto_Similarity
    
    def _generate_scaffold(smiles):
        smiles = rdMolStandardize.StandardizeSmiles(smiles)
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Could not generate mol from smiles: {smiles}")
        
        mol_scaffold_generic = MurckoScaffold.GetScaffoldForMol(mol)
        
        return Chem.CanonSmiles(Chem.MolToSmiles(mol_scaffold_generic))
    
    # Make sure to remove structures with no smiles strings since we can't split them
    train_summary = pd.read_csv(input_train_csv)
    train_summary = train_summary.loc[~ train_summary.Smiles.isna()]

    test_summary = pd.read_csv(input_test_csv)
    test_summary = test_summary.loc[~ test_summary.Smiles.isna()]
    
    output_test_csv_path  = input_test_csv[:-4] + f"_structural.csv"
    output_test_mgf_path  = input_test_mgf[:-4] + f"_structural.mgf"
       
    structural_similarities = pd.read_csv(structural_similarities)
    
    # Remove obvious ones where spectrum_id is shared
    all_test_ids = set(test_summary.spectrum_id.values)
    train_summary = train_summary.loc[~train_summary['spectrum_id'].isin(all_test_ids)]
    
    # Remove identical scaffolds from the training set
    # train_smiles = train_summary.Smiles.unique()
    # train_scaffold_map = {smiles: _generate_scaffold(smiles) for smiles in train_smiles}
    
    # test_smiles = test_summary.Smiles.unique()
    # test_scaffold_map = {smiles: _generate_scaffold(smiles) for smiles in test_smiles}
    
    # test_summary['scaffold'] = test_summary.Smiles.apply(lambda x: test_scaffold_map[x])
    # train_summary['scaffold'] = train_summary.Smiles.apply(lambda x: train_scaffold_map[x])
    
    test_summary.to_csv(output_test_csv_path)
    synchronize_spectra(input_test_mgf, output_test_mgf_path, test_summary)
    
    # # Remove shared scaffolds
    # shared_scaffolds = set(test_summary.scaffold.values).intersection(set(train_summary.scaffold.values))
    # print(f"Shared number of scaffolds: {len(shared_scaffolds)}")
    # removed_from_train = train_summary.loc[train_summary['scaffold'].isin(shared_scaffolds)]
    # print(f"Number of spectra with scaffolds in common with test {len(removed_from_train)}")
    
    # Get the maximum number of data poitns to be dropped in order to hold train size fixed.
    min_threshold = min(similarity_thresholds)
    # max_num_removed = len(set(structural_similarities[structural_similarities['Tanimoto_Similarity'] >= min_threshold].spectrumid1.unique())|set(removed_from_train.spectrum_id.unique()))
    max_num_removed = len(structural_similarities[structural_similarities['Tanimoto_Similarity'] >= min_threshold].spectrumid1.unique())
    min_num_points = len(train_summary) - max_num_removed
    
    # print(f"Length of Train Summary before Removing Intersecting Scaffolds: {len(train_summary)}")
    # train_summary = train_summary.loc[~train_summary['scaffold'].isin(shared_scaffolds)]
    # print(f"Length of Train Summary after Removing Intersecting Scaffolds: {len(train_summary)}")
    
    for threshold in sorted(similarity_thresholds, reverse=True): # Sort in descending order to remove the most similar spectra first
        output_train_csv_path = input_train_csv[:-4] + f"_structural_{threshold}.csv"
        output_train_mgf_path = input_train_mgf[:-4] + f"_structural_{threshold}.mgf"
        
        to_drop = structural_similarities.loc[structural_similarities['Tanimoto_Similarity'] >= threshold].spectrumid1.values
        train_summary = train_summary.loc[~train_summary['spectrum_id'].isin(to_drop)]
        
        
        # Subsample summary for saving
        temp = train_summary.sample(n=min_num_points, replace=False)
        temp.to_csv(output_train_csv_path)
        synchronize_spectra(input_train_mgf, output_train_mgf_path, temp)
        
def split_data_double_structural(input_train_csv, input_train_mgf, input_test_csv, input_test_mgf, structural_similarities, similarity_thresholds, progress_bar=True):
    # Column names are spectrumid1 spectrumid2 Tanimoto_Similarity
    
    # Make sure to remove structures with no smiles strings since we can't split them
    train_summary = pd.read_csv(input_train_csv)
    train_summary = train_summary.loc[~ train_summary.Smiles.isna()]

    test_summary = pd.read_csv(input_test_csv)
    test_summary = test_summary.loc[~ test_summary.Smiles.isna()]
    
    output_test_csv_path  = input_test_csv[:-4] + f"_double_structural.csv"
    output_test_mgf_path  = input_test_mgf[:-4] + f"_double_structural.mgf"
       
    structural_similarities = pd.read_csv(structural_similarities)
    
    # Compute the tanimoto similarities on daylight fingerprints
    unique_train_smiles = train_summary.Smiles.unique()
    unique_test_smiles = test_summary.Smiles.unique()
    
    # Generate the fingerprints
    train_fps = [Chem.RDKFingerprint(Chem.MolFromSmiles(smiles)) for smiles in unique_train_smiles]
    test_fps = [Chem.RDKFingerprint(Chem.MolFromSmiles(smiles)) for smiles in unique_test_smiles]
    
    # Compute the similarities
    similarities = np.zeros((len(train_fps), len(test_fps)))
    for i, train_fp in enumerate(train_fps):
        for j in range(i+1, len(test_fps)):
            test_fp = test_fps[j]
            similarities[i,j] = DataStructs.FingerprintSimilarity(train_fp, test_fp)
            similarities[j,i] = similarities[i,j]
    
    # Make into dataframe for convenience
    daylight_similarities = pd.DataFrame(similarities, index=unique_train_smiles, columns=unique_test_smiles)
    # Get max by row
    daylight_similarities = daylight_similarities.max(axis=1)
    
    # Remove obvious ones where spectrum_id is shared
    all_test_ids = set(test_summary.spectrum_id.values)
    train_summary = train_summary.loc[~train_summary['spectrum_id'].isin(all_test_ids)]
    
    
    test_summary.to_csv(output_test_csv_path)
    synchronize_spectra(input_test_mgf, output_test_mgf_path, test_summary)
    
    # Get the maximum number of data poitns to be dropped in order to hold train size fixed.
    min_threshold = min(similarity_thresholds)
    
    disallowed_smiles_daylight = daylight_similarities[daylight_similarities >= min_threshold].index.values
    disallowed_spectra_daylight = train_summary.loc[train_summary['Smiles'].isin(disallowed_smiles_daylight)]
    # max_num_removed = len(set(structural_similarities[structural_similarities['Tanimoto_Similarity'] >= min_threshold].spectrumid1.unique())|set(removed_from_train.spectrum_id.unique()))
    max_num_removed = len(set(structural_similarities[structural_similarities['Tanimoto_Similarity'] >= min_threshold].spectrumid1.unique()) | set(disallowed_spectra_daylight.spectrum_id.unique()))
    min_num_points = len(train_summary) - max_num_removed
    
    # print(f"Length of Train Summary before Removing Intersecting Scaffolds: {len(train_summary)}")
    # train_summary = train_summary.loc[~train_summary['scaffold'].isin(shared_scaffolds)]
    # print(f"Length of Train Summary after Removing Intersecting Scaffolds: {len(train_summary)}")
    
    for threshold in sorted(similarity_thresholds, reverse=True): # Sort in descending order to remove the most similar spectra first
        output_train_csv_path = input_train_csv[:-4] + f"_double_structural_{threshold}.csv"
        output_train_mgf_path = input_train_mgf[:-4] + f"_double_structural_{threshold}.mgf"
        
        to_drop = structural_similarities.loc[structural_similarities['Tanimoto_Similarity'] >= threshold].spectrumid1.values
        # Add the disallowed spectra from the daylight similarity
        disallowed_smiles_daylight = daylight_similarities[daylight_similarities >= threshold].index.values
        disallowed_spectra_daylight = train_summary.loc[train_summary['Smiles'].isin(disallowed_smiles_daylight)]
        to_drop = set(to_drop) | set(disallowed_spectra_daylight.spectrum_id.values)
        to_save = train_summary.loc[~train_summary['spectrum_id'].isin(to_drop)]
        
        
        # Subsample summary for saving
        to_save = to_save.sample(n=min_num_points, replace=False)
        to_save.to_csv(output_train_csv_path)
        synchronize_spectra(input_train_mgf, output_train_mgf_path, to_save)

def split_data_mcs(input_train_csv, input_train_mgf, input_test_csv, input_test_mgf, similarity_thresholds, progress_bar=True):
    """This function computes the maximum common substructure (MCS) between each pair of structures and removes all instances in the training set
    that have a MCS similarity above the specified thresholds. The final output is a set of csv and mgf files for the training and test sets.
    
    Truncate the training set to the same size?
    """
    train_summary = pd.read_csv(input_train_csv)
    train_summary = train_summary.loc[~ train_summary.Smiles.isna()]
    test_summary = pd.read_csv(input_test_csv)
    test_summary = test_summary.loc[~ test_summary.Smiles.isna()]
    
    output_test_csv_path  = input_test_csv[:-4] + f"_mcs.csv"
    output_test_mgf_path  = input_test_mgf[:-4] + f"_mcs.mgf"
    
    test_summary.to_csv(output_test_csv_path)
    synchronize_spectra(input_test_mgf, output_test_mgf_path, test_summary)
    
    def _generate_mcs(smiles1, smiles2):
        """Generate the jaccard index between two sets of molecular structures. From the SIMILE paper."""
        mol1 = Chem.MolFromSmiles(smiles1)
        mol2 = Chem.MolFromSmiles(smiles2)
        # Default parameters are ringMatchesRingOnly=False, and timeout=3600
        # SIMILE uses ringMatchesRingOnly=False, and timeout=1800
        res = rdFMCS.FindMCS([mol1, mol2], ringMatchesRingOnly=False, timeout=3600)
        ab = res.numBonds
        a  = mol1.GetNumBonds()
        b  = mol2.GetNumBonds()
        jaccard_similarity = ab / (a + b - ab)
        return jaccard_similarity
    
    train_smiles = train_summary.Smiles.unique()
    test_smiles = test_summary.Smiles.unique()
    
    # Calculate all pairs of MCS and store the results in a  dataframe
    sims = np.ones((len(train_smiles), len(test_smiles)))
    
    if progress_bar:
        p_bar = tqdm(train_smiles)
    else:
        p_bar = train_smiles
    
    for i, train_smile in enumerate(p_bar):
        for j in range(i+1, len(test_smiles)):
            test_smile = test_smiles[j]
            sim = _generate_mcs(train_smile, test_smile)
            sims[i,j] = sim
            sims[j,i] = sim
            
    sims = pd.DataFrame(sims, index=train_smiles, columns=test_smiles)  # Train is rows, test is columns
    
    # Get the maximum number of data poitns to be dropped in order to hold train size fixed.
    min_threshold = min(similarity_thresholds)
    max_num_removed = len(sims[sims.max(axis=1) >= min_threshold].sum())    # Count all train entries where the max train-test similarity exceeds the threshold
    min_num_points = len(train_summary) - max_num_removed
    
    for threshold in sorted(similarity_thresholds, reverse=True): # Sort in descending order to remove the most similar spectra first
        output_train_csv_path = input_train_csv[:-4] + f"_structural_{threshold}.csv"
        output_train_mgf_path = input_train_mgf[:-4] + f"_structural_{threshold}.mgf"
        
        to_drop = sims.loc[sims.max(axis=1) >= threshold]
        to_save = train_summary.loc[~train_summary['spectrum_id'].isin(to_drop)]
        
        
        # Subsample summary for saving
        temp = to_save.sample(n=min_num_points, replace=False)
        temp.to_csv(output_train_csv_path)
        synchronize_spectra(input_train_mgf, output_train_mgf_path, temp)


def split_data_structural_scaffold(input_train_csv, input_train_mgf, input_test_csv, input_test_mgf, structural_similarities, similarity_thresholds, progress_bar=True):
    # Column names are spectrumid1 spectrumid2 Tanimoto_Similarity
    
    def _generate_scaffold(smiles):
        smiles = rdMolStandardize.StandardizeSmiles(smiles)
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Could not generate mol from smiles: {smiles}")
        
        mol_scaffold_generic = MurckoScaffold.GetScaffoldForMol(mol)
        
        if mol_scaffold_generic is None:
            raise ValueError(f"Could not generate scaffold from smiles: {smiles}")
        
        try:
            return Chem.CanonSmiles(Chem.MolToSmiles(mol_scaffold_generic))
        except Exception as e:
            print(f"Failed to generate scaffold from smiles: {smiles}")
            return None
    
    # Make sure to remove structures with no smiles strings since we can't split them
    train_summary = pd.read_csv(input_train_csv)
    train_summary = train_summary.loc[~ train_summary.Smiles.isna()]

    test_summary = pd.read_csv(input_test_csv)
    test_summary = test_summary.loc[~ test_summary.Smiles.isna()]
    
    output_test_csv_path  = input_test_csv[:-4] + f"_structural_scaffold.csv"
    output_test_mgf_path  = input_test_mgf[:-4] + f"_structural_scaffold.mgf"
       
    structural_similarities = pd.read_csv(structural_similarities)
    
    # Remove obvious ones where spectrum_id is shared
    all_test_ids = set(test_summary.spectrum_id.values)
    train_summary = train_summary.loc[~train_summary['spectrum_id'].isin(all_test_ids)]
    
    # Remove identical scaffolds from the training set
    train_smiles = train_summary.Smiles.unique()
    train_scaffold_map = {smiles: _generate_scaffold(smiles) for smiles in train_smiles}
    
    test_smiles = test_summary.Smiles.unique()
    test_scaffold_map = {smiles: _generate_scaffold(smiles) for smiles in test_smiles}
    
    # Remove all failed scaffolds
    failed_smiles = [k for k, v in train_scaffold_map.items() if v is None] + [k for k, v in test_scaffold_map.items() if v is None]
    test_summary = test_summary.loc[~test_summary['Smiles'].isin(failed_smiles)]
    train_summary = train_summary.loc[~train_summary['Smiles'].isin(failed_smiles)]
    
    test_summary['scaffold'] = test_summary.Smiles.apply(lambda x: test_scaffold_map[x])
    train_summary['scaffold'] = train_summary.Smiles.apply(lambda x: train_scaffold_map[x])
    
    test_summary.to_csv(output_test_csv_path)
    synchronize_spectra(input_test_mgf, output_test_mgf_path, test_summary)
    
    # Remove shared scaffolds
    shared_scaffolds = set(test_summary.scaffold.values).intersection(set(train_summary.scaffold.values))
    print(f"Shared number of scaffolds: {len(shared_scaffolds)}")
    removed_from_train = train_summary.loc[train_summary['scaffold'].isin(shared_scaffolds)]
    print(f"Number of spectra with scaffolds in common with test {len(removed_from_train)}")
    
    # Get the maximum number of data poitns to be dropped in order to hold train size fixed.
    min_threshold = min(similarity_thresholds)
    max_num_removed = len(set(structural_similarities[structural_similarities['Tanimoto_Similarity'] >= min_threshold].spectrumid1.unique())|set(removed_from_train.spectrum_id.unique()))
    # max_num_removed = len(structural_similarities[structural_similarities['Tanimoto_Similarity'] >= min_threshold].spectrumid1.unique())
    min_num_points = len(train_summary) - max_num_removed
    
    print(f"Length of Train Summary before Removing Intersecting Scaffolds: {len(train_summary)}")
    train_summary = train_summary.loc[~train_summary['scaffold'].isin(shared_scaffolds)]
    print(f"Length of Train Summary after Removing Intersecting Scaffolds: {len(train_summary)}")
    
    for threshold in sorted(similarity_thresholds, reverse=True): # Sort in descending order to remove the most similar spectra first
        output_train_csv_path = input_train_csv[:-4] + f"_structural_scaffold_{threshold}.csv"
        output_train_mgf_path = input_train_mgf[:-4] + f"_structural_scaffold_{threshold}.mgf"
        
        to_drop = structural_similarities.loc[structural_similarities['Tanimoto_Similarity'] >= threshold].spectrumid1.values
        train_summary = train_summary.loc[~train_summary['spectrum_id'].isin(to_drop)]
        
        
        # Subsample summary for saving
        temp = train_summary.sample(n=min_num_points, replace=False)
        temp.to_csv(output_train_csv_path)
        synchronize_spectra(input_train_mgf, output_train_mgf_path, temp)
        
def split_data_scaffold(input_train_csv, input_train_mgf, input_test_csv, input_test_mgf, progress_bar=True):
    """This function removes all overlapping data points from the training set based on the general Murcko scaffold. 
    For mroe information on this approach, see here: 
    https://www.blopig.com/blog/2021/06/out-of-distribution-generalisation-and-scaffold-splitting-in-molecular-property-prediction/
    """
    train_summary = pd.read_csv(input_train_csv)
    train_summary = train_summary.loc[~ train_summary.Smiles.isna()]

    test_summary = pd.read_csv(input_test_csv)
    test_summary = test_summary.loc[~ test_summary.Smiles.isna()]
    
    output_test_csv_path  = input_test_csv[:-4] + f"_scaffold.csv"
    output_test_mgf_path  = input_test_mgf[:-4] + f"_scaffold.mgf"
    
    # Generate a dictionary from the smiles to the generalized murcko scaffold
    def _generate_scaffold(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Could not generate mol from smiles: {smiles}")
        
        mol_scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        mol_scaffold_generic = MurckoScaffold.MakeScaffoldGeneric(mol_scaffold)
        
        return Chem.CanonSmiles(Chem.MolToSmiles(mol_scaffold_generic))

    train_smiles = train_summary.Smiles.unique()
    train_scaffold_map = {smiles: _generate_scaffold(smiles) for smiles in train_smiles}
    
    test_smiles = test_summary.Smiles.unique()
    test_scaffold_map = {smiles: _generate_scaffold(smiles) for smiles in test_smiles}
    
    test_summary['scaffold'] = test_summary.Smiles.apply(lambda x: test_scaffold_map[x])
    train_summary['scaffold'] = train_summary.Smiles.apply(lambda x: train_scaffold_map[x])
    
    test_summary.to_csv(output_test_csv_path)
    synchronize_spectra(input_test_mgf, output_test_mgf_path, test_summary)
    
    # Remove shared scaffolds
    shared_scaffolds = set(test_summary.scaffold.values).intersection(set(train_summary.scaffold.values))
    train_summary = train_summary.loc[~train_summary['scaffold'].isin(shared_scaffolds)]
    
    output_train_csv_path = input_train_csv[:-4] + f"_scaffold.csv"
    output_train_mgf_path = input_train_mgf[:-4] + f"_scaffold.mgf"
    
    train_summary.to_csv(output_train_csv_path)
    synchronize_spectra(input_train_mgf, output_train_mgf_path, train_summary)

def split_data_murcko_histogram(input_train_csv, input_train_mgf, input_test_csv, input_test_mgf, similarity_thresholds=[3,4,5,6], progress_bar=True,):
    """This function removes all overlapping data points from the training set based on the Murcko Scaffold Histogram 
    introduced in "Emergence of molecular structures from self-supervised learning on mass spectra" by Bushuiev et al.
    """
    train_summary = pd.read_csv(input_train_csv)
    train_summary = train_summary.loc[~ train_summary.Smiles.isna()]
    test_summary = pd.read_csv(input_test_csv)
    test_summary = test_summary.loc[~ test_summary.Smiles.isna()]
    
    output_test_csv_path  = input_test_csv[:-4] + f"_murcko_histogram.csv"
    output_test_mgf_path  = input_test_mgf[:-4] + f"_murcko_histogram.mgf"
    
    # Calculate the Murcko Scaffold Histogram for the training and test sets
    train_smiles_hist_mapping = {smiles: calculate_murcko_histogram(smiles) for smiles in train_summary.Smiles.unique()}
    test_smiles_hist_mapping = {smiles: calculate_murcko_histogram(smiles) for smiles in test_summary.Smiles.unique()}
       
    # Generate the pairwise distance for each of the histograms using calculate_murcko_histogram_distance(histogram1, histogram2, min_rings=4, min_distance=5)
    
    distance_matrix = np.zeros((len(train_smiles_hist_mapping), len(test_smiles_hist_mapping)))
    train_smiles = train_smiles_hist_mapping.keys()
    test_smiles = test_smiles_hist_mapping.keys()
    
    p_bar = tqdm(train_smiles_hist_mapping.items()) if progress_bar else train_smiles_hist_mapping.items()
    
    for i, (smiles1, hist1) in enumerate(p_bar):
        for j, (smiles2, hist2) in enumerate(test_smiles_hist_mapping.items()):
            distance_matrix[i,j] = calculate_murcko_histogram_distance(hist1, hist2)

    distance_df = pd.DataFrame(distance_matrix, index=train_smiles, columns=test_smiles)
    # DEBUG
    distance_df.to_csv('murcko_histogram_distance.csv')
    
    # Remove data points from the training set that have overlapping Murcko Scaffold Histograms with the test set
    for threshold in similarity_thresholds:
        output_train_csv_path = input_train_csv[:-4] + f"_murcko_histogram_{threshold}.csv"
        output_train_mgf_path = input_train_mgf[:-4] + f"_murcko_histogram_{threshold}.mgf"
        
        to_drop = distance_df.loc[distance_df.min(axis=1) <= threshold]
        to_save = train_summary.loc[~train_summary['Smiles'].isin(to_drop.index)]
        
        to_save.to_csv(output_train_csv_path)
        synchronize_spectra(input_train_mgf, output_train_mgf_path, to_save)
    

def split_data(input_train_csv, input_train_mgf, input_test_csv, input_test_mgf, spectral_similarities=None, structural_similarities=None, similarity_thresholds=[1.0]):
    """This function removes all instances of scans listed in the similarity file from the input csv and mgf files.
    If a similarity threshold is provided, only similarities above the listed threshold will be removed.    
    """
    
    # Up here for debug
    # split_data_murcko_histogram(input_train_csv, input_train_mgf, input_test_csv, input_test_mgf, similarity_thresholds=[3,4,5,6], progress_bar=True)
    
    if spectral_similarities is not None:
        split_data_spectral(input_train_csv, input_train_mgf, input_test_csv, input_test_mgf, spectral_similarities, similarity_thresholds)
    if structural_similarities is not None:
        split_data_double_structural(input_train_csv, input_train_mgf, input_test_csv, input_test_mgf, structural_similarities, similarity_thresholds)
        split_data_structural(input_train_csv, input_train_mgf, input_test_csv, input_test_mgf, structural_similarities, similarity_thresholds)
        # split_data_mcs(input_train_csv, input_train_mgf, input_test_csv, input_test_mgf, similarity_thresholds)split_data_structural_scaffold
        # split_data_structural_scaffold(input_train_csv, input_train_mgf, input_test_csv, input_test_mgf, structural_similarities, similarity_thresholds)

def main():
    parser = argparse.ArgumentParser(description='Split Data')
    parser.add_argument('--input_train_csv', type=str, help='Input Train CSV')
    parser.add_argument('--input_train_mgf', type=str, help='Input Test MGF')
    parser.add_argument('--input_test_csv', type=str, help='Input Test CSV')
    parser.add_argument('--input_test_mgf', type=str, help='Input Test MGF')
    parser.add_argument('--spectral_similarities', help='Similarity file from fasst search containing indices to remove from training set.')
    parser.add_argument('--strucutral_similarities', help='Similarity file from fasst search containing indices to remove from training set.')
    parser.add_argument('--similarity_thresholds', help='Similarity threshold for removing spectra from training set.', nargs='+', type=float)
    args = parser.parse_args()
    split_data(args.input_train_csv, args.input_train_mgf, args.input_test_csv, args.input_test_mgf, args.spectral_similarities, args.strucutral_similarities,  args.similarity_thresholds)

if __name__ == "__main__":
    main()
    
def test_split_data_structural():
    # Create test data
    import os
    import pandas as pd
    import numpy as np
    
    if not os.path.exists("./test_data"):
        os.makedirs("./test_data")
        
    train_csv = "./test_data/train.csv"
    train_mgf = "./test_data/train.mgf"
    test_csv = "./test_data/test.csv"
    test_mgf = "./test_data/test.mgf"
    
    spectral_similarities = "./test_data/spectral_similarities.csv"
    structural_similarities = "./test_data/structural_similarities.csv"
    
    train_csv_df = pd.DataFrame({'spectrum_id': [f'spectrum_id_{x}' for x in np.arange(100)], 'scan': np.arange(100), 'Charge': np.ones(100), 'Smiles': ['C1CCCCC1' for x in np.arange(100)]})
    train_csv_df.to_csv(train_csv, index=False)
    test_csv_df = pd.DataFrame({'spectrum_id': [f'spectrum_id_{x}' for x in np.arange(50,100)], 'scan': np.arange(100, 150), 'Charge': np.ones(50), 'Smiles': ['C1CCCCC1' for x in np.arange(50)]})
    test_csv_df.to_csv(test_csv, index=False)
    
    with open(train_mgf, 'w') as f:
        for i in range(100):
            f.write(f"BEGIN IONS\nPEPMASS=100.0\nCHARGE=1+\nTITLE=spectrum_id_{i}\nSCANS={i}\n100.0 100.0\nEND IONS\n")
    with open(test_mgf, 'w') as f:
        for i in range(50,100):
            f.write(f"BEGIN IONS\nPEPMASS=100.0\nCHARGE=1+\nTITLE=spectrum_id_{i}\nSCANS={i}\n100.0 100.0\nEND IONS\n")
            
    structural_similarities_df = pd.DataFrame({'spectrumid1': [f'spectrum_id_{x}' for x in np.arange(50)], 'spectrumid2': [f'spectrum_id_{x}' for x in np.arange(50,100)], 'Tanimoto_Similarity': np.concatenate([np.ones(25), np.zeros(25)])})
    structural_similarities_df.to_csv(structural_similarities, index=False)
    
    split_data_structural(train_csv, train_mgf, test_csv, test_mgf, structural_similarities, [0.7], progress_bar=True)
    
    output_train_csv = "./test_data/train_structural_0.7.csv"
    output_train_mgf = "./test_data/train_structural_0.7.mgf"
    
    assert os.path.exists(output_train_csv)
    assert os.path.exists(output_train_mgf)
    
    output_train_csv_df = pd.read_csv(output_train_csv)
    output_train_mgf = IndexedMGF(output_train_mgf)
    
    # We start with 100, we lose 50 to ID collision, we lose 25 to structural similarity
    assert len(output_train_csv_df) == 25 
    assert len(output_train_mgf) == 25
    
    # Check the output ids are correct
    assert set(output_train_csv_df.spectrum_id.values) == set([f'spectrum_id_{x}' for x in np.arange(25, 50)])

    print(output_train_csv_df.spectrum_id.values)
    raise ValueError("Test Complete")