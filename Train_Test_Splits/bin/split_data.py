import argparse
import pandas as pd
from pyteomics.mgf import IndexedMGF
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
            output_mgf.write("PEPMASS={}\n".format(float(spectra['params']['pepmass'][0])))
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
        
        spectral_similarities = spectral_similarities.loc[spectral_similarities['Cosine'] <= threshold]
        
        to_drop = spectral_similarities['Query_spectrumid'].values
        summary = summary.loc[~summary['spectrum_id'].isin(to_drop)]
        
        # Subsample summary for saving
        temp = summary.sample(n=min_num_points, replace=False)
        temp.to_csv(output_csv_path)
        synchronize_spectra(input_train_mgf, output_mgf_path, temp)
        

def split_data_structural(input_train_csv, input_train_mgf, input_test_csv, input_test_mgf, structural_similarities, similarity_thresholds, progress_bar=True):
    # Column names are spectrumid1 spectrumid2 Tanimoto_Similarity
    
    # Make sure to remove structures with no smiles strings since we can't split them
    train_summary = pd.read_csv(input_train_csv)
    train_summary = train_summary.loc[~ train_summary.Smiles.isna()]

    test_summary = pd.read_csv(input_test_csv)
    test_summary = test_summary.loc[~ test_summary.Smiles.isna()]
    
    output_test_csv_path  = input_test_csv[:-4] + f"_structural.csv"
    output_test_mgf_path  = input_test_mgf[:-4] + f"_structural.mgf"
    
    test_summary.to_csv(output_test_csv_path)
    synchronize_spectra(input_test_mgf, output_test_mgf_path, test_summary)
    
    structural_similarities = pd.read_csv(structural_similarities)
    
    # Get the maximum number of data poitns to be dropped in order to hold train size fixed.
    min_threshold = min(similarity_thresholds)
    max_num_removed = len(structural_similarities[structural_similarities['Tanimoto_Similarity'] >= min_threshold].spectrumid1.unique())
    min_num_points = len(train_summary) - max_num_removed
    
    for threshold in sorted(similarity_thresholds, reverse=True): # Sort in descending order to remove the most similar spectra first
        output_train_csv_path = input_train_csv[:-4] + f"_structural_{threshold}.csv"
        output_train_mgf_path = input_train_mgf[:-4] + f"_structural_{threshold}.mgf"

        structural_similarities = structural_similarities.loc[structural_similarities['Tanimoto_Similarity'] <= threshold]
        
        to_drop = structural_similarities['spectrumid1'].values
        train_summary = train_summary.loc[~train_summary['spectrum_id'].isin(to_drop)]
        
        
        # Subsample summary for saving
        temp = train_summary.sample(n=min_num_points, replace=False)
        temp.to_csv(output_train_csv_path)
        synchronize_spectra(input_train_mgf, output_train_mgf_path, temp)

def split_data(input_train_csv, input_train_mgf, input_test_csv, input_test_mgf, spectral_similarities=None, structural_similarities=None, similarity_thresholds=[1.0]):
    """This function removes all instances of scans listed in the similarity file from the input csv and mgf files.
    If a similarity threshold is provided, only similarities above the listed threshold will be removed.    
    """
    if spectral_similarities is not None:
        split_data_spectral(input_train_csv, input_train_mgf, input_test_csv, input_test_mgf, spectral_similarities, similarity_thresholds)
    if structural_similarities is not None:
        split_data_structural(input_train_csv, input_train_mgf, input_test_csv, input_test_mgf, structural_similarities, similarity_thresholds)

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