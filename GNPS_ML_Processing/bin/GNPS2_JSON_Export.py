import argparse
import vaex
import pandas as pd
import json
from time import time
from tqdm import tqdm
import os

def generate_json_parquet(parquet_path:str, csv_path:str, output_path:str):
    """Generate a json file containing all the metadata in the specified csv and the spectrum in 
        parquet.

    Args:
        parquet_path (str): Path to the input parquet file
        csv_path (str): Path to the input csv file.
        output_path (str): Path to output the complete json file.
    """
    print("Opening Summary")
    start_time = time()
    summary = pd.read_csv(csv_path)
    print(f"Completed in {time() - start_time} seconds.")
    print("Opening Spectra")
    start_time = time()
    # spectra = vaex.open(parquet_path)
    spectra = pd.read_parquet(parquet_path)
    print(f"Completed in {time() - start_time} seconds.")
    
    print("Writing JSON")
    start_time = time()
    with open(output_path, 'w') as json_file:
        json_file.write("[")
        # Iterate over metadata table
        counter = 0
        for _, row in tqdm(summary.iterrows(), total=summary.shape[0]):
            row = row.to_dict()
                
            this_spectra = spectra[spectra.spectrum_id == row['spectrum_id']]
            if len(this_spectra) == 0:
                continue
            
            # A list of peaks: [(mz, i), (mz, i), ...]
            row['peaks_json'] = list(zip(this_spectra.mz.values.to_pylist(), this_spectra.i.values.to_pylist()))
            # Convert dict to dict of strings (really just adding double quotes)
            for k, v in row.items():
                row[k] = str(v)
                
            json.dump(row, json_file)
            # If this is the last entry, don't add the ","
            if counter != len(summary)-1:
                json_file.write(",")
            counter += 1
        json_file.write("]")
    print(f"Completed in {time() - start_time} seconds.")
        
def generate_json_mgf(mgf_path:str, csv_path:str, output_path:str, progress_bar=True, per_library=False):
    """Generates a gnps-style json file containing all of the spectra in the mgf and the corresponding metadata in the csv

    Args:
        mgf_path (str): A path to an mgf file.
        csv_path (str): A path to a csv file.
        output_path (str): Path to output the json file to.
        progress_bar (bool): Whether to show a progress bar. Defaults to True.
        per_library (bool): It True, exports additional json files for each GNPS library. Defaults to False.
    """
    from pyteomics.mgf import IndexedMGF
    
    mgf_file = IndexedMGF(mgf_path, index_by_scans=True)
    
    summary = pd.read_csv(csv_path)
    
    if output_path[-4:] == "json":
        if per_library:
            # A file output is provided, we cannot do per-library
            raise ValueError("Please provide an output path when per_libary is specified, rather than an output filename")
        # We have a named json output for our full library, no need to change
        main_output_path = output_path
    else: 
        # We assume a directory is provided:
        if not os.path.isdir(output_path):
            os.makedirs(output_path)
        main_output_path = os.path.join(output_path, "ALL_GNPS_cleaned.json")
    
    try:
        # Dict of files for all GNPS_libraries
        if per_library:
            all_libraries = summary.GNPS_library_membership.unique()
            per_library_jsons = {library_name: open(os.path.join(output_path, library_name+"_cleaned.json"), 'w') for library_name in all_libraries}
        
        print("Writing JSON")
        start_time = time()
        with open(main_output_path, 'w') as json_file:
            json_file.write("[")
            if per_library:
                for lib_json_file in per_library_jsons.values():
                    lib_json_file.write("[")
            # Iterate over mgf file
            counter = 0
            
            if progress_bar:
                mgf_file = tqdm(mgf_file)
            
            for _, spectra in enumerate(mgf_file):
                ccms_id = spectra['params']['title']
                row = summary.loc[summary['spectrum_id'] == ccms_id].iloc[0]
                if len(row) == 0:
                    continue
                row = row.to_dict()
                
                # A list of peaks: [(mz, i), (mz, i), ...]
                row['peaks_json'] = list(zip(spectra['m/z array'], spectra['intensity array']))
                # Convert dict to dict of strings (really just adding double quotes)
                for k, v in row.items():
                    row[k] = str(v)
                    
                json.dump(row, json_file)
                if per_library: 
                    json.dump(row, per_library_jsons[row['GNPS_library_membership']])
                # If this is the last entry, don't add the ","
                if counter != len(summary)-1:
                    json_file.write(",")
                counter += 1
            json_file.write("]")
            if per_library:
                for lib_json_file in per_library_jsons.values():
                    lib_json_file.write("]")
    finally:
        if per_library:
            for lib_json_file in per_library_jsons.values():
                lib_json_file.close()
    print(f"Completed in {time() - start_time} seconds.")

def main():
    parser = parser = argparse.ArgumentParser(description='Generate JSON Files.')
    parser.add_argument("--input_parquet_path", required=False)
    parser.add_argument("--input_mgf_path", required=False)
    parser.add_argument("--input_csv_path", required=True)
    parser.add_argument("--output_path", required=True)
    parser.add_argument("--progress", action="store_true")
    parser.add_argument("--per_library", action="store_true")
    args = parser.parse_args()
    
    if not args.input_parquet_path and not args.input_mgf_path:
        raise Exception("Must specify either input_parquet_path or input_mgf_path")
    if args.input_parquet_path and args.input_mgf_path:
        raise Exception("Must specify either input_parquet_path or input_mgf_path")
    if args.input_parquet_path and args.per_library:
        raise Exception("Per Library outputs are not supported for parquet files")
    
    if args.input_parquet_path:
        generate_json_parquet(args.input_parquet_path, args.input_csv_path, args.output_path)
    if args.input_mgf_path:
        generate_json_mgf(args.input_mgf_path, args.input_csv_path, args.output_path, args.progress, args.per_library)

if __name__ == "__main__":
    main()
    
def test_json_writing():
    temp_csv_path = './json_export_test.csv'
    temp_parquet_path = './json_export_test.parquet'
    temp_json_path = './json_export_test.json'
    
    try:
        import requests
        import json
        import numpy as np
        import os
        
        url = 'httpssacc://gnps.ucsd.edu/ProteoSAFe/SpectrumCommentServlet?SpectrumID=CCMSLIB00000071738'
        response = requests.get(url, timeout=10)

        assert response.status_code == 200

        summary = pd.DataFrame(json.loads(response.text)["annotations"][0], index=[0])
        spectra = pd.DataFrame(json.loads(response.text)["spectruminfo"], index=[0])
        
        # Rename library provenance column to mimic our naming:
        summary['spectrum_id'] = summary['SpectrumID']
        
        spectra_parquet = []
        level_0 = 0
        for spectrum_id in spectra.spectrum_id.unique():
            this_spectra = spectra[spectra.spectrum_id == spectrum_id]
            peaks = np.array(json.loads(this_spectra['peaks_json'].values.item()))
            mz_array = peaks[:,0]
            intensity_array = peaks[:,1]
            for index, (mz, intensity) in enumerate(zip(mz_array, intensity_array)):
                    spectra_parquet.append({'spectrum_id':spectrum_id, 'level_0': level_0, 'index':index, 'i':intensity, 'mz':mz})
                    level_0 += 1      
                    
        # Write temporary file to dataframe
        summary.to_csv(temp_csv_path)
        pd.DataFrame(spectra_parquet).to_parquet(temp_parquet_path)
        
        generate_json_mgf(temp_parquet_path, temp_csv_path, temp_json_path)
        
        
        # Begin actual testing:
        assert os.path.isfile(temp_json_path)
        
        with open('/home/user/SourceCode/GNPS_ML_Processing_Workflow/bin/json_export_test.json') as f:
            test_file = json.load(f)
            assert len(test_file) == 1
            assert len(test_file[0]['peaks_json']) > 10 # It's 11252, but this is just a sanity check
            assert test_file[0]['spectrum_id'] == 'CCMSLIB00000071738'
        
    finally:
        if os.path.isfile(temp_csv_path):
            os.remove(temp_csv_path)
        if os.path.isfile(temp_parquet_path):
            os.remove(temp_parquet_path)
        if os.path.isfile(temp_json_path):
            os.remove(temp_json_path)
        
    