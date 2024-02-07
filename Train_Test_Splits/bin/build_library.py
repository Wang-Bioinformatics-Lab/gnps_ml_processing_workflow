from pyteomics.mgf import IndexedMGF
import argparse
from pathlib import Path
import os
from tqdm import tqdm
import pandas as pd
import subprocess

def create_library(input_mgf, output_folder):
    build_input_folder = Path(output_folder).joinpath("input_library")
    build_input_folder.mkdir(exist_ok=True,parents=True)

    build_index_folder = Path(output_folder).joinpath("build_index")
    build_index_folder.mkdir(exist_ok=True,parents=True)

    file_level_metadata_path = Path(output_folder).joinpath("files.tsv")
    spectrum_level_metadata_path = Path(output_folder).joinpath("spectra.tsv")

    file_level_metadata = []
    spectrum_level_metadata = []

    input_mgf_basename = os.path.basename(input_mgf)
    temp_query_mgf = build_input_folder.joinpath("library_{}.mgf".format(os.path.splitext(input_mgf_basename)[0]))

    library_idx = 0 # Only one library

    file_level_metadata.append({
        'Idx':library_idx,
        'Filename':temp_query_mgf
    })

    library_spectra = IndexedMGF(input_mgf, index_by_scans=True)

    with open(temp_query_mgf, "w") as o:
        for i, spectrum in tqdm(enumerate(library_spectra)):
            # scan = i + 1  # This was previously what was used, but we don't know that scans are contiguous, so we'll use the spectrum scans
            scan = spectrum['params']['scans']
            
            spectrum_level_metadata.append({
                'FileIdx':library_idx,
                'Scan':scan,
                'GNPSLibraryAccession': spectrum['params']['title'],
            })

            # charge = int(spectrum['params']['charge'][0])
            charge =1

            o.write("BEGIN IONS\n")
            o.write("SEQ=*..*\n")
            o.write("SCANS={}\n".format(scan))
            o.write("CHARGE={}+\n".format(charge))
            o.write("MSLEVEL=2\n")
            o.write("PEPMASS={}\n".format(spectrum['params']['pepmass'][0]))
            for peak in zip(spectrum['m/z array'], spectrum['intensity array']):
                o.write("{} {}\n".format(peak[0], peak[1]))
            o.write("END IONS\n")


    df_file_level_metadata = pd.DataFrame(file_level_metadata)
    df_spectrum_level_metadata = pd.DataFrame(spectrum_level_metadata)

    df_file_level_metadata.to_csv(file_level_metadata_path, sep="\t", index=False)
    df_spectrum_level_metadata.to_csv(spectrum_level_metadata_path, sep="\t", index=False)
    
    # Return build_index_folder, spectrum_level_metadata_path, file_level_metadata_path to bash
    print(build_input_folder, build_index_folder, spectrum_level_metadata_path, file_level_metadata_path)
    
    # # Build Index
    # cmd = [
    #     "./bin/GNPS_FastSearch_Library/bin/main_execmodule", "ExecIndex", "./bin/generic_params",
    #     "-autoINDEX_TYPE", "mxc",
    #     "-autoINPUT_SPECS",  str(build_input_folder),
    #     "-autoOUTPUT_INDEX",  str(build_index_folder),
    #     "-autoSNR", "0",
    #     "-autoMIN_MZ_PEAK", "0",
    #     "-autoWINDOW_FILTER_RANK", "6",
    #     "-autoWINDOW_FILTER_SIZE", "50",
    #     "-autoMIN_PEAKS_AFTER_PROCESSING", "1"
    # ]
    
    # build_proc = subprocess.run(cmd)

    # # Add Metadata To Index
    # cmd = [
    #     "./bin/GNPS_FastSearch_Library/bin/main_execmodule", "ExecIndex", "./bin/generic_params",
    #     "-autoINPUT_INDEX",  str(build_index_folder),
    #     "-autoINPUT_METADATA_ANNOTATIONS", str(spectrum_level_metadata_path),
    #     "-autoINPUT_METADATA_FILENAMES", str(file_level_metadata_path)
    # ]
    # annotate_proc = subprocess.run(cmd)

def main():
    parser = argparse.ArgumentParser(description='Rewrite MGF file to fasst format.')
    parser.add_argument('--input_mgf', type=str, help='Input MGF')
    parser.add_argument('--output_folder', type=str, help='Folder to output library.')
    args = parser.parse_args()
    
    return create_library(args.input_mgf, args.output_folder)

if __name__ == "__main__":
    main()