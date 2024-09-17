import os
import argparse
import pandas as pd
import matchms
import logging
from time import time
from tqdm import tqdm
import sys
from utils import synchronize_spectra, synchronize_spectra_to_json
import pyteomics.mgf as py_mgf


"""
This script performs basic filtering and common processing ofund in machine learning pipelines. 

1) All spectra in the GNPS-LIBRARY are removed from the dataset (to ensure maximal data quality).add()
2) The dataset is split based on the ion mode (positive or negative), sepcfied by the "ion_mode" argument.
3) The dataset is filtered based on the following criteria:
    a) The precursor m/z value is less than or equal to 1500.0 (by default).
    b) The precursor ppm error is less than or equal to 50 (by default).
    c) The dataset must have associated structures (by default).
4) The spectra are cleaned using standard MatchMS filters:
    a) Normalize intensities
    b) Select peaks by m/z (10-2000 by default)
    c) Filter peaks by relative intensity (0.001 by default)
    d) Require at least 5 peaks with at least 2.0 relative intensity (defalt values).
5) The filtered summary and spectra are written to the files "selected_summary.csv" and "selected_spectra.mgf" respectively.
"""

def missing_structure_check(summary:pd.DataFrame):
    """Check if there are missing structures in the summary DataFrame.

    Args:
        summary (pd.DataFrame): The DataFrame containing the summary information.

    Returns:
        bool: True if there are missing structures, False otherwise.
    """    
    if len(summary.loc[summary.Smiles.isna()]) > 0 or \
        len(summary.loc[summary.INCHI.isna()]) > 0 or \
        len(summary.loc[summary.InChIKey_smiles.isna()]) > 0:
        return True
    else:
        return False

def clean_spectra(summary:pd.DataFrame, spectra:list,
                  min_peaks=(5, 2.0),
                  mz_range=(10, 2000),
                  relative_intensity_filter=0.001):
    """This fuction applies standard MatchMS filters to clean the spectra. Spectra that do not pass cleaning
    are removed from the list. The function also removes the corresponding rows from the summary dataframe.
    The filters that are applied are:
    1. Normalize intensities
    2. Select by mz
    3. Filter by relative intensity
    4. Require at least min_peaks[0] peaks with at least min_peaks[1] relative intensity.
    

    Args:
        summary (pd.DataFrame): The summary dataframe with metadata for the spectra
        spectra (list): The list of spectra to be cleaned
        min_peaks (tuple, optional): The first value is the minimum number of peaks, the second value is the minimum relative intensity to be counted. Defaults to (5, 2.0).
        mz_range (tuple, optional): The range of peaks considered. Defaults to (10, 2000).
        relative_intensity_filter (float, optional): The minimum relative intensity to be included in the output spectrum. Defaults to 0.001.
    """    
    def _cleaning_helper(spectrum):
        logging_dict = {}
        
        # Normalize intensities
        if spectrum is not None:
            spectrum = matchms.filtering.normalize_intensities(spectrum)
        
        # Select by mz
        if spectrum is not None:
            original_size = len(spectrum.peaks)
            spectrum = matchms.filtering.select_by_mz(spectrum, mz_from=mz_range[0], mz_to=mz_range[1])
            peaks_removed = original_size - len(spectrum.peaks)
            logging_dict['peaks_removed'] = peaks_removed
            logging.debug("Removed %s peaks from spectrum %s", peaks_removed, spectrum.metadata.get('spectrum_id'))
        
        # Filter by relative intensity
        if spectrum is not None:
            original_size = len(spectrum.peaks)
            spectrum = matchms.filtering.select_by_relative_intensity(spectrum, intensity_from=relative_intensity_filter)
            peaks_removed = original_size - len(spectrum.peaks)
            logging.debug("Removed %s peaks from spectrum %s", peaks_removed, spectrum.metadata.get('spectrum_id'))
            logging_dict['peaks_removed_relative_intensity'] = peaks_removed
        
        # Require at least 10 peaks
        if spectrum is not None:
            existed_previously = spectrum is not None
            spectrum = matchms.filtering.require_minimum_of_high_peaks(spectrum, no_peaks=min_peaks[0], intensity_percent=min_peaks[1])
            removed_by_peak_count = False
            if spectrum is None and existed_previously:
                removed_by_peak_count = True
            logging_dict['removed_by_peak_count'] = removed_by_peak_count
            if removed_by_peak_count:
                logging.debug("Removed spectrum by peak count. Params: {#Peaks: %s, %%Int: %s}", str(min_peaks[0]), str(min_peaks[1]))
            
        return spectrum, logging_dict
    
    spectra = [_cleaning_helper(x) for x in tqdm(spectra)]
    
    logging_dicts = [x[1] for x in spectra]
    spectra = [x[0] for x in spectra]
   
    # Generate a single logging dict
    peaks_removed_list = []
    peaks_removed_int_list = []
    spectra_removed_by_peak_count = 0
    
    for log_dict in logging_dicts:
        peaks_removed_list.append(log_dict.get('peaks_removed', 0))
        peaks_removed_int_list.append(log_dict.get('peaks_removed_relative_intensity', 0))
        spectra_removed_by_peak_count += log_dict.get('removed_by_peak_count', 0)
        
    # Write info to logs
    logging.info("Initial Spectra Count spectra: %s", len(spectra))
    spectra = [x for x in spectra if x is not None]   # Remove spectra that failed cleaning
    logging.info("Total peaks removed due to m/z range: %s", sum(peaks_removed_list))
    logging.info("Total peaks removed by intensity: %s", sum(peaks_removed_int_list))
    logging.info("Total spectra removed by peak count: %s", spectra_removed_by_peak_count)
    logging.info("Final Number of Spectra: %s", len(spectra))
    
    allowed_ids = [x.metadata.get("spectrum_id") for x in spectra]
    summary = summary.loc[summary.spectrum_id.isin(allowed_ids)]
    
    return summary, spectra

def basic_spectral_filter(summary:pd.DataFrame, spectra:list,
                            max_precursor_ppm_error=50,
                            max_precursor_mz=1500.0,
                            require_structure=True):
    """
    Filter the summary and spectra based on specified criteria.

    Args:
        summary (pd.DataFrame): The summary dataframe containing information about the spectra.
        spectra (list): The list of spectra to be filtered.
        max_precursor_ppm_error (int, optional): The maximum allowed precursor ppm error. Defaults to 50.
        max_precursor_mz (float, optional): The maximum allowed precursor m/z value. Defaults to 1500.0.
        require_structure (bool, optional): Flag indicating whether the spectra must have associated structures. Defaults to True.

    Raises:
        ValueError: Raised if there are missing structures in the summary when require_structure is True.

    Returns:
        tuple: A tuple containing the filtered summary dataframe and the filtered spectra list.
    """
    # Apply filters that apply to summary
    time_start = time()
    org_len = len(summary)
    logging.info("Original Summary Count: %s", org_len)
    if require_structure:
        summary = summary.loc[~summary.Smiles.isna()]
        summary = summary.loc[~summary.INCHI.isna()]
        summary = summary.loc[~summary.InChIKey_smiles.isna()]
    logging.info("Filtered missing structures in %s seconds", time() - time_start, )
    logging.info("New Summary Count: %s Length difference: %s", len(summary), org_len - len(summary))

    time_start = time()
    org_len = len(summary)
    if max_precursor_mz is not None:
        summary = summary.loc[summary.Precursor_MZ <= max_precursor_mz]
    logging.info("Filtered by precursor_mz in %s seconds", time() - time_start)
    logging.info("New Summary Count: %s Length difference: %s", len(summary), org_len - len(summary))

    time_start = time()
    org_len = len(summary)
    if max_precursor_ppm_error is not None:
        if require_structure and missing_structure_check(summary):
            raise ValueError("""Found missing structures in summary, but max_precursor_ppm_error is not None.
                             Please remove missing structures or set require_structure=True.""")
        summary = summary.loc[summary.ppmBetweenExpAndThMass <= max_precursor_ppm_error]
    logging.info("Filtered by precursor_mz error in %s seconds", time() - time_start)
    logging.info("New Summary Count: %s Length difference: %s", len(summary), org_len - len(summary))

    allowed_ids = set(summary.loc[:, "spectrum_id"].values)

    time_start = time()
    org_len = len(spectra)
    spectra = [x for x in spectra if x.metadata.get("spectrum_id") in allowed_ids]
    logging.info("Removed filtered spectra from spectra list in %s seconds", time() - time_start)
    
    return summary, spectra

def select_data(input_csv_path:str, input_mgf_path:str, ion_mode:str):
    """
    Selects and processes data based on the given input CSV file path, input MGF file path, and ion mode.

    Args:
        input_csv_path (str): The file path of the input CSV file containing the data.
        input_mgf_path (str): The file path of the input MGF file containing the spectra.
        ion_mode (str): The ion mode to filter the data.

    Returns:
        None
    """
    # Load data
    logging.info("Loading data")
    summary = pd.read_csv(input_csv_path)

    # Load mgf files
    logging.info("Loading mgf files")
    spectra = list(matchms.importing.load_from_mgf(input_mgf_path))
    
    # Remove GNPS-LIBRARY from the train/test data
    logging.info("Removing GNPS-LIBRARY from the data")
    original_len = len(summary)
    summary = summary.loc[summary['GNPS_library_membership'] != 'GNPS-LIBRARY']
    logging.info("Original Summary Count: %s", original_len)
    logging.info("Number of Spectra not in GNPS-LIBRARY: %s", len(summary))

    # Remove BMDMS from the train/test data
    logging.info("Removing BMDMS from the data")
    original_len = len(summary)
    summary = summary.loc[summary['GNPS_library_membership'] != 'BMDMS-NP']
    summary = summary.loc[summary['GNPS_library_membership'] != 'MSMS-Pos-bmdms-np_20200811'] # If RIKEN is imported
    logging.info("Original Summary Count: %s", original_len)
    logging.info("Number of Spectra not in BMDMS: %s", len(summary))
    
    # Split based on ion mode
    logging.info('Selected Ion Mode: %s', ion_mode)
    original_len = len(summary)
    selected_summary = summary.loc[summary['Ion_Mode'] == ion_mode]
    logging.info("Original Summary Count: %s", original_len)
    logging.info("Number of Spectra With Correct Ion Mode: %s", len(selected_summary))
    
    # Filter spectra by metadata requirements
    logging.info("Filtering Spectra")
    start_time = time()
    orignal_count = len(selected_summary)

    selected_summary, spectra = basic_spectral_filter(selected_summary, spectra)
    logging.info("Filtered Spectra in %s seconds", time() - start_time)
    logging.info("Original Spectra Count: %s", orignal_count)
    logging.info("New Spectra Count: %s", len(selected_summary))
    
    logging.info("Cleaning Spectra")
    start_time = time()

    selected_summary, spectra = clean_spectra(selected_summary, spectra)
    logging.info("Cleaned Spectra in %s seconds", time() - start_time)
    
    ids_to_keep = set(selected_summary.loc[:, "spectrum_id"].values)
    logging.info("Found %s spectrum ids from the summary to keep", len(ids_to_keep))
    
    selected_spectra = [x for x in spectra if x.metadata.get("spectrum_id") in ids_to_keep]
    logging.info("Found %s spectra from spectra list to keep", len(selected_spectra))
    if len(selected_spectra) != len(selected_summary):
        raise ValueError("Spectra and Summary lengths do not match")
    
    # Write positive and negative spectra
    logging.info("Writing Summary")
    start_time = time()
    selected_summary.to_csv('selected_summary.csv', index=False)
    logging.info("Wrote selected_summary.csv in %s seconds", time() - start_time)
    
    logging.info("Writing Spectra")
    start_time = time()
    
    # matchms.exporting.save_as_mgf(selected_spectra, 'selected_spectra.mgf')
    
    # For some reason matchms.exporting.save_as_mgf will only output one spectra
    # This is a workaround for that code, but the above line should be used once working
    def output_gen():
        for spectrum in selected_spectra:
            spectrum_dict = {"m/z array": spectrum.peaks.mz,
                                "intensity array": spectrum.peaks.intensities,
                                "params": spectrum.metadata}
            if 'fingerprint' in spectrum_dict["params"]:
                del spectrum_dict["params"]["fingerprint"]
            yield spectrum_dict
    py_mgf.write(output_gen(), 'selected_spectra.mgf', file_mode="w")
    
    logging.info("Wrote selected_spectra.mgf in %s seconds", time() - start_time)

    logging.info("Writing Spectra to JSON")
    start_time = time()
    synchronize_spectra_to_json('selected_spectra.mgf', 'selected_spectra.json')
    logging.info("Wrote selected_spectra.json in %s seconds", time() - start_time)

def main():
    parser = argparse.ArgumentParser(description='Create Parallel Parameters.')
    parser.add_argument('--input_csv', type=str, help='Input CSV')
    parser.add_argument('--input_mgf', type=str, help='Input MGF')
    parser.add_argument('--ion_mode', type=str, choices=['positive','negative'], help='Ion Mode')
    args = parser.parse_args()
    
    logging.basicConfig(format='[%(levelname)s]: %(message)s', stream=sys.stdout)
    logging.getLogger().setLevel(logging.INFO)
    
    select_data(args.input_csv, args.input_mgf, args.ion_mode)

if __name__ == "__main__":
    main()