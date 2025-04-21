import argparse
import os
from pathlib import Path

from tqdm import tqdm
from joblib import Parallel, delayed
import tempfile

import numpy  as np
import pandas as pd
import dask.dataframe as dd
import time

import pickle

from data_generators import FilteredPairsGenerator

def build_pairs( main_inchikey,
                inchikey_list,
                tanimoto_df,
                summary,
                train_test_sim_dict,
                pairs_generator=None,
                temp_file_dir=None,
                buffer_size=7_900_000):
    # Buffer Size Calculation:
    # 80 gb (leaving 40 for overhead) / 32 / 254 bytes per row * 0.75 for extra safety

    hdf_store = None
    hdf_path = None
    try:
        # Tempfile for hdf storage
        hdf_path = tempfile.NamedTemporaryFile(dir=temp_file_dir, delete=False)

        # Calculate the similarity list in a memory-efficent manner
        hdf_store = pd.HDFStore(hdf_path.name)

        spectrum_id_list_i = summary.loc[main_inchikey, ['spectrum_id']].values
        if len(spectrum_id_list_i.shape) == 2:
            spectrum_id_list_i = spectrum_id_list_i.squeeze()
        
        max_main_train_test_sim  = train_test_sim_dict[main_inchikey]['max']
        min_main_train_test_sim  = train_test_sim_dict[main_inchikey]['min']
        mean_main_train_test_sim = train_test_sim_dict[main_inchikey]['mean']

        columns=['inchikey1', 'inchikey2','spectrumid1', 'spectrumid2', 'precursor_ppm_diff', 'ground_truth_similarity', 'inchikey1_max_test_sim', 'inchikey2_max_test_sim', 'mean_max_train_test_sim', 'mean_mean_train_test_sim',
                'max_max_train_test_sim', 'max_mean_train_test_sim', 'max_min_train_test_sim']

        
        spectrum_id_pm_dict = summary.loc[:, ['spectrum_id', 'Precursor_MZ']].set_index('spectrum_id').to_dict()['Precursor_MZ']

        output_list = []
        curr_buffer_size = buffer_size
        for inchikey_j in inchikey_list:

            max_j_train_test_sim  = train_test_sim_dict[inchikey_j]['max']
            min_j_train_test_sim  = train_test_sim_dict[inchikey_j]['min']
            mean_j_train_test_sim = train_test_sim_dict[inchikey_j]['mean']
            mean_mean_train_test_sim = np.mean([mean_main_train_test_sim, mean_j_train_test_sim])
            mean_max_train_test_sim = np.mean([max_main_train_test_sim, max_j_train_test_sim])
            max_max_train_test_sim = max(max_main_train_test_sim, max_j_train_test_sim)
            max_mean_train_test_sim = max(mean_main_train_test_sim, mean_j_train_test_sim)
            max_min_train_test_sim = max(min_main_train_test_sim, min_j_train_test_sim)
            gt = tanimoto_df.loc[main_inchikey, inchikey_j]

            spectrum_id_list_j = summary.loc[inchikey_j, ['spectrum_id']]
            if len(spectrum_id_list_j.shape) == 2:
                spectrum_id_list_j = spectrum_id_list_j.values.squeeze()            

            valid_pairs = None
            if pairs_generator is not None:
                # Filter by valid pairs only
                # Returns list of spectrum_id1, spectrumid2
                valid_pairs = pairs_generator.get_spectrum_pair_with_inchikey(main_inchikey, inchikey_j,
                                                                              return_type='ids',
                                                                              return_all=True)
                if valid_pairs == (None, None):
                    # No valid pairs, skip
                    continue

                # Structure as dataframe
                valid_pairs = pd.DataFrame(valid_pairs, columns=['spectrum_id_1', 'spectrum_id_2'])
                # Use a set since it's significantly faster than querying the index
                valid_spectra_1_ids = set(valid_pairs['spectrum_id_1'].values)
                valid_pairs = valid_pairs.set_index('spectrum_id_1')

            # For each pair of spectra in i, j, calculate the similarity
            for spectrum_id_i in spectrum_id_list_i:
                if valid_pairs is not None:
                    if spectrum_id_i not in valid_spectra_1_ids:
                        continue
                spectrum_i_pm = float(spectrum_id_pm_dict[spectrum_id_i])

                for spectrum_id_j in spectrum_id_list_j:
                    if valid_pairs is not None:
                        if spectrum_id_j not in valid_pairs.loc[spectrum_id_i, ['spectrum_id_2']].values:
                            continue
                    spectrum_j_pm = float(spectrum_id_pm_dict[spectrum_id_j])
                    # The only time this computation will cross below the main diagonal.
                    # (in terms of spectrum_ids) is when the inchikeys are the same.
                    # When this happens, we only want to compute the similarity once so
                    # these cases are not weighted differently
                    if main_inchikey == inchikey_j:
                        if spectrum_id_i < spectrum_id_j:
                            continue
                    
                    mz_diff = abs(spectrum_i_pm - spectrum_j_pm)
                    smaller_mz = min(spectrum_i_pm, spectrum_j_pm)
                    pm_ppm_diff = (mz_diff / smaller_mz) * 1e6

                    output_list.append((main_inchikey,
                                        inchikey_j,
                                        spectrum_id_i,
                                        spectrum_id_j,
                                        pm_ppm_diff,
                                        gt,
                                        max_main_train_test_sim,
                                        max_j_train_test_sim,
                                        mean_max_train_test_sim,
                                        mean_mean_train_test_sim,
                                        max_max_train_test_sim,
                                        max_mean_train_test_sim,
                                        max_min_train_test_sim,))
                    # print(f"Appending took {time() - start_time} seconds")

                    curr_buffer_size -= 1
                    if curr_buffer_size == 0:
                        # Store the similarity at /main_inchikey/inchikey_j
                        # Should have spectrumid1, spectrumid2, ground_truth_similarity, predicted_similarity

                        output_frame = pd.DataFrame(output_list, columns=columns)

                        hdf_store.put(f"{main_inchikey}", output_frame, format='table', append=True, track_times=False,
                                        min_itemsize={'inchikey1':14, 'inchikey2':14, 'spectrumid1': 50, 'spectrumid2': 50})
                        # print(f"Storing took {time() - start_time} seconds")
                        curr_buffer_size = buffer_size
                        output_list = []

        # Dump the remaining content
        output_frame = pd.DataFrame(output_list, columns=columns)

        hdf_store.put(f"{main_inchikey}", output_frame, format='table', append=True, track_times=False,
                        min_itemsize={'inchikey1':14, 'inchikey2':14, 'spectrumid1': 50, 'spectrumid2': 50})
        curr_buffer_size = buffer_size
        output_list = []

        hdf_store.close()
        hdf_path.close()
    except Exception as e:
        # Close the hdf store if needed and exists
        if hdf_store is not None and hdf_store.is_open:
            hdf_store.close()

        # Delete the hdf file if needed
        if os.path.exists(hdf_path.name):
            os.remove(hdf_path.name)

        raise e
    except KeyboardInterrupt as ki:
        # Close the hdf store if needed and exists
        if hdf_store is not None and hdf_store.is_open:
            hdf_store.close()

        # Delete the hdf file if needed
        if os.path.exists(hdf_path.name):
            os.remove(hdf_path.name)

        raise ki

    return hdf_path.name # Cleanup in serial process after concat

def build_test_pairs_list(  test_summary_path:str,
                            pairwise_similarities_path:str,
                            train_test_similarity_path:str,
                            output_path:str,
                            n_jobs:int=1,
                            merge_on_lst=None,
                            mass_analyzer_lst=None,
                            collision_energy_thresh=5.0,
                            filter=False,
                            temp_file_dir=None,
                            inital_index=0,
                            skip_merge=False,
                            strict_collision_energy=False):

    if not filter:
        if merge_on_lst is not None or mass_analyzer_lst is not None or collision_energy_thresh != 5.0:
            print("Warning: Filtering is not enabled, but merge_on_lst, mass_analyzer_lst, or collision_energy_thresh is set. Ignoring filtering parameters.")
    # Load summary
    print("Loading summary...", flush=True)
    test_summary = pd.read_csv(test_summary_path)
    if 'InChIKey_smiles_14' not in test_summary.columns:
        # Try to use InChIKey_smiles of InChIKey
        if 'InChIKey' in test_summary.columns:
            test_summary['InChIKey_smiles_14'] = test_summary['InChIKey'].str[:14]
        elif 'InChIKey_smiles' in test_summary.columns:
            test_summary['InChIKey_smiles_14'] = test_summary['InChIKey_smiles'].str[:14]
        else:
            raise ValueError("No InChIKey, InChIKey_smiles, or InChIKey_smiles_14 column found in test summary file.")

    # Index by inchikey_14
    test_summary = test_summary.set_index('InChIKey_smiles_14')

    # Load train_test similarities
    print("Loading train-test similarities...", flush=True)
    train_test_similarities = pd.read_csv(train_test_similarity_path, index_col=0)

    # Load pairwise similarities
    tanimoto_df = pd.read_csv(pairwise_similarities_path, index_col=0)

    unique_test_inchikeys = test_summary.index.unique()
    
    # Remove anything that doesn't have an inchikey in the pairwise similarities, it didn't make it through pre-processing
    unique_test_inchikeys = [inchikey for inchikey in unique_test_inchikeys if inchikey in tanimoto_df.index]
    test_summary = test_summary.loc[unique_test_inchikeys]

    # Remove bits we can with filtering
    if filter:
        local_mass_analyzer_lst = mass_analyzer_lst.split(';') if mass_analyzer_lst is not None else None
        if mass_analyzer_lst is not None:
            test_summary = test_summary[test_summary['msMassAnalyzer'].isin(local_mass_analyzer_lst)]
        local_merge_on_lst = merge_on_lst.split(';') if merge_on_lst is not None else []
        if merge_on_lst is None:
            local_merge_on_lst = ['ms_mass_analyzer', 'ms_ionisation', 'adduct']
        if 'ms_mass_analyzer' in local_merge_on_lst:
            test_summary = test_summary.dropna(subset=['msMassAnalyzer'])
        if 'ms_ionisation' in local_merge_on_lst:
            test_summary = test_summary.dropna(subset=['msIonisation'])
        if 'adduct' in local_merge_on_lst:
            test_summary = test_summary.dropna(subset=['Adduct'])

    # Remove the filtered parts + remove any inchikeys that have no spectra in the summary
    unique_test_inchikeys = test_summary.index.unique()
    tanimoto_df = tanimoto_df.loc[unique_test_inchikeys, unique_test_inchikeys]
    
    train_test_sim_dict = {}
    print("Precomputing train-test similarities...", flush=True)
    for inchikey in tqdm(unique_test_inchikeys):
        main_train_test_sim = train_test_similarities.loc[:, inchikey]
        train_test_sim_dict[inchikey] = {'max': main_train_test_sim.max(),
                                        'min': main_train_test_sim.min(),
                                        'mean': main_train_test_sim.mean()}

    # If we're filtering, create a pairs_generator
    if filter:
        # Rename col to expected
        test_summary['inchikey_14'] = test_summary.index
        test_summary['collisionEnergy'] = test_summary['collision_energy']
        pairs_generator = FilteredPairsGenerator(test_summary,
                                                tanimoto_df,
                                                ignore_equal_pairs=False,   # Will explicitly be handled in build_pairs
                                                merge_on_lst=merge_on_lst,
                                                mass_analyzer_lst=mass_analyzer_lst,
                                                collision_energy_thresh=collision_energy_thresh,
                                                strict_collision_energy=strict_collision_energy)
    else:
        pairs_generator = None

    # Compute parallel scores using _parallel_cosine
    output_hdf_files = []
    hdf_store = None
    temp_store = None
    try:
        print("Generating pairs of test data exhausively...", flush=True)
        output_hdf_files = Parallel(n_jobs=n_jobs)(delayed(build_pairs)(unique_test_inchikeys[i],
                                                                        unique_test_inchikeys[i:],
                                                                        tanimoto_df,
                                                                        test_summary,
                                                                        train_test_sim_dict,
                                                                        pairs_generator=pairs_generator,
                                                                        temp_file_dir=temp_file_dir,)
                                                                        for i in tqdm(range(inital_index, len(unique_test_inchikeys))))

        print(output_hdf_files)
        start_time = time.time()
        dask_df = dd.read_hdf(output_hdf_files, key='/*')
        print(f"Reading the hdf file lazily took {time.time() - start_time:.2f} seconds")
        start_time = time.time()
        # Will concatenate to single table, but use parallelism
        dask_df.to_parquet(output_path, overwrite=True, compute=True)
        print(f"Writing the hdf file took {time.time() - start_time:.2f} seconds")
    finally:
        if hdf_store is not None and hdf_store.is_open:
            hdf_store.close()
        
        if temp_store is not None and temp_store.is_open:
            temp_store.close()

        for hdf_file in output_hdf_files:
            if os.path.exists(hdf_file):
                os.remove(hdf_file)

def main():
    # Args
    parser = argparse.ArgumentParser()
    parser.add_argument("--n_jobs", type=int, default=1)
    parser.add_argument("--metadata_path", type=str, help="Path to the metadata csv file.")
    parser.add_argument("--pairwise_similarities_path", type=str, help="Path to the pairwise similarities csv file.")
    parser.add_argument("--train_test_similarities_path", help="Path to the train-test similarity csv file.")
    parser.add_argument("--output_path", type=str, help="Path to the output directory.")
    parser.add_argument("--temp_file_dir", type=str, help="Path to the temporary directory for storing intermediate files.", default=None)
    parser.add_argument("--inital_index", type=int, help="The initial index to start from.", default=0)
    parser.add_argument("--skip_merge", action='store_true', help="Whether to skip merging the output files.")
    parser.add_argument("--filter", action='store_true', help="Whether to filter the pairs.")
    # Filtering parameters
    parser.add_argument('--merge_on_lst', type=str, help='A semicolon delimited list of criteria to merge on. \
                                                            Options are: ["ms_mass_analyzer", "ms_ionisation", "adduct", "library"].\
                                                            Default behavior is to filter on ["ms_mass_analyzer", "ms_ionisation", "adduct"]. \
                                                            Ignored if --filter is not set.',
                                                    default=None)
    parser.add_argument('--mass_analyzer_lst', type=str, default=None,
                                                help='A semicolon delimited list of allowed mass analyzers. All mass analyzers are \
                                                      allowed when not specified. Ignored if --filter is not set.')
    parser.add_argument('--collision_energy_thresh', type=float, default=5.0,
                                                help='The maximum difference between collision energies of two spectra to be considered\
                                                    a pair. Default is <= 5. "-1.0" means collision energies are not filtered. \
                                                    If no collision enegy is available for either spectra, both will be included. \
                                                    Ignored if --filter is not set.')
    parser.add_argument('--strict_collision_energy', action='store_true', help='Whether to strictly filter by collision energy. \
                                                                                Unless set, pairs with NaN collision energy will be included.')
    args = parser.parse_args()

    if args.temp_file_dir == '':
        temp_file_dir = None
    else:
        temp_file_dir = args.temp_file_dir

    metadata_path = Path(args.metadata_path)
    assert metadata_path.exists(), f"Metadata file {metadata_path} does not exist."
    pairwise_similarities_path = Path(args.pairwise_similarities_path)    
    assert pairwise_similarities_path.exists(), f"Pairwise similarities file {pairwise_similarities_path} does not exist."
    train_test_similarities_path = Path(args.train_test_similarities_path)
    assert train_test_similarities_path.exists(), f"Train-test similarities file {train_test_similarities_path} does not exist."
    output_path = Path(args.output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Build the test pairs list
    build_test_pairs_list(metadata_path,
                          pairwise_similarities_path,
                          train_test_similarities_path,
                          output_path,
                          n_jobs=args.n_jobs,
                          merge_on_lst=args.merge_on_lst,
                          mass_analyzer_lst=args.mass_analyzer_lst,
                          collision_energy_thresh=args.collision_energy_thresh,
                          filter=args.filter,
                          temp_file_dir=temp_file_dir,
                          inital_index=args.inital_index,
                          skip_merge=args.skip_merge,
                          strict_collision_energy=args.strict_collision_energy)



if __name__ == "__main__":
    main()