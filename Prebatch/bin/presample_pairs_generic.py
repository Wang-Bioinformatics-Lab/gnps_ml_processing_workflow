import argparse
import gc
import os
import logging
import pandas as pd
import numpy as np
import pickle
from tqdm import tqdm
from datetime import datetime
from typing import Iterator, List, NamedTuple, Optional
from data_generators import DataGeneratorAllInchikeys, FilteredPairsGenerator
import warnings

def main():
    parser = argparse.ArgumentParser(description='Train MS2DeepScore on the original data')
    parser.add_argument('--metadata', type=str, help='Path to the training data')
    parser.add_argument('--tanimoto_scores_path', type=str, help='Path to the tanimoto scores', default=None)
    parser.add_argument('--save_dir', type=str, help='Path to the save directory')
    parser.add_argument('--save_format', type=str, choices=['pickle', 'hdf5'], default='hdf5')
    parser.add_argument('--num_epochs', type=int, default=1800)
    parser.add_argument('--batch_size', type=int, default=32)
    parser.add_argument('--num_turns', type=int, default=2) # Default for MS2DeepScore Training is 2, Validation is 10
    parser.add_argument('--mode', type=str, choices=['standard', 'filter', 'triplet'], default='standard')
    parser.add_argument('--filter', action='store_true', help='Use memory efficient mode')
    parser.add_argument('--seed', type=int, default=42, help='Random seed for reproducibility')
    # Parameters for filtering, only applys in memory efficent mode
    parser.add_argument('--merge_on_lst', type=str, help='A semicolon delimited list of criteria to merge on. \
                                                            Options are: ["ms_mass_analyzer", "ms_ionisation", "adduct", "library"]. Only applys in memory efficent mode.',
                                                    default=None)
    parser.add_argument('--mass_analyzer_lst', type=str, default=None,
                                                help='A semicolon delimited list of allowed mass analyzers. All mass analyzers are \
                                                      allowed when not specified. Only applys in memory efficent mode.')
    parser.add_argument('--collision_energy_thresh', type=float, default=5.0,
                                                help='The maximum difference between collision energies of two spectra to be considered\
                                                    a pair. Default is <= 5. "-1.0" means collision energies are not filtered. \
                                                    If no collision enegy is available for either spectra, both will be included. Only applys in memory efficent mode.')
    parser.add_argument('--exponential_bins', help='Exponent for equal probability adjustment.', default=None, type=float)
    parser.add_argument('--num_bins', help='Number of bins for equal probability adjustment.', default=10, type=int)
    parser.add_argument('--strict_collision_energy', action='store_true', help="Require all pairs to have an associated collision energy.", default=False)
    args = parser.parse_args()

    if not args.mode =='filter':
        # Ensure None of the filtering flags were set
        if not args.merge_on_lst is None:
            raise ValueError("merge_on_lst is only compatible with --filter set.")
        if not args.mass_analyzer_lst is None:
            raise ValueError("mass_analyzer_lst is only compatible with --filter set.")
        if args.collision_energy_thresh != 5.0:
            raise ValueError("collision_energy_thresh is only compatible with --filter set.")

    current_time = datetime.now().strftime("%d_%m_%Y_%H_%M_%S")
    print(os.path.join(args.save_dir, 'logs', f"presampling_{current_time}.log"), flush=True)
    
    logging.basicConfig(format='[%(levelname)s]: %(message)s',
                    level=logging.DEBUG,
                    handlers=[logging.FileHandler(os.path.join(args.save_dir, 'logs', f"presampling_{current_time}.log"),
                                                mode='w'),
                            logging.StreamHandler()]
                    )
    
    logging.info('All arguments:')
    logging.info('Note: merge_on_lst, mass_analyzer_lst, collision_energy_thresh, no_pm_requirement are only used in memory efficent mode.')
    for arg in vars(args):
        logging.info('%s: %s', arg, getattr(args, arg))
    logging.info('')
    
    num_epochs = int(args.num_epochs)
    batch_size = int(args.batch_size)
    num_turns  = int(args.num_turns)

    if args.exponential_bins is not None:
        all_bins = np.linspace(0, 1.0, args.num_bins+1)**args.exponential_bins # Will spend 30% of it's time in the top 20% of the bins
        lower_bins = all_bins[:-1]
        upper_bins = all_bins[1:]
        same_prob_bins = list(zip(lower_bins, upper_bins))
    else:
        same_prob_bins = list(zip(np.linspace(0, 0.9, args.num_bins), np.linspace(0.1, 1, args.num_bins)))
    
    logging.info("Loading data...")
    metadata = pd.read_csv(args.metadata)
    metadata['inchikey_14'] = metadata['InChIKey_smiles'].str[:14]
    all_inchikeys = metadata['inchikey_14'].values
    unique_inchikeys = metadata['inchikey_14'].unique()
    spectrum_ids = metadata['spectrum_id'].unique()

    reference_scores_df = pd.read_csv(args.tanimoto_scores_path, index_col=0)
    print(reference_scores_df.iloc[:5, :5])
    reference_scores_df = reference_scores_df.loc[unique_inchikeys, unique_inchikeys]

    headers = None
    if args.mode == 'filter':
        logging.info("Creating FilteredPairsGenerator... (Computing Filtered Pairs on the Fly)")
        headers = ['spectrumid1', 'spectrumid2', 'inchikey1', 'inchikey2', 'score']
        training_generator = FilteredPairsGenerator(metadata,
                                                    reference_scores_df,
                                                    same_prob_bins=same_prob_bins,
                                                    shuffle=True,
                                                    random_seed=args.seed,
                                                    num_turns=num_turns,
                                                    batch_size=batch_size,
                                                    use_fixed_set=False,
                                                    ignore_equal_pairs=True,
                                                    merge_on_lst=args.merge_on_lst,
                                                    mass_analyzer_lst=args.mass_analyzer_lst,
                                                    collision_energy_thresh=args.collision_energy_thresh,
                                                    strict_collision_energy=args.strict_collision_energy)
    elif args.mode == 'triplet':
        logging.info("Creating DataGeneratorTriplets... (All-Pairs)")
        raise NotImplementedError()
        headers = ['anchor', 'positive', 'negative', 'anchor_inchikey', 'positive_inchikey', 'negative_inchikey', 'anchor_positive_score', 'anchor_negative_score']
        # training_gnerator = DataGeneratorTriplets(all_inchikeys,

    else:
        logging.info("Creating DataGeneratorAllInchikeys... (All-Pairs)")
        headers = ['spectrumid1', 'spectrumid2', 'inchikey1', 'inchikey2', 'score']
        training_generator = DataGeneratorAllInchikeys(all_inchikeys,
                                                        spectrum_ids,
                                                        reference_scores_df,
                                                        same_prob_bins=same_prob_bins,
                                                        shuffle=True,
                                                        random_seed=args.seed,
                                                        num_turns=num_turns,
                                                        batch_size=batch_size,
                                                        use_fixed_set=False,
                                                        ignore_equal_pairs=True)
    hdf_path = os.path.join(args.save_dir, 'data.hdf5')
    with pd.HDFStore(hdf_path, 'w') as store:
        for epoch_num in range(num_epochs):
            group_name = f'epoch_{epoch_num}'
            batch_accumulator = []
            for batch_num, batch in enumerate(tqdm(training_generator, desc=f"Epoch {epoch_num}")):
                # Save the batch DataFrame to the hdf5 file
                batch_accumulator.extend(list(batch))
            store.put(group_name, pd.DataFrame(batch_accumulator, columns=['spectrumid1', 'spectrumid2', 'inchikey1', 'inchikey2', 'score']), format='table')


if __name__ == "__main__":

    main()
