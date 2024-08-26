import argparse
import os
import sys
from glob import glob
import numpy as np
import pandas as pd
from pyteomics import mgf
from pyteomics.mgf import IndexedMGF
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.DataStructs import cDataStructs

from tqdm import tqdm
from time import time
from utils import missing_structure_check#, greedy_subset
from GNPS2_JSON_Export import synchronize_spectra_to_json
import matplotlib.pyplot as plt
import logging
from joblib import Parallel, delayed
import gc

np.random.seed(42)

def compute_similarity(fp1, fp2, i, j):
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    return i, j, similarity

def compute_chunk(chunk, fps):
    result_chunk = []
    for i, j in chunk:
        similarity = DataStructs.TanimotoSimilarity(fps[i], fps[j])
        result_chunk.append((i, j, similarity))
    return result_chunk

def parallel_compute_similarity_matrix(fps, unique_structures, n_jobs=8):
    # Pre-generate the pairs of indices for the upper triangle
    logging.debug("Computing Similarity Matrix")
    logging.debug("Pre-generating pairs")
    pairs = [(i, j) for i in tqdm(range(len(fps))) for j in range(i + 1, len(fps))]

    # Split pairs into chunks for parallel processing
    if n_jobs == -1:
        n_jobs = os.cpu_count()
    chunks = np.array_split(pairs, n_jobs*100)

    # Compute similarities in parallel
    logging.debug("Computing similarities")
    results = Parallel(n_jobs=n_jobs)(delayed(compute_chunk)(chunk, fps) for chunk in tqdm(chunks))
    logging.debug("Finished computing similarities. Reshaping into a dataframe")

    # Initialize the similarity matrix
    sim_matrix = np.zeros((len(fps), len(fps)))

    # Fill in the similarity matrix with computed values
    for chunk in results:
        for i, j, similarity in chunk:
            sim_matrix[i, j] = similarity
            sim_matrix[j, i] = similarity  # Mirror the value for the symmetric matrix

    # Set diagonal values to -1
    np.fill_diagonal(sim_matrix, -1)

    # Convert to DataFrame
    sim_matrix_df = pd.DataFrame(sim_matrix, index=unique_structures, columns=unique_structures)

    return sim_matrix_df


def sample_spectra(summary, num_test_points):
    return summary.sample(n=num_test_points)

def sample_structures(summary, num_test_points, coverage_threshold=2.6):
    if  missing_structure_check(summary):
        raise ValueError("Found missing structures in summary.")
    summary['InChIKey_smiles_14'] = summary['InChIKey_smiles'].str[:14]
    # counts = summary['InChIKey_smiles_14'].value_counts()
    # unique_test_inchis = set(greedy_subset(counts, num_test_points))
    
    # Randomly sample 14 values from the InChIKey_smiles_14 column
    unique_test_inchis = np.random.choice(summary['InChIKey_smiles_14'].unique(), num_test_points, replace=False)
    
    test_rows = summary.loc[summary['InChIKey_smiles_14'].isin(unique_test_inchis)]
    train_rows = summary.loc[~summary['InChIKey_smiles_14'].isin(unique_test_inchis)]
    unique_train_inchis = train_rows['InChIKey_smiles_14'].unique()
    unique_structures = summary['InChIKey_smiles_14'].unique()

    # Log all data
    logging.info("Number of test rows: %i", len(test_rows))
    logging.info("Number of test structures: %i", len(unique_test_inchis))
    logging.info("Number of training rows: %i", len(train_rows))
    logging.info("Number of training structures: %i", len(unique_train_inchis))

    logging.info("Computing Similarity Matrix")
    # sim_matrix = parallel_compute_similarity_matrix(fps, unique_structures)
    fps = [Chem.RDKFingerprint(Chem.MolFromSmiles(summary.loc[summary.InChIKey_smiles_14 == x, 'Smiles'].iloc[0])) for x in tqdm(unique_structures)]
    sim_matrix = np.zeros((len(fps), len(fps)))
    for i, fp1 in enumerate(tqdm(fps)):
        sim_matrix[i, i] = -1
        for j in range(i+1, len(fps)):
            fp2 = fps[j]
            sim_matrix[i, j] = DataStructs.TanimotoSimilarity(fp1, fp2)
            sim_matrix[j, i] = sim_matrix[i, j]
    sim_matrix = pd.DataFrame(sim_matrix, index=unique_structures, columns=unique_structures)
    
    logging.info("Calculating test pairs")
    test_pairs = get_pairs(unique_train_inchis, unique_test_inchis, test_rows, sim_matrix)
    logging.info("Number of test structure pairs: %i", len(test_pairs))

    # Create linear bins for the histogram
    linear_bins = np.linspace(0.2, 1, 21)

    # Compute the 2D histogram
    logging.info("Computing 2D histogram")

    # Create random_sampling_heatmaps directory if it doesn't exist
    if not os.path.exists('./random_sampling_heatmaps'):
        os.makedirs('./random_sampling_heatmaps', exist_ok=True)

    hist, xedges, yedges = np.histogram2d(test_pairs['train_test_similarity'], test_pairs['pairwise_similarity'], bins=[linear_bins, linear_bins])

    # Apply logarithmic transformation to the histogram values (avoid log(0) by adding a small constant)
    hist_log = np.log10(hist + 1)  # log10 is log base 10, adding 1 to handle zero values

    # Save unthresholded heatmap
    plt.figure(figsize=(10, 10))
    plt.title(f"Pairwise Similarity vs Train-Test Similarity \n \
                Random Sampling {num_test_points} structures")
    # Plot the pairwise similarity vs avg train-test similarity using imshow
    plt.imshow(hist_log.T, cmap='viridis', origin='lower', extent=[0.2, 1, 0.2, 1])

    plt.ylabel('Pairwise Similarity')
    plt.xlabel('Train-Test Similarity')

    plt.colorbar(label='Log Count')
    plt.savefig(f'./random_sampling_heatmaps/random_sampling_{num_test_points}_structures_unthresholded.png')

    # Save unthresholded heatmap to csv
    np.savetxt(f'./random_sampling_heatmaps/random_sampling_{num_test_points}_structures_unthresholded.csv', hist_log, delimiter=',')
    # Save no log heatmap to csv
    np.savetxt(f'./random_sampling_heatmaps/random_sampling_{num_test_points}_structures_unthresholded_no_log.csv', hist, delimiter=',')

    # Apply coverage threshold
    if coverage_threshold is not None:
        hist_log[hist_log < coverage_threshold] = np.nan

        # Get percent coverage
        percent_coverage = np.sum(~np.isnan(hist_log)) / hist_log.size
    else:
        percent_coverage = None

    # Plot the histogram
    logging.info("Plotting histogram")
    logging.getLogger().setLevel(logging.INFO)
    plt.figure(figsize=(10, 10))
    plt.title(f"Pairwise Similarity vs Train-Test Similarity \n \
                Random Sampling {num_test_points} structures Coverage: {coverage_threshold} \n % Coverage: {percent_coverage}")
    # Plot the pairwise similarity vs avg train-test similarity using imshow
    plt.imshow(hist_log.T, cmap='viridis', origin='lower', extent=[0.2, 1, 0.2, 1])

    plt.ylabel('Pairwise Similarity')
    plt.xlabel('Train-Test Similarity')

    plt.colorbar(label='Log Count')
    plt.savefig(f'./random_sampling_heatmaps/random_sampling_{num_test_points}_structures.png')

    
    return test_rows, train_rows

    
def random_tail_hybrid(summary, datapoints_per_bin:int=30, num_bins:int=20, pairwise_similarity_compensation=False):
    if missing_structure_check(summary):
        raise ValueError("Found missing structures in summary.")
    
    maximum_training_set_reduction = 0.80 # If we lose more than 0.xx(|train_set|), stop reassigning
    inital_margin = 1.8 # The inital margin for the count in each similarity bing
    
    # Check arg types 
    assert isinstance(datapoints_per_bin, int), "datapoints_per_bin should be an integer."
    assert isinstance(num_bins, int), "num_bins should be an integer."
   
    logging.info("Computing Test Set By Sampling the Tail of the Similarity Matrix and Iteratively Adding Points")
    logging.info("Computing Fingerprints")
    
    unique_structures = list(summary['Smiles'].unique())
    
    # Compute Fingerprints (RDK only)
    fps = [Chem.RDKFingerprint(Chem.MolFromSmiles(x)) for x in tqdm(unique_structures)]

    # Create similarity matrix
    logging.info("Computing Similarity Matrix")
    # sim_matrix = parallel_compute_similarity_matrix(fps, unique_structures)
    sim_matrix = np.zeros((len(fps), len(fps)))
    for i, fp1 in enumerate(tqdm(fps)):
        sim_matrix[i, i] = -1
        for j in range(i+1, len(fps)):
            fp2 = fps[j]
            sim_matrix[i, j] = DataStructs.TanimotoSimilarity(fp1, fp2)
            sim_matrix[j, i] = sim_matrix[i, j]
    sim_matrix = pd.DataFrame(sim_matrix, index=unique_structures, columns=unique_structures)

    """Generate the inital test set using the "freebies" or points 
    already have a maximum similarity to other points corresponding to our bins"""
    logging.info("Generating initial test set")
    # Bin all structures
    bins = pd.cut(sim_matrix.max(axis=1), np.linspace(0, 1, 21))

    # For bins with less than datapoints_per_bin structures, sample all
    sampled_structures = []
    for bin_, group in sim_matrix.groupby(bins):
        if len(group) > int(datapoints_per_bin * inital_margin):    # Add a bit of padding since these values will shift around once we start removing training structures
            sampled_structures.extend(group.sample(int(datapoints_per_bin * inital_margin)).index)
        else:
            sampled_structures.extend(group.index)
            
    test_set = sampled_structures
    training_set = list(set(unique_structures) - set(test_set))
    
    # To the test set, add some randomly sampled points
    training_set = list(set(unique_structures) - set(test_set))
    additional_test_points = np.random.choice(training_set, 250, replace=False)
    test_set = list(np.concatenate([np.unique(test_set), np.unique(additional_test_points)]))
    training_set = list(set(unique_structures) - set(test_set))
    
    
    """Iteratively mine out structures from the training set until each bin is at the desired size"""
    logging.info("Iteratively mining out structures")
    train_test_similarity = sim_matrix.loc[test_set, training_set]

    train_test_bins = pd.cut(train_test_similarity.max(axis=1), np.linspace(0, 1, 21))
    grouped_train_test_sim_matrix = train_test_similarity.max(axis=1).groupby(train_test_bins)
    smallest_bin = train_test_similarity.max(axis=1).groupby(train_test_bins).count().idxmin()

    inital_bins = grouped_train_test_sim_matrix.count()

    initial_train_size = len(training_set)

    total_num_structures_lost = 0
    
    train_test_similarity = sim_matrix.loc[test_set, training_set]
    
    train_test_bins = pd.cut(train_test_similarity.max(axis=1), np.linspace(bins[0], bins[1], bins[2]))
    grouped_train_test_sim_matrix = train_test_similarity.max(axis=1).groupby(train_test_bins)
    groups = sorted(list(grouped_train_test_sim_matrix.groups.keys()), key=lambda x: x.right)
    # current_bin_sizes = [len(grouped_train_test_sim_matrix.get_group(group)) for group in groups]
    current_bin_sizes = grouped_train_test_sim_matrix.size()
    bin_to_optimize = 0

    inital_bins = grouped_train_test_sim_matrix.count()

    # train_sims = sim_matrix.loc[training_set, training_set]
    initial_train_size = len(training_set)

    total_num_structures_lost = 0
    
    while min(current_bin_sizes) < datapoints_per_bin and len(training_set) > maximum_training_set_reduction * initial_train_size and bin_to_optimize < len(current_bin_sizes) - 1:
            logging.debug("Optimizing bin: %s", bin_to_optimize)
            logging.debug("Bin size: %i", current_bin_sizes.iloc[bin_to_optimize])
            
            # bins_over_threshold = [x - datapoints_per_bin for x in current_bin_sizes if x > datapoints_per_bin]
            # available_to_sample = sum(bins_over_threshold)
            bin_to_sample_from = -1
            for i in range(bin_to_optimize+1, len(current_bin_sizes)):
                if current_bin_sizes.iloc[i] != 0:
                    bin_to_sample_from = i
                    break
            if bin_to_sample_from == -1:
                logging.warning("No bins have any structures left.")
                break    
            
            available_to_sample = min(datapoints_per_bin - current_bin_sizes.iloc[bin_to_optimize], current_bin_sizes.iloc[bin_to_sample_from])
            logging.debug("%i structures available to sample from bin %s", available_to_sample, bin_to_sample_from)

            # Iteratively shift the highest similarity bins to the left until we have enough structures
            upper_bound = groups[bin_to_optimize].right
            # First we will try with just moving the bin to the right (later we may implement using the largest bin)
            to_move = grouped_train_test_sim_matrix.get_group(groups[bin_to_sample_from]).sample(available_to_sample).index.values

            logging.debug('Moving %i structures from %s to %s', len(to_move), bin_to_sample_from, bin_to_optimize)
            
            # "Move the structures"
            # Remove all of it's neighbors greater than the upper_bound from the trianing set
            structures_to_remove = ((train_test_similarity >= upper_bound).loc[to_move] >= upper_bound).any(axis='index')
            structures_to_remove = structures_to_remove[structures_to_remove].index.values
            training_set = list(set(training_set) - set(structures_to_remove))
            total_num_structures_lost += len(set(structures_to_remove))
            logging.debug("Removing %i structures", len(set(structures_to_remove)))

            # Reinitailize the similarity matrix
            train_test_similarity = sim_matrix.loc[test_set, training_set]
            train_test_bins = pd.cut(train_test_similarity.max(axis=1),  np.linspace(bins[0], bins[1], bins[2]))
            grouped_train_test_sim_matrix = train_test_similarity.max(axis=1).groupby(train_test_bins)
            groups = sorted(list(grouped_train_test_sim_matrix.groups.keys()), key=lambda x: x.right)
            logging.debug("Current Groups: %s", str(grouped_train_test_sim_matrix.groups.keys()))
            # current_bin_sizes = [len(grouped_train_test_sim_matrix.get_group(group)) if group in grouped_train_test_sim_matrix.groups.keys() else 0 for group in groups]
            current_bin_sizes = grouped_train_test_sim_matrix.size()


            print(grouped_train_test_sim_matrix.count())
            
            if len(grouped_train_test_sim_matrix.groups.get(groups[bin_to_optimize])) >= datapoints_per_bin:
                bin_to_optimize += 1

    final_bins = grouped_train_test_sim_matrix.count()
    
    test_rows = summary.loc[summary['Smiles'].isin(test_set)]
    train_rows = summary.loc[summary['Smiles'].isin(training_set)]
    logging.info("Number of test rows: %i", len(test_rows))
    logging.info("Number of test structures: %i", len(test_set))
    logging.info("Number of training rows: %i", len(train_rows))
    logging.info("Number of training structures: %i", len(training_set))
    logging.info("Total number of training structures removed: %i", total_num_structures_lost)
    
    if len(training_set) < maximum_training_set_reduction * initial_train_size:
        logging.warning("Stopped early due to losing more than %f of the inital training set.", maximum_training_set_reduction)
    
    # Check if logging level is debug
    if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
        logging.debug("Initial Bins:")
        logging.debug(inital_bins)
        logging.debug("Final Bins:")
        logging.debug(final_bins)
        
        # Change logging level back at last minute to supress matplotlib
        logging.getLogger().setLevel(logging.INFO)        
        
        # Save a plot of the bins
        plt.figure(figsize=(20, 5))
        plt.title("Initial vs Final Bins of Test Set")
        plt.bar(inital_bins.index.astype(str), inital_bins.values, alpha=0.5, label='Initial')
        plt.bar(final_bins.index.astype(str), final_bins.values, alpha=0.5, label='Final')
        plt.xlabel('Similarity Bin (Max Similarity to Train Set)')
        plt.ylabel('Number of Structures')
        plt.legend()
        plt.savefig('initial_vs_final_bins.png')
        
        # Sae a plot of the pairwise similarities
        plt.figure(figsize=(20, 5))
        plt.title("Pairwise Similarities of Test Set")
        test_sims = sim_matrix.loc[test_set, test_set]
        plt.hist(test_sims.values.flatten(), bins=100)
        plt.xlabel('Pairwise Similarity')
        plt.ylabel('Number of Pairs')
        plt.savefig('pairwise_similarities.png')
        
        logging.getLogger().setLevel(logging.DEBUG)

        # Plot train-test similarity vs pairwise similarity
        test_pairs = test_rows[['spectrum_id', 'InChIKey_smiles_14']]
        test_pairs = test_pairs.set_index('spectrum_id')
        # Cartesian product to form all pairs
        test_pairs = test_pairs.assign(key=1).merge(test_pairs.assign(key=1), on='key').drop('key', 1)
        # Get the similarity of the to spectra to the training set
        test_pairs['train_test_similarity'] = sim_matrix.loc[test_pairs['InChIKey_smiles_14_x'], test_pairs['InChIKey_smiles_14_y']].values
        # Average the similarity
        test_pairs['train_test_similarity'] = test_pairs['train_test_similarity'].apply(lambda x: np.mean(x))
        # Get the pairwise similarity
        test_pairs['pairwise_similarity'] = test_sims.loc[test_pairs['InChIKey_smiles_14_x'], test_pairs['InChIKey_smiles_14_y']].values

        plt.figure(figsize=(10,10))
        plt.title("Pairwise Similarity vs Train-Test Similarity")

        train_test_bins = pd.cut(train_test_similarity.max(axis=1),  np.linspace(bins[0], bins[1], bins[2]))
        grouped_train_test_sim_matrix = train_test_similarity.max(axis=1).groupby(train_test_bins)
        
        # Plot the the pairwise similarity vs avg train-test similarity using imshow
        plt.imshow(np.histogram2d(test_pairs['pairwise_similarity'], test_pairs['train_test_similarity'], bins=20)[0], cmap='viridis', origin='lower', extent=[0.2, 1, 0.2, 1])

        plt.xlabel('Pairwise Similarity')
        plt.ylabel('Train-Test Similarity')
        plt.savefig('pairwise_vs_train_test_similarity.png')
    
    return test_rows, train_rows

def get_pairs(train_inchis, test_inchis, test_rows, sim_matrix, ):
    # Plot train-test similarity vs pairwise similarity
    test_pairs = test_rows[['InChIKey_smiles_14']]

    # For plotting, remove duplicates
    test_pairs = test_pairs.drop_duplicates()

    train_test_similarity = sim_matrix.loc[test_inchis, train_inchis]

    max_train_test_sim = train_test_similarity.max(axis=1)

    gc.collect()

    # Cartesian product to form all pairs
    logging.info("Forming all pairs")
    test_pairs = test_pairs.assign(key=1).merge(test_pairs.assign(key=1), on='key').drop('key', axis=1)

    # Set inchikeys to categorical
    test_pairs['InChIKey_smiles_14_x'] = test_pairs['InChIKey_smiles_14_x'].astype('category')
    test_pairs['InChIKey_smiles_14_y'] = test_pairs['InChIKey_smiles_14_y'].astype('category')

    tqdm.pandas()

    # Get the similarity of the to spectra to the training set
    logging.info("Getting similarity")
    similarity_x = test_pairs['InChIKey_smiles_14_x'].map(max_train_test_sim).astype(float)
    similarity_y = test_pairs['InChIKey_smiles_14_y'].map(max_train_test_sim).astype(float)

    # test_pairs['train_test_similarity'] = test_pairs.apply(lambda x: [max_train_test_sim.at[x['InChIKey_smiles_14_x']], max_train_test_sim.at[x['InChIKey_smiles_14_y']]], axis=1)
    # Average the similarity
    logging.info("Averaging similarity")
    # test_pairs['train_test_similarity'] = test_pairs['train_test_similarity'].apply(lambda x: np.mean(x))
    test_pairs['train_test_similarity'] = (similarity_x + similarity_y) / 2
    del similarity_x, similarity_y
    gc.collect()

    # Get the pairwise similarity
    logging.info("Getting pairwise similarity")
    test_pairs['pairwise_similarity'] = test_pairs.progress_apply(lambda x: sim_matrix.at[x['InChIKey_smiles_14_x'], x['InChIKey_smiles_14_y']], axis=1)

    return test_pairs

def basic_sampling_scheme(summary,
                          n_lowest_sim_points = 3000, max_sim_threshold = 0.5, coverage_threshold = 2.6,
                          n_closest_neighbor=7):

    if  missing_structure_check(summary):
        raise ValueError("Found missing structures in summary.")
    summary['InChIKey_smiles_14'] = summary['InChIKey_smiles'].str[:14]

    # Now that we have our "test set" we want remove training points to force the max train-test similarity
    unique_structures = list(summary['InChIKey_smiles_14'].unique())

    logging.info("Computing Similarity Matrix")
    # sim_matrix = parallel_compute_similarity_matrix(fps, unique_structures)
    fps = [Chem.RDKFingerprint(Chem.MolFromSmiles(summary.loc[summary.InChIKey_smiles_14 == x, 'Smiles'].iloc[0])) for x in tqdm(unique_structures)]
    sim_matrix = np.zeros((len(fps), len(fps)))
    for i, fp1 in enumerate(tqdm(fps)):
        sim_matrix[i, i] = -1
        for j in range(i+1, len(fps)):
            fp2 = fps[j]
            sim_matrix[i, j] = DataStructs.TanimotoSimilarity(fp1, fp2)
            sim_matrix[j, i] = sim_matrix[i, j]

    sim_matrix = pd.DataFrame(sim_matrix, index=unique_structures, columns=unique_structures)

    truncated_sim = sim_matrix.copy()
    truncated_sim[truncated_sim < 0.2] = 0

    logging.info("Computing average and max similarity")
    nan_truncated_sim = truncated_sim.copy()
    nan_truncated_sim[nan_truncated_sim == 0] = np.nan
    average_similarity = nan_truncated_sim.mean(axis=1)
    max_similarity = nan_truncated_sim.max(axis=1)
    
    # Get points with max similarity greater than 0.5, and lowest average similarity possible
    logging.info("Getting test set")

    test_inchis = average_similarity[(max_similarity > max_sim_threshold)].nsmallest(n_lowest_sim_points).index

    # Sample their highest similarity neighbor
    logging.info("Sampling highest similarity neighbor")
    test_inchis_w_highest_sim_neighbors = []
    for inchi in test_inchis:
        test_inchis_w_highest_sim_neighbors.append(truncated_sim.loc[inchi].nlargest(n_closest_neighbor).index.tolist())

    # Flatten
    test_inchis_w_highest_sim_neighbors = [item for sublist in test_inchis_w_highest_sim_neighbors for item in sublist]

    test_inchis_w_highest_sim_neighbors = list(set(test_inchis_w_highest_sim_neighbors) | set(test_inchis))

    test_inchis = list(set(test_inchis_w_highest_sim_neighbors) | set(test_inchis))

    del test_inchis_w_highest_sim_neighbors
    gc.collect()

    # Get the test set
    test_rows = summary.loc[summary.InChIKey_smiles_14.isin(test_inchis)]
    unique_test_inchis = test_rows.InChIKey_smiles_14.unique()

    # Get the training set
    train_rows = summary.loc[~summary.InChIKey_smiles_14.isin(test_inchis)]
    unique_train_inchis = train_rows.InChIKey_smiles_14.unique()

    # Log all data
    logging.info("Number of test rows: %i", len(test_rows))
    logging.info("Number of test structures: %i", len(unique_test_inchis))
    logging.info("Number of training rows: %i", len(train_rows))
    logging.info("Number of training structures: %i", len(unique_train_inchis))
    
    logging.info("Calculating test pairs")
    test_pairs = get_pairs(unique_train_inchis, unique_test_inchis, test_rows, sim_matrix)
    logging.info("Number of test structure pairs: %i", len(test_pairs))

    # Create linear bins for the histogram
    linear_bins = np.linspace(0.2, 1, 21)

    # Compute the 2D histogram
    logging.info("Computing 2D histogram")

    # Create basic_sampling_heatmaps directory if it doesn't exist
    if not os.path.exists('./basic_sampling_heatmaps'):
        os.makedirs('./basic_sampling_heatmaps', exist_ok=True)

    hist, xedges, yedges = np.histogram2d(test_pairs['train_test_similarity'], test_pairs['pairwise_similarity'], bins=[linear_bins, linear_bins])

    # Apply logarithmic transformation to the histogram values (avoid log(0) by adding a small constant)
    hist_log = np.log10(hist + 1)  # log10 is log base 10, adding 1 to handle zero values

    # Save unthresholded heatmap
    plt.figure(figsize=(10, 10))
    plt.title(f"Pairwise Similarity vs Train-Test Similarity \n \
                n_lowest_sim_points: {n_lowest_sim_points}, Max Sim: {max_sim_threshold}, Coverage: {coverage_threshold} Closest: {n_closest_neighbor}")
    # Plot the pairwise similarity vs avg train-test similarity using imshow
    plt.imshow(hist_log.T, cmap='viridis', origin='lower', extent=[0.2, 1, 0.2, 1])

    plt.ylabel('Pairwise Similarity')
    plt.xlabel('Train-Test Similarity')

    plt.colorbar(label='Log Count')
    plt.savefig(f'./basic_sampling_heatmaps/n_lowest_{n_lowest_sim_points}_max_{max_sim_threshold}_cov_thresh_{coverage_threshold}_closest_{n_closest_neighbor}_unthresholded.png')

    # Save unthresholded heatmap to csv
    np.savetxt(f'./basic_sampling_heatmaps/n_lowest_{n_lowest_sim_points}_max_{max_sim_threshold}_cov_thresh_{coverage_threshold}_closest_{n_closest_neighbor}_unthresholded.csv', hist_log, delimiter=',')
    # Save no log heatmap to csv
    np.savetxt(f'./basic_sampling_heatmaps/n_lowest_{n_lowest_sim_points}_max_{max_sim_threshold}_cov_thresh_{coverage_threshold}_closest_{n_closest_neighbor}_unthresholded_no_log.csv', hist, delimiter=',')

    # Apply coverage threshold
    if coverage_threshold is not None:
        hist_log[hist_log < coverage_threshold] = np.nan

        # Get percent coverage
        percent_coverage = np.sum(~np.isnan(hist_log)) / hist_log.size
    else:
        percent_coverage = None

    # Plot the histogram
    logging.info("Plotting histogram")
    logging.getLogger().setLevel(logging.INFO)
    plt.figure(figsize=(10, 10))
    plt.title(f"Pairwise Similarity vs Train-Test Similarity \n \
                n_lowest_sim_points: {n_lowest_sim_points}, Max Sim: {max_sim_threshold}, Coverage: {coverage_threshold} \n % Coverage: {percent_coverage} Closest: {n_closest_neighbor}")
    # Plot the pairwise similarity vs avg train-test similarity using imshow
    plt.imshow(hist_log.T, cmap='viridis', origin='lower', extent=[0.2, 1, 0.2, 1])

    plt.ylabel('Pairwise Similarity')
    plt.xlabel('Train-Test Similarity')

    plt.colorbar(label='Log Count')
    plt.savefig(f'./basic_sampling_heatmaps/n_lowest_{n_lowest_sim_points}_max_{max_sim_threshold}_cov_thresh_{coverage_threshold}_closest_{n_closest_neighbor}.png')
    # plt.show()

    return test_rows, train_rows

def sample_structures_smart_inchikey(summary,
                                    similarity_matrix,
                                    num_test_roots=500,
                                    test_path_len=3,
                                    test_edge_cutoff=0.70,
                                    bins=(0.2, 1.0, 17),
                                    maximum_training_set_reduction=0.80, 
                                    datapoints_per_bin=143,
                                    move_structures=False):
    """
    Randomly sample num_test_roots structures and then iteratively remove structures from the training set
    until each bin has at least 30 structures.
    
    Parameters:
    summary: pd.DataFrame
        The summary dataframe
    similarity_matrix: pd.DataFrame
        The similarity matrix of the structures with index and columns as the InChIKey_smiles_14
    num_test_roots: int
        The number of structures to randomly sample. Used as the root for the random walks
    test_path_len: int
        The length of the random walk from each root
    test_edge_cutoff: float 
        The minimum similarity to be considered a neighbor in the random walk
    bins: tuple
        The bins to use for the similarity matrix
    maximum_training_set_reduction: float
        The maximum fraction of the training set that can be removed
    datapoints_per_bin: int
        The number of datapoints to have in each bin
    move_structures: bool
        Whether to move structures from the training set to the test set or remove them entirely.
    """
    
    if  missing_structure_check(summary):
        raise ValueError("Found missing structures in summary.")
    summary['InChIKey_smiles_14'] = summary['InChIKey_smiles'].str[:14]
    
    # Now that we have our "test set" we want remove training points to force the max train-test similarity
    # Get unique structures by inchikey14, get corresponding smiles
    inchikey_smiles_mapping = summary.groupby('InChIKey_smiles_14')['Smiles'].agg(lambda x: pd.Series.mode(x).iloc[0])
    unique_inchikeys = inchikey_smiles_mapping.index
    unique_smiles = inchikey_smiles_mapping.values

    logging.info("Reading Precomputed Similarity Matrix")
    sim_matrix = pd.read_csv(similarity_matrix, index_col=0)

    if not all([x in sim_matrix.index.values for x in unique_inchikeys]):
        raise ValueError("Not all unique inchikeys are in the similarity matrix index.")
    if not all([x in sim_matrix.columns.values for x in unique_inchikeys]):
        raise ValueError("Not all unique inchikeys are in the similarity matrix columns.")

    logging.info("Computing Test Set By Sampling the Tail of the Similarity Matrix")
    
    # For bins with less than datapoints_per_bin structures, sample all
    test_set = []
    training_set = unique_inchikeys
    train_test_similarity = sim_matrix.loc[test_set, training_set]
    train_test_bins = pd.cut(train_test_similarity.max(axis=1), np.linspace(bins[0], bins[1], bins[2]))
    tail_structures = []
    sim_matrix_bins = pd.cut(sim_matrix.max(axis=1), np.linspace(bins[0], bins[1], bins[2]))
    for bin_, group in sim_matrix.groupby(sim_matrix_bins):
        if len(group) > int(datapoints_per_bin):
            tail_structures.extend(group.sample(int(datapoints_per_bin)).index)
        else:
            tail_structures.extend(group.index)

    logging.info("Sampled %i tail structures", len(tail_structures))
    
    # Random walk matrix
    random_walk_matrix = sim_matrix.copy(deep=True)
    random_walk_matrix[random_walk_matrix < test_edge_cutoff] = 0
    random_walk_matrix[random_walk_matrix >= test_edge_cutoff] = 1
    
    # Randomly sample values from the InChIKey_smiles_14 column
    logging.debug("Sampling %i random structures to start the test set.", num_test_roots)
    structure_options = list(set(summary['InChIKey_smiles_14'].unique()) - set(tail_structures))
    rw_roots = np.random.choice(structure_options, num_test_roots, replace=False)           # structure_options not defined

    random_walk_structures = []
    logging.debug("Performing random walks of length %i from the test set roots.", test_path_len)
    for root in rw_roots:
        # Random Walk
        current = root
        for i in range(test_path_len):
            neighbors = random_walk_matrix.loc[current, :]
            neighbors = neighbors[neighbors == 1]
            if len(neighbors) == 0:
                break
            current = np.random.choice(neighbors.index)
            random_walk_structures.append(current)

    # Concatenate random walk roots and samples nodes
    random_walk_structures = np.concatenate([rw_roots, random_walk_structures])
            
    logging.info("Sampled %i nodes from random walks.", len(random_walk_structures))
              
    test_set = np.unique(np.concatenate([random_walk_structures, tail_structures]))
    
    # Get training structures by removing the test structures
    training_set = summary.loc[~summary['InChIKey_smiles_14'].isin(test_set), 'InChIKey_smiles_14'].unique()
        
    train_test_similarity = sim_matrix.loc[test_set, training_set]
    
    train_test_bins = pd.cut(train_test_similarity.max(axis=1), np.linspace(bins[0], bins[1], bins[2]))
    grouped_train_test_sim_matrix = train_test_similarity.max(axis=1).groupby(train_test_bins)
    groups = sorted(list(grouped_train_test_sim_matrix.groups.keys()), key=lambda x: x.right)
    current_bin_sizes = grouped_train_test_sim_matrix.size()
    # DEBUG
    print("Current Bin Sizes:")
    print(current_bin_sizes)
    bin_to_optimize = 0

    inital_bins = grouped_train_test_sim_matrix.count()

    initial_train_size = len(training_set)

    total_num_structures_lost = 0
    
    while min(current_bin_sizes) < datapoints_per_bin and len(training_set) > maximum_training_set_reduction * initial_train_size:
            if bin_to_optimize == len(current_bin_sizes) - 1:
                break
                # This needs more thinking, there may be no structures left to move to the bin
                bin_to_optimize = 0
            
            logging.debug("Optimizing bin: %s", bin_to_optimize)
            logging.debug("Bin size: %i", current_bin_sizes.iloc[bin_to_optimize])
            
            # bins_over_threshold = [x - datapoints_per_bin for x in current_bin_sizes if x > datapoints_per_bin]
            # available_to_sample = sum(bins_over_threshold)
            bin_to_sample_from = -1
            for i in range(bin_to_optimize+1, len(current_bin_sizes)):
                if current_bin_sizes.iloc[i] != 0:
                    bin_to_sample_from = i
                    break
            if bin_to_sample_from == -1:
                logging.warning("No bins have any structures left.")
                break    
            
            available_to_sample = min(datapoints_per_bin - current_bin_sizes.iloc[bin_to_optimize], current_bin_sizes.iloc[bin_to_sample_from])
            logging.debug("%i structures available to sample from bin %s", available_to_sample, bin_to_sample_from)

            # Iteratively shift the highest similarity bins to the left until we have enough structures
            upper_bound = groups[bin_to_optimize].right
            # First we will try with just moving the bin to the right (later we may implement using the largest bin)
            to_move = grouped_train_test_sim_matrix.get_group(groups[bin_to_sample_from]).sample(available_to_sample).index.values

            logging.debug('Moving %i structures from %s to %s', len(to_move), bin_to_sample_from, bin_to_optimize)
            
            # "Move the structures"
            # Remove all of it's neighbors greater than the upper_bound from the trianing set
            structures_to_remove = ((train_test_similarity >= upper_bound).loc[to_move] >= upper_bound).any(axis='index')
            structures_to_remove = structures_to_remove[structures_to_remove].index.values
            if not move_structures:
                training_set = list(set(training_set) - set(structures_to_remove))
                total_num_structures_lost += len(set(structures_to_remove))
                logging.debug("Removing %i structures", len(set(structures_to_remove)))
            else:
                test_set.extend(structures_to_remove)
                training_set = list(set(training_set) - set(structures_to_remove))
                logging.debug("Moving %i structures", len(set(structures_to_remove)))
            
            # Removing structures that go outside of our range, to help reduce confusion later
            structures_to_remove = train_test_similarity.loc[train_test_similarity.max(axis=1) < bins[0]].index.values
            test_set = list(set(test_set) - set(structures_to_remove))
            logging.debug("Removing %i structures from the test set that are below the minimum similarity threshold", len(structures_to_remove))

            # Reinitailize the similarity matrix
            train_test_similarity = sim_matrix.loc[test_set, training_set]
            train_test_bins = pd.cut(train_test_similarity.max(axis=1),  np.linspace(bins[0], bins[1], bins[2]))
            grouped_train_test_sim_matrix = train_test_similarity.max(axis=1).groupby(train_test_bins)
            groups = sorted(list(grouped_train_test_sim_matrix.groups.keys()), key=lambda x: x.right)
            logging.debug("Current Groups: %s", str(grouped_train_test_sim_matrix.groups.keys()))
            # current_bin_sizes = [len(grouped_train_test_sim_matrix.get_group(group)) if group in grouped_train_test_sim_matrix.groups.keys() else 0 for group in groups]
            current_bin_sizes = grouped_train_test_sim_matrix.size()


            print(grouped_train_test_sim_matrix.count())
            
            if len(grouped_train_test_sim_matrix.groups.get(groups[bin_to_optimize])) >= datapoints_per_bin:
                bin_to_optimize += 1

    final_bins = grouped_train_test_sim_matrix.count()
    
    test_rows = summary.loc[summary['InChIKey_smiles_14'].isin(test_set)]
    # TEMP: As a temporary measure, remove the BMDMS-NP from the test set
    # The wide range of collision energies is causing issues
    test_rows = test_rows.loc[test_rows['GNPS_library_membership'] != 'BMDMS-NP']
    
    
    train_rows = summary.loc[summary['InChIKey_smiles_14'].isin(training_set)]
    logging.info("Number of test rows: %i", len(test_rows))
    logging.info("Number of test structures: %i", len(test_set))
    logging.info("Number of training rows: %i", len(train_rows))
    logging.info("Number of training structures: %i", len(training_set))
    logging.info("Total number of training structures removed: %i", total_num_structures_lost)
    
    if len(training_set) < maximum_training_set_reduction * initial_train_size:
        logging.warning("Stopped early due to losing more than %f of the inital training set.", maximum_training_set_reduction)
    
    # Check if logging level is debug
    if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
        logging.debug("Initial Bins:")
        logging.debug(inital_bins)
        logging.debug("Final Bins:")
        logging.debug(final_bins)
        
        # Change logging level back at last minute to supress matplotlib
        logging.getLogger().setLevel(logging.INFO)        
        
        # Save a plot of the bins
        plt.figure(figsize=(20, 5))
        plt.title("Initial vs Final Bins of Test Set")
        plt.bar(inital_bins.index.astype(str), inital_bins.values, alpha=0.5, label='Initial')
        plt.bar(final_bins.index.astype(str), final_bins.values, alpha=0.5, label='Final')
        plt.xlabel('Similarity Bin (Max Similarity to Train Set)')
        plt.ylabel('Number of Structures')
        plt.legend()
        plt.savefig('initial_vs_final_bins.png')
        
        # Sae a plot of the pairwise similarities
        plt.figure(figsize=(20, 5))
        plt.title("Pairwise Similarities of Test Set")
        test_sims = sim_matrix.loc[test_set, test_set]
        plt.hist(test_sims.values.flatten(), bins=100)
        plt.xlabel('Pairwise Similarity')
        plt.ylabel('Number of Pairs')
        plt.savefig('pairwise_similarities.png')
        
        # # Save a plot of hte pairwise similarities of the test set spectra
        ### TODO: Fix this plot
        # plt.figure(figsize=(20, 5))
        # plt.title("Pairwise Similarities of Test Set Spectra")
        # bin_edges = np.linspace(0,1.0,11)
        # test_sims_bins = np.zeros(10)
        
        # smiles_counts = summary.Smiles.value_counts()
        # for smiles in pd.Series(test_set).unique():
        #     test_sims_bins += np.histogram(sim_matrix.loc[smiles, test_set].values.flatten(), bins=bin_edges)[0] * smiles_counts[smiles]
        
                
        # plt.bar(np.linspace(0, 1, 11)[1:], test_sims_bins, align='edge')
        # plt.xlabel('Pairwise Similarity')
        # plt.ylabel('Number of Pairs')
        # plt.savefig('pairwise_similarities_all_spectra.png')
        
        logging.getLogger().setLevel(logging.DEBUG)
    
    return test_rows, train_rows
    
def sample_structures_smart_asms(summary, num_test_roots=500, test_path_len=3, test_edge_cutoff=0.70, bins=(0.2, 1.0, 17), maximum_training_set_reduction=0.80, datapoints_per_bin=143):
    """
    Randomly sample num_test_roots structures and then iteratively remove structures from the training set
    until each bin has at least 30 structures.
    
    Parameters:
    summary: pd.DataFrame
        The summary dataframe
    num_test_roots: int
        The number of structures to randomly sample. Used as the root for the random walks
    test_path_len: int
        The length of the random walk from each root
    test_edge_cutoff: float 
        The minimum similarity to be considered a neighbor in the random walk
    bins: tuple
        The bins to use for the similarity matrix
    maximum_training_set_reduction: float
        The maximum fraction of the training set that can be removed
    datapoints_per_bin: int
        The number of datapoints to have in each bin
    """
    
    # Remove BMDMS-NP from contention
    summary = summary.loc[summary['GNPS_library_membership'] != 'BMDMS-NP']
    
    if  missing_structure_check(summary):
        raise ValueError("Found missing structures in summary.")
    summary['InChIKey_smiles_14'] = summary['InChIKey_smiles'].str[:14]
    
    # Now that we have our "test set" we want remove training points to force the max train-test similarity
    unique_structures = list(summary['Smiles'].unique())
    logging.info("Computing Test Set By Sampling the Tail of the Similarity Matrix")
    logging.info("Computing Fingerprints")
    fps = [Chem.RDKFingerprint(Chem.MolFromSmiles(x)) for x in tqdm(unique_structures)]

    # Create similarity matrix
    logging.info("Computing Similarity Matrix")
    sim_matrix = np.zeros((len(fps), len(fps)))
    for i, fp1 in enumerate(tqdm(fps)):
        sim_matrix[i, i] = -1
        for j in range(i+1, len(fps)):
            fp2 = fps[j]
            sim_matrix[i, j] = DataStructs.TanimotoSimilarity(fp1, fp2)
            sim_matrix[j, i] = sim_matrix[i, j]
    sim_matrix = pd.DataFrame(sim_matrix, index=unique_structures, columns=unique_structures)
    
    # For bins with less than datapoints_per_bin structures, sample all
    test_set = []
    training_set = unique_structures
    train_test_similarity = sim_matrix.loc[test_set, training_set]
    train_test_bins = pd.cut(train_test_similarity.max(axis=1), np.linspace(bins[0], bins[1], bins[2]))
    tail_structures = []
    sim_matrix_bins = pd.cut(sim_matrix.max(axis=1), np.linspace(bins[0], bins[1], bins[2]))
    for bin_, group in sim_matrix.groupby(sim_matrix_bins):
        if len(group) > int(datapoints_per_bin):    # Add a bit of padding since these values will shift around once we start removing training structures
            tail_structures.extend(group.sample(int(datapoints_per_bin)).index)
        else:
            tail_structures.extend(group.index)
    
    # Potential Test Smiles (We Will Require that they're not from MassBank to reduce diversity, boosting pairwise similarity in the test set)
    structure_options = summary.loc[['MSBNK' not in x for x in summary['spectrum_id']]].InChIKey_smiles_14
    
    # Random walk matrix
    random_walk_matrix = sim_matrix.copy(deep=True)
    random_walk_matrix[random_walk_matrix < test_edge_cutoff] = 0
    random_walk_matrix[random_walk_matrix >= test_edge_cutoff] = 1
    
    
    # Randomly sample values from the InChIKey_smiles_14 column
    logging.debug("Sampling %i random structures to start the test set.", num_test_roots)
    test_inchi_14 = np.random.choice(structure_options, num_test_roots, replace=False)
    corresponding_smiles = summary.loc[summary['InChIKey_smiles_14'].isin(test_inchi_14)]['Smiles'].unique()
    random_walk_structures = []
    logging.debug("Performing random walks of length %i from the test set roots.", test_path_len)
    for root in corresponding_smiles:
        # Random Walk
        current = root
        for i in range(test_path_len):
            neighbors = random_walk_matrix.loc[current, :]
            neighbors = neighbors[neighbors == 1]
            if len(neighbors) == 0:
                break
            current = np.random.choice(neighbors.index)
            random_walk_structures.append(current)
    # Convert random walk smiles to inchi_14
    random_walk_structures = summary.loc[summary['Smiles'].isin(random_walk_structures)]['InChIKey_smiles_14'].unique()
    test_inchi_14 = np.concatenate([test_inchi_14, random_walk_structures])
            
    logging.info("Sampled %i nodes from random walks.", len(random_walk_structures))
    logging("Sampled %i structures from random walks.", len(random_walk_structures))
              
    tail_inchi_14 = summary.loc[summary['Smiles'].isin(tail_structures)]['InChIKey_smiles_14'].unique()
    test_inchi_14 = np.unique(np.concatenate([test_inchi_14, tail_inchi_14]))
    
    # Get corresponding smiles strings
    test_set = summary.loc[summary['InChIKey_smiles_14'].isin(test_inchi_14)]['Smiles'].unique()
    training_set = summary.loc[~summary['InChIKey_smiles_14'].isin(test_inchi_14)]['Smiles'].unique()
    
    # # To the test set, add 0.05 * num_test_roots of high similarity inchis to ensure these pairs are sampled
    # high_similarity_smiles = sim_matrix.loc[test_set, training_set].max(axis=1).sort_values(ascending=False).head(int(0.05 * num_test_roots)).index.values
    # logging.debug("Adding an additional %i high similarity structures to the test set", len(high_similarity_smiles))
    # test_inchi_14 = np.concatenate([test_inchi_14, summary.loc[summary['Smiles'].isin(high_similarity_smiles)]['InChIKey_smiles_14'].unique()])
    
    # Get corresponding smiles strings
    test_set = summary.loc[summary['InChIKey_smiles_14'].isin(test_inchi_14)]['Smiles'].unique()
    training_set = summary.loc[~summary['InChIKey_smiles_14'].isin(test_inchi_14)]['Smiles'].unique()
    
    train_test_similarity = sim_matrix.loc[test_set, training_set]
    
    train_test_bins = pd.cut(train_test_similarity.max(axis=1), np.linspace(bins[0], bins[1], bins[2]))
    grouped_train_test_sim_matrix = train_test_similarity.max(axis=1).groupby(train_test_bins)
    groups = sorted(list(grouped_train_test_sim_matrix.groups.keys()), key=lambda x: x.right)
    # current_bin_sizes = [len(grouped_train_test_sim_matrix.get_group(group)) for group in groups]
    current_bin_sizes = grouped_train_test_sim_matrix.size()
    bin_to_optimize = 0

    inital_bins = grouped_train_test_sim_matrix.count()

    # train_sims = sim_matrix.loc[training_set, training_set]
    initial_train_size = len(training_set)

    total_num_structures_lost = 0
    
    while min(current_bin_sizes) < datapoints_per_bin and len(training_set) > maximum_training_set_reduction * initial_train_size:
            if bin_to_optimize == len(current_bin_sizes) - 1:
                break
                # This needs more thinking, there may be no structures left to move to the bin
                bin_to_optimize = 0
            
            logging.debug("Optimizing bin: %s", bin_to_optimize)
            logging.debug("Bin size: %i", current_bin_sizes.iloc[bin_to_optimize])
            
            # bins_over_threshold = [x - datapoints_per_bin for x in current_bin_sizes if x > datapoints_per_bin]
            # available_to_sample = sum(bins_over_threshold)
            bin_to_sample_from = -1
            for i in range(bin_to_optimize+1, len(current_bin_sizes)):
                if current_bin_sizes.iloc[i] != 0:
                    bin_to_sample_from = i
                    break
            if bin_to_sample_from == -1:
                logging.warning("No bins have any structures left.")
                break    
            
            available_to_sample = min(datapoints_per_bin - current_bin_sizes.iloc[bin_to_optimize], current_bin_sizes.iloc[bin_to_sample_from])
            logging.debug("%i structures available to sample from bin %s", available_to_sample, bin_to_sample_from)

            # Iteratively shift the highest similarity bins to the left until we have enough structures
            upper_bound = groups[bin_to_optimize].right
            # First we will try with just moving the bin to the right (later we may implement using the largest bin)
            to_move = grouped_train_test_sim_matrix.get_group(groups[bin_to_sample_from]).sample(available_to_sample).index.values

            logging.debug('Moving %i structures from %s to %s', len(to_move), bin_to_sample_from, bin_to_optimize)
            
            # "Move the structures"
            # Remove all of it's neighbors greater than the upper_bound from the trianing set
            structures_to_remove = ((train_test_similarity >= upper_bound).loc[to_move] >= upper_bound).any(axis='index')
            structures_to_remove = structures_to_remove[structures_to_remove].index.values
            training_set = list(set(training_set) - set(structures_to_remove))
            total_num_structures_lost += len(set(structures_to_remove))
            logging.debug("Removing %i structures", len(set(structures_to_remove)))
            
            # Removing structures that go outside of our range, to help reduce confusion later
            structures_to_remove = train_test_similarity.loc[train_test_similarity.max(axis=1) < bins[0]].index.values
            test_set = list(set(test_set) - set(structures_to_remove))
            logging.debug("Removing %i structures from the test set that are below the minimum similarity threshold", len(structures_to_remove))

            # Reinitailize the similarity matrix
            train_test_similarity = sim_matrix.loc[test_set, training_set]
            train_test_bins = pd.cut(train_test_similarity.max(axis=1),  np.linspace(bins[0], bins[1], bins[2]))
            grouped_train_test_sim_matrix = train_test_similarity.max(axis=1).groupby(train_test_bins)
            groups = sorted(list(grouped_train_test_sim_matrix.groups.keys()), key=lambda x: x.right)
            logging.debug("Current Groups: %s", str(grouped_train_test_sim_matrix.groups.keys()))
            # current_bin_sizes = [len(grouped_train_test_sim_matrix.get_group(group)) if group in grouped_train_test_sim_matrix.groups.keys() else 0 for group in groups]
            current_bin_sizes = grouped_train_test_sim_matrix.size()


            print(grouped_train_test_sim_matrix.count())
            
            if len(grouped_train_test_sim_matrix.groups.get(groups[bin_to_optimize])) >= datapoints_per_bin:
                bin_to_optimize += 1

    final_bins = grouped_train_test_sim_matrix.count()
    
    test_rows = summary.loc[summary['Smiles'].isin(test_set)]
    # TEMP: As a temporary measure, remove the BMDMS-NP from the test set
    # The wide range of collision energies is causing issues
    test_rows = test_rows.loc[test_rows['GNPS_library_membership'] != 'BMDMS-NP']
    
    
    train_rows = summary.loc[summary['Smiles'].isin(training_set)]
    logging.info("Number of test rows: %i", len(test_rows))
    logging.info("Number of test structures: %i", len(test_set))
    logging.info("Number of training rows: %i", len(train_rows))
    logging.info("Number of training structures: %i", len(training_set))
    logging.info("Total number of training structures removed: %i", total_num_structures_lost)
    
    if len(training_set) < maximum_training_set_reduction * initial_train_size:
        logging.warning("Stopped early due to losing more than %f of the inital training set.", maximum_training_set_reduction)
    
    # Check if logging level is debug
    if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
        logging.debug("Initial Bins:")
        logging.debug(inital_bins)
        logging.debug("Final Bins:")
        logging.debug(final_bins)
        
        # Change logging level back at last minute to supress matplotlib
        logging.getLogger().setLevel(logging.INFO)        
        
        # Save a plot of the bins
        plt.figure(figsize=(20, 5))
        plt.title("Initial vs Final Bins of Test Set")
        plt.bar(inital_bins.index.astype(str), inital_bins.values, alpha=0.5, label='Initial')
        plt.bar(final_bins.index.astype(str), final_bins.values, alpha=0.5, label='Final')
        plt.xlabel('Similarity Bin (Max Similarity to Train Set)')
        plt.ylabel('Number of Structures')
        plt.legend()
        plt.savefig('initial_vs_final_bins.png')
        
        # Sae a plot of the pairwise similarities
        plt.figure(figsize=(20, 5))
        plt.title("Pairwise Similarities of Test Set")
        test_sims = sim_matrix.loc[test_set, test_set]
        plt.hist(test_sims.values.flatten(), bins=100)
        plt.xlabel('Pairwise Similarity')
        plt.ylabel('Number of Pairs')
        plt.savefig('pairwise_similarities.png')
        
        # # Save a plot of hte pairwise similarities of the test set spectra
        ### TODO: Fix this plot
        # plt.figure(figsize=(20, 5))
        # plt.title("Pairwise Similarities of Test Set Spectra")
        # bin_edges = np.linspace(0,1.0,11)
        # test_sims_bins = np.zeros(10)
        
        # smiles_counts = summary.Smiles.value_counts()
        # for smiles in pd.Series(test_set).unique():
        #     test_sims_bins += np.histogram(sim_matrix.loc[smiles, test_set].values.flatten(), bins=bin_edges)[0] * smiles_counts[smiles]
        
                
        # plt.bar(np.linspace(0, 1, 11)[1:], test_sims_bins, align='edge')
        # plt.xlabel('Pairwise Similarity')
        # plt.ylabel('Number of Pairs')
        # plt.savefig('pairwise_similarities_all_spectra.png')
        
        logging.getLogger().setLevel(logging.DEBUG)
    
    return test_rows, train_rows

def sample_structures_tail(summary, threshold:float=0.5)->pd.DataFrame:
    if missing_structure_check(summary):
        raise ValueError("Found missing structures in summary.")
    
    if threshold < 0 or threshold > 1:
        raise ValueError("Threshold should be between 0 and 1.")
    
    logging.info("Computing Test Set By Sampling the Tail of the Similarity Matrix")
    logging.info("Computing Fingerprints")
    
    unique_structures = list(summary['Smiles'].unique())
    
    # Compute Fingerprints
    all_morgans = [list(Chem.rdMolDescriptors.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(x), 2, 2048)) for x in tqdm(unique_structures)]
    all_daylights = [list(Chem.RDKFingerprint(Chem.MolFromSmiles(x))) for x in tqdm(unique_structures)]
    # Practically, the similarity of two concatenated fingerprints is the average of the similarities of the two fingerprints
    joined_fps = [DataStructs.CreateFromBitString(str(a + b)) for a, b in zip(all_morgans, all_daylights)]
    

    # Create similarity matrix
    logging.info("Computing Similarity Matrix")
    sim_matrix = np.zeros((len(joined_fps), len(joined_fps)))
    for i, fp1 in enumerate(tqdm(joined_fps)):
        sim_matrix[i, i] = -1
        for j in range(i+1, len(joined_fps)):
            fp2 = joined_fps[j]
            sim_matrix[i, j] = DataStructs.TanimotoSimilarity(fp1, fp2)
            sim_matrix[j, i] = sim_matrix[i, j]
    sim_matrix = pd.DataFrame(sim_matrix, index=unique_structures, columns=unique_structures)
    
    test_structures = set(sim_matrix.loc[sim_matrix.max(axis=1) <= threshold].index.values)
    logging.info("Number of test structures: %i", len(test_structures))
    
    # Get original rows
    test_rows = summary.loc[summary['Smiles'].isin(test_structures)]
    logging.info("Number of test rows: %i", len(test_rows))
    
    return test_rows

def sample_structures_smart_tail(summary, datapoints_per_bin:int=30, num_bins:int=20, pairwise_similarity_compensation=True):
    if missing_structure_check(summary):
        raise ValueError("Found missing structures in summary.")
    
    maximum_training_set_reduction = 0.80 # If we lose more than 0.xx(|train_set|), stop reassigning
    inital_margin = 1.8 # The inital margin for the count in each similarity bing
    
    # Check arg types 
    assert isinstance(datapoints_per_bin, int), "datapoints_per_bin should be an integer."
    assert isinstance(num_bins, int), "num_bins should be an integer."
   
    logging.info("Computing Test Set By Sampling the Tail of the Similarity Matrix and Iteratively Adding Points")
    logging.info("Computing Fingerprints")
    
    unique_structures = list(summary['Smiles'].unique())
    
    # Compute Fingerprints (RDK only)
    fps = [Chem.RDKFingerprint(Chem.MolFromSmiles(x)) for x in tqdm(unique_structures)]

    # Create similarity matrix
    logging.info("Computing Similarity Matrix")
    sim_matrix = np.zeros((len(fps), len(fps)))
    for i, fp1 in enumerate(tqdm(fps)):
        sim_matrix[i, i] = -1
        for j in range(i+1, len(fps)):
            fp2 = fps[j]
            sim_matrix[i, j] = DataStructs.TanimotoSimilarity(fp1, fp2)
            sim_matrix[j, i] = sim_matrix[i, j]
    sim_matrix = pd.DataFrame(sim_matrix, index=unique_structures, columns=unique_structures)

    """Generate the inital test set using the "freebies" or points 
    already have a maximum similarity to other points corresponding to our bins"""
    logging.info("Generating initial test set")
    # Bin all structures
    bins = pd.cut(sim_matrix.max(axis=1), np.linspace(0, 1, 21))

    # For bins with less than datapoints_per_bin structures, sample all
    sampled_structures = []
    for bin_, group in sim_matrix.groupby(bins):
        if len(group) > int(datapoints_per_bin * inital_margin):    # Add a bit of padding since these values will shift around once we start removing training structures
            sampled_structures.extend(group.sample(int(datapoints_per_bin * inital_margin)).index)
        else:
            sampled_structures.extend(group.index)
            
    test_set = sampled_structures
    training_set = list(set(unique_structures) - set(test_set))
    
    # To the test set, add 0.05 * num_test_points of high similarity inchis to ensure these pairs are sampled
    high_similarity_smiles = sim_matrix.loc[test_set, training_set].max(axis=1).sort_values(ascending=False).head(int(datapoints_per_bin/2)).index.values
    logging.debug("Adding an additional %i high similarity structures to the test set", len(high_similarity_smiles))
    test_set = list(np.concatenate([np.unique(test_set), np.unique(high_similarity_smiles)]))
    training_set = list(set(unique_structures) - set(test_set))
    
    
    """Iteratively mine out structures from the training set until each bin is at the desired size"""
    logging.info("Iteratively mining out structures")
    train_test_similarity = sim_matrix.loc[test_set, training_set]

    train_test_bins = pd.cut(train_test_similarity.max(axis=1), np.linspace(0, 1, 21))
    grouped_train_test_sim_matrix = train_test_similarity.max(axis=1).groupby(train_test_bins)
    smallest_bin = train_test_similarity.max(axis=1).groupby(train_test_bins).count().idxmin()

    inital_bins = grouped_train_test_sim_matrix.count()

    train_sims = sim_matrix.loc[training_set, training_set]
    initial_train_size = len(training_set)

    total_num_structures_lost = 0
    
    while len(grouped_train_test_sim_matrix.groups.get(smallest_bin)) < datapoints_per_bin and len(training_set) > maximum_training_set_reduction * initial_train_size:
            logging.debug("Smallest bin: %s", smallest_bin)
            logging.debug("Smallest bin size: %i", len(grouped_train_test_sim_matrix.groups.get(smallest_bin)))
            
            number_to_mine = datapoints_per_bin - len(grouped_train_test_sim_matrix.groups.get(smallest_bin))  # Desired number of structures
            
            # Get the set of structures in the training set with the lowest average similairty (implies least number of additional training points removed)
            # that is greater than the upper bound of our similarity bin
            upper_bound = smallest_bin.right
            candidates_to_move = (train_sims.loc[train_sims.max(axis=1) > upper_bound, :]  > upper_bound).sum(axis=1).sort_values(ascending=True)
            to_move = candidates_to_move.head(number_to_mine).index.values
            
            if pairwise_similarity_compensation:
                # Sample some of the high similarity neighbors of the points that we're moving.
                high_sim_neighbors = train_sims.loc[candidates_to_move.head(number_to_mine).index.values, :]
                high_sim_neighbors = high_sim_neighbors.loc[:, high_sim_neighbors.max(axis=0) > 0.8]
                num_to_sample = int(max(len(high_sim_neighbors), 5))    # Bring a friend :)
                high_sim_neighbors = high_sim_neighbors.sample(num_to_sample)
                if len(high_sim_neighbors) > 0:
                    logging.info("Adding %i high similarity neighbors to the test set", len(high_sim_neighbors))
                else:
                    logging.warning("No high similarity neighbors found but pairwise_similarity_compensation=True was specified.")
                
                to_move = list(set(to_move) | set(high_sim_neighbors.index))

            logging.debug('Moving %i structures to %s', len(to_move), smallest_bin)
            
            # "Move the structures"
            # Remove all of it's neighbors greater than the upper_bound from the trianing set
            structures_to_remove = ((train_sims >= upper_bound).loc[to_move] >= upper_bound).any(axis='index')
            structures_to_remove = structures_to_remove[structures_to_remove].index.values
            training_set = list(set(training_set) - set(structures_to_remove))
            total_num_structures_lost += len(set(structures_to_remove) - set(to_move))
            logging.debug("Removing %i structures", len(set(structures_to_remove) - set(to_move)))
            # Add the new structures to the test set
            test_set.extend(to_move)

            # Reinitailize the similarity matrix
            train_test_similarity = sim_matrix.loc[test_set, training_set]
            train_test_bins = pd.cut(train_test_similarity.max(axis=1), np.linspace(0, 1, 21))
            grouped_train_test_sim_matrix = train_test_similarity.max(axis=1).groupby(train_test_bins)
            smallest_bin = train_test_similarity.max(axis=1).groupby(train_test_bins).count().idxmin()

            train_sims = sim_matrix.loc[training_set, training_set]
            
            print(grouped_train_test_sim_matrix.count())

    final_bins = grouped_train_test_sim_matrix.count()
    
    test_rows = summary.loc[summary['Smiles'].isin(test_set)]
    train_rows = summary.loc[summary['Smiles'].isin(training_set)]
    logging.info("Number of test rows: %i", len(test_rows))
    logging.info("Number of test structures: %i", len(test_set))
    logging.info("Number of training rows: %i", len(train_rows))
    logging.info("Number of training structures: %i", len(training_set))
    logging.info("Total number of training structures removed: %i", total_num_structures_lost)
    logging.info("Smallest bin size: %i", len(grouped_train_test_sim_matrix.groups.get(smallest_bin)))
    
    if len(training_set) < maximum_training_set_reduction * initial_train_size:
        logging.warning("Stopped early due to losing more than %f of the inital training set.", maximum_training_set_reduction)
    
    # Check if logging level is debug
    if logging.getLogger().getEffectiveLevel() == logging.DEBUG:
        logging.debug("Initial Bins:")
        logging.debug(inital_bins)
        logging.debug("Final Bins:")
        logging.debug(final_bins)
        
        # Change logging level back at last minute to supress matplotlib
        logging.getLogger().setLevel(logging.INFO)        
        
        # Save a plot of the bins
        plt.figure(figsize=(20, 5))
        plt.title("Initial vs Final Bins of Test Set")
        plt.bar(inital_bins.index.astype(str), inital_bins.values, alpha=0.5, label='Initial')
        plt.bar(final_bins.index.astype(str), final_bins.values, alpha=0.5, label='Final')
        plt.xlabel('Similarity Bin (Max Similarity to Train Set)')
        plt.ylabel('Number of Structures')
        plt.legend()
        plt.savefig('initial_vs_final_bins.png')
        
        # Sae a plot of the pairwise similarities
        plt.figure(figsize=(20, 5))
        plt.title("Pairwise Similarities of Test Set")
        test_sims = sim_matrix.loc[test_set, test_set]
        plt.hist(test_sims.values.flatten(), bins=100)
        plt.xlabel('Pairwise Similarity')
        plt.ylabel('Number of Pairs')
        plt.savefig('pairwise_similarities.png')
        
        logging.getLogger().setLevel(logging.DEBUG)
    
    return test_rows, train_rows

def sample_structures_umbrella(summary, datapoints_per_bin:int=30, num_bins:int=20):
    if missing_structure_check(summary):
        raise ValueError("Found missing structures in summary.")
    
    maximum_training_set_reduction = 0.60 # If we lose more than 0.xx(|train_set|), stop reassigning
    inital_margin = 1.8 # The inital margin for the count in each similarity bing
    
    # Check arg types 
    assert isinstance(datapoints_per_bin, int), "datapoints_per_bin should be an integer."
    assert isinstance(num_bins, int), "num_bins should be an integer."
   
    logging.info("Computing Test Set By Sampling the Tail of the Similarity Matrix and Iteratively Adding Points")
    logging.info("Computing Fingerprints")
    
    unique_structures = list(summary['Smiles'].unique())
    
    # Compute Fingerprints (RDK only)
    fps = [Chem.RDKFingerprint(Chem.MolFromSmiles(x)) for x in tqdm(unique_structures)]

    # Create similarity matrix
    logging.info("Computing Similarity Matrix")
    sim_matrix = np.zeros((len(fps), len(fps)))
    for i, fp1 in enumerate(tqdm(fps)):
        sim_matrix[i, i] = -1
        for j in range(i+1, len(fps)):
            fp2 = fps[j]
            sim_matrix[i, j] = DataStructs.TanimotoSimilarity(fp1, fp2)
            sim_matrix[j, i] = sim_matrix[i, j]
    sim_matrix = pd.DataFrame(sim_matrix, index=unique_structures, columns=unique_structures)

    # Candidate test points are those whose similarity histogram contains a point in every bin
    sim_histogram = {}
    candidate_test_points = set()
    for index, row in sim_matrix.iterrows():
        sim_histogram[index] = np.histogram(row, np.linspace(0, 1, num_bins+1))[0]
        if sum(sim_histogram[index] == 0)/num_bins < 0.5:   # We want these data points to fulfill at least half of the bins
            candidate_test_points.add(index)
    

    # Bin all structures
    bins = pd.cut(sim_matrix.max(axis=1), np.linspace(0, 1, 21))

    # For bins with less than datapoints_per_bin structures, sample all
    sampled_structures = []
    for bin_, group in sim_matrix.groupby(bins):
        if len(group) > int(datapoints_per_bin * inital_margin):    # Add a bit of padding since these values will shift around once we start removing training structures
            sampled_structures.extend(group.sample(int(datapoints_per_bin * inital_margin)).index)
        else:
            sampled_structures.extend(group.index)
            
    test_set = sampled_structures
    training_set = list(set(unique_structures) - set(test_set))        
    # # Initialize empty test_set
    # test_set = []
    # # Initialize complete training_test
    # training_set = list(sim_matrix.index)
    
    """Iteratively mine out structures from the training set until each bin is at the desired size"""
    logging.info("Iteratively mining out structures")
    train_test_similarity = sim_matrix.loc[test_set, training_set]

    train_test_bins = pd.cut(train_test_similarity.max(axis=1), np.linspace(0, 1, num_bins+1))
    grouped_train_test_sim_matrix = train_test_similarity.max(axis=1).groupby(train_test_bins)
    bin_counts = train_test_similarity.max(axis=1).groupby(train_test_bins).count()
    # Start with the highest last bin, this prevents starting with the most different bin
    min_count = bin_counts.min()
    smallest_bin = np.arange(len(bin_counts))[bin_counts == min_count][-1]

    inital_bins = grouped_train_test_sim_matrix.count()

    train_sims = sim_matrix.loc[training_set, training_set]
    initial_train_size = len(training_set)

    total_num_structures_lost = 0
    
    while len(grouped_train_test_sim_matrix.groups.get(smallest_bin)) < datapoints_per_bin and len(training_set) > maximum_training_set_reduction * initial_train_size:
            logging.debug("Smallest bin: %s", smallest_bin)
            logging.debug("Smallest bin size: %i", len(grouped_train_test_sim_matrix.groups.get(smallest_bin)))
            
            number_to_mine = datapoints_per_bin - len(grouped_train_test_sim_matrix.groups.get(smallest_bin))  # Desired number of structures
            
            # Get the set of structures in the training set with the lowest average similairty (implies least number of additional training points removed)
            # that is greater than the upper bound of our similarity bin
            upper_bound = smallest_bin.right
            candidates_to_move = (train_sims.loc[train_sims.max(axis=1) > upper_bound, :]  > upper_bound).median(axis=1).sort_values(ascending=True)
            # Only move structures that are candidates (have good neighbor prosepects)
            logging.debug("Number of candidate test points: %i", len(candidate_test_points))
            candidates_to_move = candidates_to_move[candidates_to_move.index.isin(candidate_test_points)]
            logging.debug("Number of candidate test points in candidates_to_move: %i", len(candidates_to_move))
            to_move = candidates_to_move.head(number_to_mine).index.values
            
            # Sample the neighbor structures of the points that we're moving.
            sampled_neighbors = []
            sim_index = 1
            for structure in to_move:
                while sim_index < num_bins+1:
                    neighbor_lb = np.linspace(0, 1, num_bins+1)[sim_index - 1]
                    neighbor_ub = np.linspace(0, 1, num_bins+1)[sim_index]
                    neighbors = sim_matrix.loc[structure, (sim_matrix.loc[structure] >= neighbor_lb) & (sim_matrix.loc[structure] < neighbor_ub)].index
                    to_append = sim_matrix[neighbors].mean(axis=1).idxmin()
                    if pd.isna(to_append):
                        logging.debug("No neighbors found for %s in bin %i", structure, sim_index)
                        logging.debug("Dataframe: %s", sim_matrix.loc[neighbors])
                        logging.debug("Dataframe.mean(axis=1): %s", sim_matrix.loc[neighbors].mean(axis=1))
                        logging.debug("Current sampled_neighbors: %s", sampled_neighbors)
                        sim_index += 1
                        continue
                    sampled_neighbors.append(to_append)
                    sim_index += 1
                    
            sampled_neighbors = list(set(sampled_neighbors))
            logging.debug('Moving %i neighbors', len(sampled_neighbors))
            logging.debug('Moving %i structures to %s', len(to_move), smallest_bin)
            logging.debug('Moving neighbors: %s', str(sampled_neighbors))
            logging.debug('Moving structures: %s', str(to_move))
            
    
            # to_move = list(set(to_move) | set(sampled_neighbors))
            
            # "Move the structures"
            # Remove all of it's neighbors greater than the upper_bound from the training set
            structures_to_remove = ((train_sims >= upper_bound).loc[to_move] >= upper_bound).any(axis='index')
            structures_to_remove = structures_to_remove[structures_to_remove].index.values
            training_set = list(set(training_set) - set(structures_to_remove))
            total_num_structures_lost += len(set(structures_to_remove) - set(to_move))
            logging.debug("Removing %i structures", len(set(structures_to_remove) - set(to_move)))
            # Add the new structures to the test set
            test_set.extend(to_move)

            # Reinitailize the similarity matrix
            train_test_similarity = sim_matrix.loc[test_set, training_set]
            train_test_bins = pd.cut(train_test_similarity.max(axis=1), np.linspace(0, 1, num_bins+1))
            grouped_train_test_sim_matrix = train_test_similarity.max(axis=1).groupby(train_test_bins)
            smallest_bin = train_test_similarity.max(axis=1).groupby(train_test_bins).count().idxmin()

            train_sims = sim_matrix.loc[training_set, training_set]
            
            print(grouped_train_test_sim_matrix.count())

    final_bins = grouped_train_test_sim_matrix.count()
    
    test_rows = summary.loc[summary['Smiles'].isin(test_set)]
    train_rows = summary.loc[summary['Smiles'].isin(training_set)]
    logging.info("Number of test rows: %i", len(test_rows))
    logging.info("Number of test structures: %i", len(test_set))
    logging.info("Number of training rows: %i", len(train_rows))
    logging.info("Number of training structures: %i", len(training_set))
    logging.info("Total number of training structures removed: %i", total_num_structures_lost)
    logging.info("Smallest bin size: %i", len(grouped_train_test_sim_matrix.groups.get(smallest_bin)))
    
    if len(training_set) < maximum_training_set_reduction * initial_train_size:
        logging.warning("Stopped early due to losing more than %f of the inital training set.", maximum_training_set_reduction)
    
    return test_rows, train_rows

def select_test(args):
    """This function takes in an input_csv, an input_mgf, and the number of test points
        to be randomly generated. The function selects these test points and writes
        them to test_rows.csv, test_rows.mgf, train_rows.csv, and train_rows.mgf.

    Args:
        input_csv (_type_): _description_
        input_similarities (_type_): _description_
        input_mgf (_type_): _description_
    """
    
    # Decompose args
    input_csv = args.input_csv
    input_mgf = args.input_mgf
    num_test_points = args.num_test_points
    sampling_strategy = args.sampling_strategy
    threshold = float(args.threshold) if args.threshold is not None else None
    
    if sampling_strategy not in ['basic_sampling_scheme', 'random', 'structure', 'tail', 'umbrella', 'structure_smart', 'random_tail_hybrid', 'sample_structures_smart_inchikey']:
        raise ValueError(f'Sampling strategy should be either "random", "structure", "tail", or "umbrella" but got {sampling_strategy}.')
    
    print("Loading Summary", flush=True)
    summary = pd.read_csv(input_csv)#.head(100) #HEAD FOR DEBUG
    print("Loading Spectra", flush=True)
    spectra  = IndexedMGF(input_mgf, index_by_scans=False)
    
    if sampling_strategy in ['random', 'structure']:
        try:
            num_test_points = float(num_test_points)
        except Exception as e:
            print("num_test_points should be an integer > 1, or a float within (0,1].")
            raise e
        
        if num_test_points <= 0:
            raise ValueError('Number of test points must be greater than 0.')
        elif num_test_points < 1:
            # Assume it's a percentage
            num_test_points = int(num_test_points * len(summary))
        else:
            # Assume it's an integer
            num_test_points = int(num_test_points)
        
        if num_test_points > len(summary):
            raise ValueError('Number of Test Points is greater than number of rows in input_csv.')
        
    if sampling_strategy not in ['random', 'structure', 'sample_structures_smart_inchikey'] and \
        args.input_similarities is not None:
            raise ValueError("Input similarities are not required for this sampling strategy.")

    # Select Test Points
    if sampling_strategy == 'random':
        logging.info("Sampling randomly")
        test_rows = sample_spectra(summary, num_test_points)
    elif sampling_strategy == 'structure':
        logging.info("Sampling by structure")
        # test_rows = sample_structures(summary.iloc[0:500], 10)    #DEBUG
        test_rows, train_rows = sample_structures(summary, num_test_points)
    elif sampling_strategy == 'structure_smart':
        logging.info("Sampling by structure smart")
        # test_rows, train_rows = sample_structures_smart(summary)
        test_rows, train_rows = sample_structures_smart_asms(summary, test_edge_cutoff=threshold)
    elif sampling_strategy == 'sample_structures_smart_inchikey':
        logging.info("Sampling by structure smart inchikey")
        test_rows, train_rows = sample_structures_smart_inchikey(summary, args.input_similarities, test_edge_cutoff=threshold)
    elif sampling_strategy == 'tail':
        logging.info("Sampling by structure tail")
        test_rows, train_rows = sample_structures_smart_tail(summary)
    elif sampling_strategy == 'umbrella':
        logging.info("Sampling by structure umbrella")
        test_rows, train_rows = sample_structures_umbrella(summary)
    elif sampling_strategy == "random_tail_hybrid":
        logging.info("Sampling by random tail hybrid")
        test_rows, train_rows = random_tail_hybrid(summary)
    elif sampling_strategy == "basic_sampling_scheme":
        test_rows, train_rows = basic_sampling_scheme(summary)
        
    logging.info("Total number of test rows: %i", len(test_rows))
        
    # Only some of the splitting methodologies currently implement a train set
    if 'train_rows' not in vars():
        train_rows = summary.loc[~ summary.spectrum_id.isin(test_rows.spectrum_id)]
        
    logging.info("Total number of potential training rows: %i", len(train_rows))
    test_ids = test_rows['spectrum_id'].values
    train_ids = train_rows['spectrum_id'].values
    
    # Get validations set from 500 inchikeys
    logging.info("Splitting into Train and Validation Sets")
    train_rows['temp_col'] = train_rows['InChIKey_smiles'].str[:14]
    random_sample = np.random.choice(train_rows['temp_col'].unique(), size=500, replace=False)
    val_rows   = train_rows.loc[train_rows['temp_col'].isin(random_sample)]
    train_rows = train_rows.loc[~train_rows['temp_col'].isin(random_sample)]
    logging.info("Number of validation rows: %i", len(val_rows))
    logging.info("Number of training rows: %i", len(train_rows))
    logging.info("Number of validation structures: %i", len(val_rows['temp_col'].unique()))
    logging.info("Number of training structures: %i", len(train_rows['temp_col'].unique()))
    val_rows.drop(columns=['temp_col'], inplace=True)
    train_rows.drop(columns=['temp_col'], inplace=True)

    all_pairs_similarities = pd.read_csv(args.input_similarities, index_col=0)

    logging.info("Writing Test Rows")
    start_time = time()
    # Write Test Rows
    test_rows.to_csv('test_rows.csv', index=False)
    test_inchis = test_rows['InChIKey_smiles'].str[:14].unique()
    mgf.write(spectra[[str(spectrum_id) for spectrum_id in test_rows['spectrum_id']]],  output='test_rows.mgf')
    # Write subsection of train-test similarity matrix
    test_similaries = all_pairs_similarities.loc[test_inchis, test_inchis]
    test_similaries.to_csv('test_similarities.csv')
    logging.info("Done in %.2f seconds.", time() - start_time)
    
    logging.info("Writing Train Rows.")
    start_time = time()
    # Write Train Rows
    train_rows.to_csv('train_rows.csv', index=False)
    mgf.write(spectra[[str(spectrum_id) for spectrum_id in train_rows['spectrum_id']]],  output='train_rows.mgf')
    # Write subsection of train-test similarity matrix
    train_similaries = all_pairs_similarities.loc[train_rows['InChIKey_smiles'].str[:14].unique(), train_rows['InChIKey_smiles'].str[:14].unique()]
    train_similaries.to_csv('train_similarities.csv')
    logging.info("Done in %.2f seconds.", time() - start_time)

    logging.info("Writing Validation Rows.")
    start_time = time()
    # Write Validation Rows
    val_rows.to_csv('val_rows.csv', index=False)
    val_inchis = val_rows['InChIKey_smiles'].str[:14].unique()
    mgf.write(spectra[[str(spectrum_id) for spectrum_id in val_rows['spectrum_id']]], output='val_rows.mgf')
    # Write subsection of train-test similarity matrix
    val_similaries = all_pairs_similarities.loc[val_inchis, val_inchis]
    val_similaries.to_csv('val_similarities.csv')
    logging.info("Done in %.2f seconds.", time() - start_time)

    # Write train-test similarity matrix
    logging.info("Writing Train-Test Similarity Matrix")
    start_time = time()
    all_pairs_similarities.loc[train_rows['InChIKey_smiles'].str[:14].unique(), test_inchis].to_csv('train_test_similarities.csv')
    logging.info("Done in %.2f seconds.", time() - start_time)
    
    # Write training and test rows to json
    logging.info("Writing Train Spectra to JSON")
    start_time = time()
    synchronize_spectra_to_json('train_rows.mgf', 'train_rows.json')
    logging.info("Done in %.2f seconds.", time() - start_time)

    logging.info("Writing Test Spectra to JSON")
    start_time = time()
    synchronize_spectra_to_json('test_rows.mgf', 'test_rows.json')
    logging.info("Done in %.2f seconds.", time() - start_time)

    logging.info("Writing Validation Spectra to JSON")
    start_time = time()
    synchronize_spectra_to_json('val_rows.mgf', 'val_rows.json')
    logging.info("Done in %.2f seconds.", time() - start_time)


def main():
    parser = argparse.ArgumentParser(description='Create Parallel Parameters.')
    parser.add_argument('--input_csv', type=str, help='Input CSV')
    parser.add_argument('--input_similarities', type=str, help='Input Similarities')
    parser.add_argument('--input_mgf', type=str, help='Input MGF')
    parser.add_argument('--num_test_points', help='Number of Test Points')
    parser.add_argument('--sampling_strategy', help='Sampling Strategy')
    parser.add_argument('--threshold', help='Threshold for sampling by structure', default=None)
    parser.add_argument('--debug', action='store_true', help='Debug')
    args = parser.parse_args()
    
    # Setup logging
    logging_level = logging.DEBUG if args.debug else logging.INFO
    logging.basicConfig(filename="test_selection.log", level=logging_level)   # TODO: Convert everything else to logging
    
    select_test(args)
    
if __name__ == '__main__':
    main()
    
class FakeArgs:
    def __init__(self, input_csv, input_mgf, num_test_points):
        self.input_csv = input_csv
        self.input_mgf = input_mgf
        self.num_test_points = num_test_points
    
def test_select_test():
    """This function tests the select_test function.
    """
    input_csv = '../data/test_data/test_summary.csv'
    input_mgf = '../data/test_data/test.mgf'
    num_test_points = 2
    
    faux_args = FakeArgs(input_csv, input_mgf, num_test_points)
    
    select_test(faux_args)
    
    test_rows = pd.read_csv('test_rows.csv')
    train_rows = pd.read_csv('train_rows.csv')
    
    assert len(test_rows) == num_test_points
    assert len(train_rows) == 8
    
    os.remove('test_rows.csv')
    os.remove('test_rows.mgf')
    os.remove('train_rows.csv')
    os.remove('train_rows.mgf')