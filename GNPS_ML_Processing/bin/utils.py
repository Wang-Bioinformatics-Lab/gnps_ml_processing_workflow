import csv
import importlib
from rdkit import Chem
from rdkit.Chem import MolFromSmiles
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.Chem import DataStructs
from rdkit.DataStructs.cDataStructs import BulkTanimotoSimilarity
from rdkit.Chem import MolStandardize, rdMolDescriptors, MolFromSmiles
import pandas as pd
from pyteomics.mgf import IndexedMGF
from tqdm import tqdm
import networkx as nx
import numpy as np
from networkx.algorithms import community
from typing import Tuple
import vaex
import os
# import dask.dataframe as dd
import tempfile
import json
from pandarallel import pandarallel
from joblib import Parallel, delayed

PARALLEL_WORKERS = 32 # Number of workers to use for parallel processing

class Wrapper(object):
    def __init__(self, method_name, module_name):
        """Because boost functions are not picklable, we create a wrapper around the functions. 
            This class takes two arguments: method_name and module_name and imports the module and 
            wraps the method.
            ---
            Example Usage: Wrapper("MolFromSmiles", "rdkit.Chem")
        Args:
            method_name (str): Method to be imported
            module_name (str): Module the method is located in.
        """
        self.method_name = method_name
        # self.module_name = module_name
        self.module = importlib.import_module(module_name)

    @property
    def method(self):
        return getattr(self.module, self.method_name)
    
    def __call__(self, *args, **kwargs):
        return self.method(*args, **kwargs)

def split_cliques(num_nodes:int,ids:list,training_fraction: float=0.8, early_stopping=None):
    """ A function returning the indices of cliques in a graph that maximize the number of total datapoints in a training set.
        Will only return a list of clique indices whose total number of elements is less than the desired fraction.
        Analogous to solving the 0/1 knapsack problem where the elements added to the knapsack have value=weight=clique size.

    Args:
        num_nodes (int): The number of nodes in each clique.
        ids (list): A unique integer ID for each clique,
        training_fraction (float, optional): Value between 0 and 1, specifies the maximum fraction of data points that 
        can be used in the training set. Defaults to 0.8.

    Returns:
        int, list: Returns a tuple of the number of training examples in the split, and the ids for the cliques containing the training examples.
    """
    maximum_training_instances = int(training_fraction * sum(num_nodes))
    n = len(num_nodes)
    
    dp = [[0,[]] for i in range(maximum_training_instances+1)]
    for i in tqdm(range(1,n+1)):
        for w in range(maximum_training_instances, 0, -1):
            if num_nodes[i-1] <= w:
                if dp[w][0] > dp[w-num_nodes[i-1]][0]+num_nodes[i-1]:
                    dp[w][0] = dp[w][0]
                    dp[w][1] = dp[w][1]
                else:
                    dp[w][0] = dp[w-num_nodes[i-1]][0]+num_nodes[i-1]
                    dp[w][1] = dp[w-num_nodes[i-1]][1] + [ids[i-1]]
            if early_stopping is not None:
                if dp[w][0] > early_stopping:
                    return dp[w][0], dp[w][1]
                    
    return dp[maximum_training_instances][0], dp[maximum_training_instances][1]

def get_splits_index_only(scans:list, similarity_table: pd.DataFrame, similarity_metrics:list, similarity_threshold: list, training_fraction: float = 0.80, max_retries: int = 3) -> Tuple[list,list]:  
    """This function returns the scan numbers that can be used to slice a training set given a precomputed similarity metric.
        If the function fails to find a split that contains 90% the desired number of training points, it employes the asynchronous 
        fluid communities algorithm to split the largest similarity clique in two up to max_retries times. Note that when this
        happens, some nodes will be removed from both the training and the test set.
        
        If the algorithm fails to find a satisfactory split, it will return a suboptimal split.

    Args:
        scans (list): List of scans in the mgf file to be split.
        similarity_table (pd.DataFrame): A table with headings spectrumid1, spectrumid1, and the similarity metrics specied in similarity_metrics.
        similarity_metrics (list): The names of the columns to be used as similarity measures.
        similarity_threshold (list, optional): Threshold at which point two data points will be considered similar based on the similarity metrics.        
        training_fraction (float, optional): Value between 0 and 1, specifies minimum number of examples used for training.. Defaults to 0.80.
        max_retries (int, optional): Number of times the largest clique will be broken apart. Defaults to 3.

    Returns:
        tuple (list, list): A tuple of lists containing the scan numbers of the training and test sets respectively.
    """    
   
    network = similarity_table

    if np.any(np.array([min(network[x]) for x in similarity_metrics]) > np.array(similarity_threshold)) and similarity_threshold is not None:
        print("Warning: The similarity threshold is below the minimum value in the similarity matrix for one or more of the columns. Has the similarity matrix been trimmed?")
    elif similarity_threshold is not None: 
        network = network[np.bitwise_or.reduce([network[metric] >= threshold for metric,threshold in zip(similarity_metrics, similarity_threshold)])]
        
    # It turns out it is much faster to reverse the adjacency list before building the network so we'll do that instead of transforming the graph driectly
    network = network.loc[:,['spectrumid1','spectrumid2'] + similarity_metrics]
    reversed_network = network.loc[:,['spectrumid2','spectrumid1'] + similarity_metrics]
    reversed_network.columns=['spectrumid1','spectrumid2'] + similarity_metrics
    network = pd.concat([network, reversed_network], axis=0)
    del reversed_network
    G = nx.from_pandas_edgelist(network, 'spectrumid1', 'spectrumid1', edge_attr=similarity_metrics)
    del network
    
    # Add nodes that may have no similarity
    G.add_nodes_from(scans)
    print("Desired Number of Training Instances:", int(len(G) * training_fraction), flush=True)    
    
    for i in range(max_retries):
            # Retry finding best split max_retries time
        graph_size = len(G)
        minimum_allowed_training_indices = 0.9 * int(graph_size * training_fraction)
        print("Minimum Allowed Training Instances:", minimum_allowed_training_indices)
        connected_components         = [c for c in sorted(nx.connected_components(G), key=len, reverse=True)]
        ids                          = list(range(len(connected_components)))
        num_nodes                    = [len(x) for x in connected_components]  # number of nodes in each clique   
        
        num_training_instances, training_ids = split_cliques(num_nodes, ids, training_fraction=training_fraction,early_stopping=None)
        if num_training_instances < minimum_allowed_training_indices:
            print("Warning: Unable to find a dataset split within an acceptable margin of the desired split. \nThis split: {} instances".format(num_training_instances))
            if i != max_retries-1:
                # Partition largest connected component into two distinct sets
                subgraph = G.subgraph(connected_components[0])
                partition = community.asyn_fluidc(subgraph, 2, seed=42)
                node_dict = {node: i for i, comm in enumerate(partition) for node in comm}
                set_a = set([k for k,v in node_dict.items() if v==0])
                set_b = set([k for k,v in node_dict.items() if v==1])
                edges_to_remove = np.array([[a,b] for a,b in nx.edge_boundary(G, set(set_a), set(set_b))])
                # We'll zig-zag down this list to remove close to as few as possible from each set
                to_remove = np.concatenate((edges_to_remove[::2][:,0], edges_to_remove[1::2][:,1]))
                print(len(np.unique(to_remove)))
                print(len(np.unique(edges_to_remove)))
                # G.remove_nodes_from(np.unique(np.array(edges_to_remove)))
                G.remove_nodes_from(to_remove)
                
        else:
            break
        
    print("Optimal Number of Training Instances:", num_training_instances)
    
    # Map back to scan numbers    
    training_indices = [y for x in training_ids for y in list(G.subgraph(connected_components[x]))]
    test_indices     = [x for x in list(G.nodes) if x not in training_indices]
    return training_indices, test_indices

def build_tanimoto_similarity_list(mgf_summary:pd.DataFrame,output_dir:str=None, similarity_threshold=0.50, smiles_column_name='Smiles',num_bits=4096) -> pd.DataFrame:
    """This function computes the all pairs tanimoto similarity on an MGF summary dataframe.

    Args:
        mgf_summary (pd.DataFrame): A dataframe containing a scan column and a smiles column.
        output_dir (str, optional): _description_. Defaults to None.
        similarity_threshold (float, optional): _description_. Defaults to 0.50.
        smiles_column_name (str, optional): _description_. Defaults to 'Smiles'.
        num_bits (int, optional): _description_. Defaults to 4096.

    Returns:
        pd.DataFrame: _description_
    """    

    smiles_list = mgf_summary[smiles_column_name]
    structure_mask = [False if x == "nan" else True for x in smiles_list]
    print("Found {}  structures in {} files".format(sum(structure_mask), len(structure_mask)))
    
    fps = [GetMorganFingerprintAsBitVect(MolFromSmiles(m), 2, num_bits) for m in smiles_list[structure_mask]]   
    # Indices are now non contiguous because the entries without structures are removed
    # This will map back to the original scan number  
    idx_mapping = {i: scan_num for i, scan_num in enumerate(mgf_summary.scan[structure_mask])}
    N = len(fps)

    id1 = []
    id2 = []
    sim = []

    for i in range(N):
        sims = BulkTanimotoSimilarity(fps[i],fps[i+1:])
        for j in range(0,len(sims)):
            if sims[j] > similarity_threshold:
                id1.append(idx_mapping[i])
                id2.append(idx_mapping[i+1+j])  #i + 1 +j because we skip the diagonal
                sim.append(sims[j])
    out = pd.DataFrame({"CLUSTERID1":id1, "CLUSTERID2":id2, 'Tanimoto_Similarity':sim},)
    if output_dir is not None: out.to_csv(output_dir, index_label=False)
    return out

def _sim_helper_csv(fps,
                comparison_fps_lst,
                indices,
                truncate,
                fieldnames,
                similarity_threshold,
                idx_mapping,
                output_file = None,
                progress=False):
    """A short helper function that computes the similarities between fp and comparison fps 
    in parallel.

    Args:
        fps (_type_): _description_
        comparison_fps (_type_): _description_
        idx (_type_): _description_

    Returns:
        _type_: _description_
    """
    if output_file is None:
        temp_file = os.path.join(tempfile.gettempdir(), f"GNPS_PROCESSING/TEMP_SIMILARITY_FILE_{indices[0]}.csv")
    else: 
        temp_file = output_file
    # Wrapping for easy parallelization
    wrapped_bulk_tanimoto_similarity = Wrapper('BulkTanimotoSimilarity', 'rdkit.DataStructs.cDataStructs')
    try:
        with open(temp_file, mode='w') as f:
            temp_writer = csv.DictWriter(f, fieldnames=fieldnames)
            
            # Add optional progress bar
            if progress:
                loop_iterator = tqdm(zip(indices, fps, comparison_fps_lst))
            else: 
                loop_iterator =  zip(indices, fps, comparison_fps_lst)
    
            for idx, fp, comparison_fps in loop_iterator:
                if truncate is not None:
                    n_bins, n_items_per_bin = truncate
                    
                    if not isinstance(n_bins, int) or not isinstance(n_items_per_bin, int):
                        raise ValueError("Expected arg 'truncate' to be None or a tuple of type (int,int).")
                    bin_size = (1.0 - similarity_threshold) / (n_bins-1)    # Reserve the last bin for identical spectra (rest of histogram is left aligned)
                    items_left_per_bin = np.full(n_bins, n_items_per_bin)
                                   
                sims = wrapped_bulk_tanimoto_similarity(fp, comparison_fps)
                for j, this_sim in enumerate(sims):                  
                    if this_sim >= similarity_threshold:
                        # Check if we have filled the bin for this entry.
                        if truncate is not None:
                            bin_index = int((this_sim - similarity_threshold) / bin_size)
                            if items_left_per_bin[bin_index] <= 0:
                                continue
                            items_left_per_bin[bin_index] -= 1
                        for edge_from in idx_mapping[idx]:
                            # Loop over identical structures
                            for edge_to in idx_mapping[idx + 1 + j]: # i + 1 + j because we skip the diagonal
                                row = {
                                    'spectrumid1': edge_from,
                                    'spectrumid2': edge_to,  
                                    'Tanimoto_Similarity': this_sim
                                }
                                temp_writer.writerow(row)
                        if truncate is not None:
                            if sum(items_left_per_bin) == 0:
                                break
                        
        return temp_file
    except Exception as suprise_exception:
        # Perform file cleanup if needed
        if os.path.isfile(temp_file):
            os.remove(temp_file)
        raise suprise_exception
       
def _sim_helper(fps,
                comparison_fps,
                indices,
                truncate,
                fieldnames,
                similarity_threshold,
                idx_mapping,
                output_file=None,
                progress=False):
    if output_file is None:
        temp_file = os.path.join(tempfile.gettempdir(), f"GNPS_PROCESSING/TEMP_SIMILARITY_FILE_{indices[0]}.json")
    else:
        temp_file = output_file
    
    # Create a dictionary to store the data grouped by 'spectrumid1'
    data_dict = {}

    # Wrapping for easy parallelization
    wrapped_bulk_tanimoto_similarity = BulkTanimotoSimilarity
    try:
        for idx, fp in zip(indices, fps):
            if truncate is not None:
                n_bins, n_items_per_bin = truncate

                if not isinstance(n_bins, int) or not isinstance(n_items_per_bin, int):
                    raise ValueError("Expected arg 'truncate' to be None or a tuple of type (int,int).")
                bin_size = (1.0 - similarity_threshold) / (n_bins - 1)  # Reserve the last bin for identical spectra (rest of histogram is left aligned)
                items_left_per_bin = np.full(n_bins, n_items_per_bin)

            sims = wrapped_bulk_tanimoto_similarity(fp, comparison_fps)
            for j, this_sim in enumerate(sims):
                if idx == j:
                    continue
                if this_sim >= similarity_threshold:
                    # Check if we have filled the bin for this entry.
                    if truncate is not None:
                        bin_index = int((this_sim - similarity_threshold) / bin_size)
                        if items_left_per_bin[bin_index] <= 0:
                            continue
                        items_left_per_bin[bin_index] -= 1
                    for edge_from in idx_mapping[idx]:
                        # Loop over identical structures
                        if data_dict.get(edge_from) is None: data_dict[edge_from] = []
                        data_dict[edge_from].append({
                                'Tanimoto_Similarity': this_sim,
                                'spectrumid2': [edge_to for edge_to in idx_mapping[j]]})
                    if truncate is not None:
                        if sum(items_left_per_bin) == 0:
                            break

        # Serialize the dictionary to a JSON file
        with open(temp_file, 'w') as json_file:
            json.dump(data_dict, json_file, indent=2)

    except Exception as e:
        raise e
        print(f"An error occurred: {e}")
    
    return temp_file

def build_tanimoto_similarity_list_precomputed(mgf_summary:pd.DataFrame,
                                               output_file:str, 
                                               similarity_threshold=0.50, 
                                               fingerprint_column_name='Morgan_2048_3', 
                                               truncate=None) -> pd.DataFrame:
    """This function computes the all pairs tanimoto similarity on an MGF summary dataframe.

    Args:
        mgf_summary (pd.DataFrame): A dataframe containing a scan column and a smiles column.
        output_file (str, optional): _description_.
        similarity_threshold (float, optional): _description_. Defaults to 0.50.
        fingerprint_column_name (str, optional): _description_. Defaults to 'Morgan_2048_3'.
        num_bits (int, optional): _description_. Defaults to 4096.
        truncate (int, int): Truncate the similarity matrix for a dataframe to (n_bins, n_items_per_bin) or if None, do not truncated, defaults to None.
        remove_duplicates bool: Remove all similarities with identical tanimoto scores, defaults to True.
    Returns:
        None: The output is written to output file. This file can be large and should be read with an out-of-core library.
    """
    if similarity_threshold < 0:
        raise ValueError("Expected arg 'similarity_threshold' to be between 0 and 1.")
    
    if not isinstance(mgf_summary, pd.DataFrame):
        raise ValueError(f"Expected a Pandas DataFrame but got {type(mgf_summary)}")
    
    org_len = len(mgf_summary)
    
    structure_mask = ~ mgf_summary[fingerprint_column_name].isna()
    mgf_summary = mgf_summary[structure_mask]
    
    grouped_mgf_summary = mgf_summary.groupby('Smiles')    
    mgf_summary = mgf_summary.drop_duplicates(subset='Smiles')
    

    
    pandarallel.initialize(progress_bar=False, nb_workers=PARALLEL_WORKERS, verbose=0)
    fps = mgf_summary[fingerprint_column_name].parallel_apply(lambda x: DataStructs.CreateFromBitString(''.join(str(y) for y in x))).values
        
    print(f"Found {len(mgf_summary)} unique structures in {org_len} files")

    # Indices are now non contiguous because the entries without structures are removed
    # This will map back to the original spectrum_id
    idx_mapping = {idx: group_df['spectrum_id'].values for idx, (_, group_df) in enumerate(grouped_mgf_summary)}
    
    # Cleanup variables to reduce memory overhead
    del mgf_summary
    del grouped_mgf_summary
    
    fieldnames = ['spectrumid1', 'spectrumid2', 'Tanimoto_Similarity']
           
    # _ = _sim_helper_csv([fps[j] for j in range(len(fps))],
    #             [fps[j+1:] for j in range(len(fps))],
    #             [j for j in range(len(fps))],
    #             truncate,
    #             fieldnames,
    #             similarity_threshold,
    #             idx_mapping,
    #             output_file = output_file,
    #             progress = True)
    _ = _sim_helper(fps,
                fps,
                [j for j in range(len(fps))],
                truncate,
                fieldnames,
                similarity_threshold,
                idx_mapping,
                output_file = output_file,
                progress = True)
    
def get_fingerprints(smiles_string):
    """Returns a list of fingerprints for a given smiles string"""
    error_value = np.array([None, None, None, None], dtype=object)
    
    if smiles_string=='nan':
        return error_value
    mol = Chem.MolFromSmiles(str(smiles_string), sanitize=True)
    if mol is None:
        return error_value
    
    # If vaex is using mmultiprocessing, it will implicty convert the return value into a numpy araray. 
    return np.array([list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, useChirality=False, nBits=2048)),
                    list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, useChirality=False, nBits=4096)),
                    list(AllChem.GetMorganFingerprintAsBitVect(mol, 3, useChirality=False, nBits=2048)),
                    list(AllChem.GetMorganFingerprintAsBitVect(mol, 3, useChirality=False, nBits=4096))], dtype=object)

def _generate_fingerprints_pandas_helper(smiles):
    if smiles is not None and smiles != 'nan':
        mol = Chem.MolFromSmiles(str(smiles), sanitize=True)
        if mol is not None:
            return smiles, {'Morgan_2048_2': list(AllChem.GetMorganFingerprintAsBitVect(mol,2,useChirality=False,nBits=2048)), 
                            'Morgan_4096_2': list(AllChem.GetMorganFingerprintAsBitVect(mol,2,useChirality=False,nBits=4096)), 
                            'Morgan_2048_3': list(AllChem.GetMorganFingerprintAsBitVect(mol,3,useChirality=False,nBits=2048)), 
                            'Morgan_4096_3':list(AllChem.GetMorganFingerprintAsBitVect(mol,3,useChirality=False,nBits=4096))}
    return smiles, {'Morgan_2048_2': None,
                    'Morgan_4096_2': None, 
                    'Morgan_2048_3': None, 
                    'Morgan_4096_3': None}

def _generate_fingerprints_pandas(summary, progress=True):
    if progress:
        unique_smiles = tqdm(summary.Smiles.unique())
    else:
        unique_smiles = summary.Smiles.unique()
    mapping_dict = dict(Parallel(n_jobs=-1)(delayed(_generate_fingerprints_pandas_helper)(smiles) for smiles in unique_smiles))
    
    summary[['Morgan_2048_2','Morgan_4096_2','Morgan_2048_3','Morgan_4096_3']] = pd.DataFrame.from_records(summary['Smiles'].map(mapping_dict))
    return summary

def _generate_fingerprints_dask(summary):
    # summary = dd.read_csv(summary, dtype={'Smiles':str,'msDetector':str,'msDissociationMethod':str,'msManufacturer':str,'msMassAnalyzer':str,'msModel':str})
    def _get_fingerprint(mol, params):
        if mol is None:
            return None
        else:
            if params == 'Morgan_2048_2':
                return str(GetMorganFingerprintAsBitVect(mol, 2, nBits=2048))
            elif params == 'Morgan_2048_3':
                return str(GetMorganFingerprintAsBitVect(mol, 3, nBits=2048))
            elif params == 'Morgan_4096_2':
                return str(GetMorganFingerprintAsBitVect(mol, 2, nBits=4096))
            elif params == 'Morgan_4096_3':
                return str(GetMorganFingerprintAsBitVect(mol, 3, nBits=4096))
            else:
                raise ValueError("Invalid fingerprint type")
            
    summary['mol'] = summary.apply(lambda x: MolFromSmiles(str(x['Smiles'])), axis=1, meta=('mol', 'object'))
    summary['Morgan_2048_2'] = summary.apply(lambda x: _get_fingerprint(x['mol'], 'Morgan_2048_2'), axis=1, meta=('Morgan_2048_2', 'str'))
    summary['Morgan_2048_3'] = summary.apply(lambda x: _get_fingerprint(x['mol'], 'Morgan_2048_2'), axis=1, meta=('Morgan_2048_3', 'str'))
    summary['Morgan_4096_2'] = summary.apply(lambda x: _get_fingerprint(x['mol'], 'Morgan_2048_2'), axis=1, meta=('Morgan_4096_2', 'str'))
    summary['Morgan_4096_3'] = summary.apply(lambda x: _get_fingerprint(x['mol'], 'Morgan_2048_2'), axis=1, meta=('Morgan_4096_3', 'str'))
    
    return summary.drop('mol', axis=1)

def generate_fingerprints(summary):
    columns_to_add = ['Morgan_2048_2', 'Morgan_4096_2', 'Morgan_2048_3', 'Morgan_4096_3']
    # Check to make sure that the columns are not already present
    if not all([x in summary.columns for x in columns_to_add]):
        if type(summary) is pd.DataFrame:
            return _generate_fingerprints_pandas(summary)
        elif type(summary) is vaex.dataframe.DataFrameLocal:
            raise NotImplementedError("Vaex not yet implemented")
            # We'll assume it's vaex
            print("Using vaex")
            summary['mol'] = summary.apply(lambda x: Chem.MolFromSmiles(x, sanitize=True), arguments=[summary.Smiles], multiprocessing=False)
            summary['Morgan_2048_2'] = summary.apply((lambda x: list(AllChem.GetMorganFingerprintAsBitVect(x,2,useChirality=False,nBits=2048))), arguments=[summary.mol])
            summary.execute()
            summary.materialize('Morgan_2048_2')
            # summary['fingerprints'] = summary.apply(get_fingerprints, arguments=[summary.Smiles])
            
            # summary['Morgan_2048_2'] = summary.apply((lambda x: x[0]), arguments=[summary.fingerprints])
            # summary['Morgan_4096_2'] = summary.apply((lambda x: x[1]), arguments=[summary.fingerprints])
            # summary['Morgan_2048_3'] = summary.apply((lambda x: x[2]), arguments=[summary.fingerprints])
            # summary['Morgan_4096_3'] = summary.apply((lambda x: x[3]), arguments=[summary.fingerprints])
            # summary.drop('fingerprints', inplace=True)
            # summary.execute()
        # elif type(summary) == dd.core.DataFrame:    
        #     if summary.Smiles.isnull().any().compute():
        #         raise ValueError("Smiles column contains null values")        
        #     _generate_fingerprints_dask(summary)
        
        else: 
            raise ValueError("Summary is not a pandas or dask dataframe but got type {}".format(type(summary)))

    return summary

def neutralize_atoms(smiles):
    """This function takes in an rdkit mol and verifies that it does not have a charge that can be
    neutralized by protonation or deprotonation. If it does, it will neutralize the charge and return the mol.
    Code Source: http://www.rdkit.org/docs/Cookbook.html#neutralizing-molecules
    Returns:
        Union[int,bool,int,rdkit.Chem.rdchem.Mol]: The number of charges removed, whether the resulting mol is 
        both pos and neg charged, the sum of the charges, and the neutralized mol.
    """    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    num_removed_charges = len(at_matches_list)
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    pos = False
    neg = False
    sum_of_charges = 0
    for atom in mol.GetAtoms():
        charge = atom.GetFormalCharge()
        if charge < 0:
            neg = True
        if charge > 0:
             pos = True
        sum_of_charges += charge
    if pos and neg:
        pos_and_neg =True
    else:
        pos_and_neg = False
    
    return num_removed_charges, pos_and_neg, sum_of_charges, Chem.MolToSmiles(mol)

# Code Credit: Yasin El Abiead
def harmonize_smiles_rdkit(smiles, tautomer_limit = 900, skip_tautomerization=False):
    if smiles is None or smiles == 'nan' or smiles == 'None' or smiles == '': return ''
    try:
        smiles = str(smiles)
        # take the largest covalently bound molecule
        smiles_largest = MolStandardize.fragment.LargestFragmentChooser(smiles).prefer_organic
        mol = Chem.MolFromSmiles(smiles_largest)

        if mol is None:
            # The files failed to parse, it should be removed
            return ''
        
        monomass = rdMolDescriptors.CalcExactMolWt(mol)
        if not skip_tautomerization:
            # standardize tautomer
            if monomass < tautomer_limit:
                smiles_largest = MolStandardize.canonicalize_tautomer_smiles(smiles_largest)
                mol = Chem.MolFromSmiles(smiles_largest)

        # remove unnecessary charges
        uc = MolStandardize.charge.Uncharger()
        uncharged_mol = uc.uncharge(mol)
        
        # standardize the molecule
        lfc = MolStandardize.fragment.LargestFragmentChooser()
        standard_mol = lfc.choose(uncharged_mol)

        # remove stereochemistry
        Chem.RemoveStereochemistry(standard_mol)

        # get the standardized SMILES
        standard_smiles = Chem.MolToSmiles(standard_mol)
        return standard_smiles

    except Exception as e:
        print(f"An error occurred with input {smiles}: {e}")
        return ''
    
def INCHI_to_SMILES(inchi):
    if inchi is None or inchi == 'nan': return ''
    try:
        mol = Chem.MolFromInchi(inchi)
        if mol is None:
            return ''
        return Chem.MolToSmiles(mol)
    except:
        return ''
    
def synchronize_spectra(input_path, output_path, summary, progress_bar=True):
    """Reads an MGF file from input_path and generates a new_mgf file in output_path with only the spectra in spectrum_ids.

    Args:
        input_path (str): Path to the input mgf.
        output_path (str): Path to save the output mgf.
        summary (pd.Dataframe): Dataframe with spectrum ids to keep.
    """   
    for col_name in ['spectrum_id', 'scan', 'Charge']:
        if col_name not in summary.columns:
            raise ValueError("Summary must contain columns 'spectrum_id', 'scan', 'charge', and 'Compund_Name'. \n \
                             Instead got columns {}".format(summary.columns))
    
    with open(output_path, 'w') as output_mgf:
        input_mgf = IndexedMGF(input_path)
        
        if progress_bar:
            print("Syncing MGF with summary")
            mapping = tqdm(summary[['spectrum_id','scan','Charge','Compund_Name']].itertuples())
        else:
            mapping = summary[['spectrum_id','scan','Charge','Compund_Name']].itertuples()
        
        for _, title, scan, charge, compound_name in mapping:
            spectra = input_mgf[title]
            if spectra['params']['title'] != title:
                raise ValueError("Sanity Check Failed. Expected specrum identifier did not match mgf spectrum identifier.")
            output_mgf.write("BEGIN IONS\n")
            output_mgf.write("PEPMASS={}\n".format(float(spectra['params']['pepmass'][0])))
            output_mgf.write("CHARGE={}\n".format(charge))
            output_mgf.write("TITLE={}\n".format(spectra['params']['title']))
            output_mgf.write("COMPOUND_NAME={}\n".format(compound_name))
            output_mgf.write("SCANS={}\n".format(scan))

            peaks = zip(spectra['m/z array'], spectra['intensity array'])
            for peak in peaks:
                output_mgf.write("{} {}\n".format(peak[0], peak[1]))

            output_mgf.write("END IONS\n")
                
def generate_parquet_df(input_mgf):
    """
    Details on output format:
    Columns will be [level_0, index, i, i_norm, mz, precmz]
    Index will be spectrum_id
    level_0 is the row index in file
    index is the row index in the spectra
    """
    output = []
    indexed_mgf = IndexedMGF(input_mgf,index_by_scans=True)
    level_0 = 0
    for m in tqdm(indexed_mgf):
        spectrum_id = m['params']['title']
        mz_array = m['m/z array']
        intensity_array = m['intensity array']
        precursor_mz = m['params']['pepmass']
        # charge = m['charge']
        for index, (mz, intensity) in enumerate(zip(mz_array, intensity_array)):
            output.append({'spectrum_id':spectrum_id, 'level_0': level_0, 'index':index, 'i':intensity, 'mz':mz, 'prec_mz':precursor_mz})
            level_0 += 1
                
    output = pd.DataFrame(output)
    return output