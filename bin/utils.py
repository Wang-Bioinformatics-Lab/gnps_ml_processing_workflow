from rdkit.DataStructs.cDataStructs import BulkTanimotoSimilarity
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.Chem import MolFromSmiles
import pandas as pd
from tqdm import tqdm
import pandas as pd
import networkx as nx
import numpy as np
from networkx.algorithms import community 
from typing import Tuple
import json
from rdkit.Chem import DataStructs

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

def build_tanimoto_similarity_list_precomputed(mgf_summary:pd.DataFrame,output_dir:str=None, similarity_threshold=0.50, fingerprint_column_name='Morgan_2048_3') -> pd.DataFrame:
    """This function computes the all pairs tanimoto similarity on an MGF summary dataframe.

    Args:
        mgf_summary (pd.DataFrame): A dataframe containing a scan column and a smiles column.
        output_dir (str, optional): _description_. Defaults to None.
        similarity_threshold (float, optional): _description_. Defaults to 0.50.
        fingerprint_column_name (str, optional): _description_. Defaults to 'Morgan_2048_3'.
        num_bits (int, optional): _description_. Defaults to 4096.

    Returns:
        pd.DataFrame: _description_
    """    
    
    structure_mask = ~ mgf_summary[fingerprint_column_name].isna()
    print("Found {}  structures in {} files".format(sum(structure_mask), len(structure_mask)))
    
    fps = [DataStructs.CreateFromBitString(x) for x in mgf_summary[fingerprint_column_name][structure_mask]]
    # Indices are now non contiguous because the entries without structures are removed
    # This will map back to the original scan number  
    idx_mapping = {i: spectrum_id for i, spectrum_id in enumerate(mgf_summary.spectrum_id[structure_mask])}
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
    out = pd.DataFrame({"spectrumid1":id1, "spectrumid2":id2, 'Tanimoto_Similarity':sim},)
    if output_dir is not None: out.to_csv(output_dir, index_label=False)
    return out