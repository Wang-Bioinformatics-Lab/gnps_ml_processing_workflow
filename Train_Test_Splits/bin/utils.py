import numpy as np
import pandas as pd
from tqdm import tqdm

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

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

def knapsack_subset(counts, num_test_points):
    n = len(counts)
    dp = [[0] * (num_test_points + 1) for _ in range(n + 1)]
    
    for i in tqdm(range(1, n + 1)):
        for j in range(1, num_test_points + 1):
            if counts.iat[i - 1] <= j:
                dp[i][j] = max(counts[i - 1] + dp[i - 1][j - counts.iat[i - 1]], dp[i - 1][j])
            else:
                dp[i][j] = dp[i - 1][j]
    
    selected = []
    i, j = n, num_test_points
    while i > 0 and j > 0:
        if dp[i][j] != dp[i - 1][j]:
            selected.append(counts.iloc[[i - 1]].index.values.item())
            j -= counts.iat[i - 1]
        i -= 1
    
    return selected

def greedy_subset(counts, num_test_points):
    """ Adds the smallest set to the test set until num_test_points is exceeeded"""
    
    selected = []
    running_sum = 0
    counts = counts.sort_values()
    for i in tqdm(range(len(counts))):
        if running_sum < num_test_points:
            selected.append(counts.iloc[[i]].index.values.item())
            running_sum += counts.iat[i]
        else:
            break
    return selected

def calculate_murcko_histogram(smiles):
    """
    Calculate the Murcko histogram for a given SMILES string as defined in Algorithm 3 of
    "Emergence of molecular structures from self-supervised learning on mass spectra" by 
    Bushuiev et al.

    Args:
        smiles (str): The SMILES string representing the molecule.

    Raises:
        ValueError: If the SMILES string cannot be parsed.

    Returns:
        dict: A dictionary representing the Murcko histogram. The keys are tuples (r, l) where r is the number of adjacent rings and l is the number of adjacent linkers. The values are the counts of molecules with the corresponding (r, l) values.
    """    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Failed to parse smiles: {smiles}")
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    
    # Get set of rings
    ring_info = scaffold.GetRingInfo()
    rings_by_atom_index = ring_info.AtomRings()
    # Remove rings if their atom count <= 3
    rings_by_atom_index = [ring for ring in rings_by_atom_index if len(ring) > 3]
    # Get atoms with a degree > 1 and not in any ring
    non_ring_atoms = [atom.GetIdx() for atom in scaffold.GetAtoms() if atom.GetDegree() > 1 and not any(atom.GetIdx() in ring for ring in rings_by_atom_index)]
    
    histogram = dict()
    
    for idx, ring in enumerate(rings_by_atom_index):
        other_rings = rings_by_atom_index[:idx] + rings_by_atom_index[idx+1:]
        # Get number of adjacent rings
        r = 0
        for other_ring in other_rings:
            if any(atom in other_ring for atom in ring):
                r += 1
                continue
        # Get number of adjacent linkers
        l = 0
        for atom_index in ring:
            for bonded_atom in scaffold.GetAtomWithIdx(atom_index).GetNeighbors():
                if bonded_atom.GetIdx() in non_ring_atoms:
                    l += 1
                    continue
        if histogram.get((r,l)) is None:
            histogram[r] = dict()
            histogram[r][l] = 1
        else:
            if histogram[r].get(l) is None:
                histogram[r][l] = 1
            else:
                histogram[r][l] += 1
            
    return histogram

def calculate_murcko_histogram_distance(histogram1, histogram2, min_rings=4):
    """
    Calculate the distance between two Murcko histograms as defined in Algorithm 4 of
    "Emergence of molecular structures from self-supervised learning on mass spectra" by
    Bushuiev et al.
    """
    
    # Check if we have enough rings to apply distance criteria
    histogram1_count = 0
    for r in histogram1.keys():
        histogram1_count += sum(list(histogram1[r].values()))
    histogram2_count = 0
    for r in histogram2.keys():
        histogram2_count += sum(list(histogram2[r].values()))
    if min(histogram1_count, histogram2_count) < min_rings:
        # Return the identity of the two dicts
        if sorted(list(histogram1.keys())) != sorted(list(histogram2.keys())):
            return np.inf
        for r in histogram1.keys():
            if histogram1[r] != histogram2[r]:
                return np.inf
            
        return 0
    
    # Calculate distance
    distance = 0
    for r in set(histogram1.keys()).union(histogram2.keys()):
        if histogram1.get(r) is None:
            hist1_sum = 0
        else:
            hist1_sum = sum(list(histogram1[r].values()))
        if histogram2.get(r) is None:
            hist2_sum = 0
        else:
            hist2_sum = sum(list(histogram2[r].values()))

        distance += abs(hist1_sum - hist2_sum)
        
    return distance