import argparse
import os
from glob import glob
import numpy as np
import pandas as pd
from pyteomics import mgf
from pyteomics.mgf import IndexedMGF
from time import time

def select_test(input_csv, input_mgf, num_test_points):
    """This function takes in an input_csv, an input_mgf, and the number of test points
        to be randomly generated. The function selects these test points and writes
        them to test_rows.csv, test_rows.mgf, train_rows.csv, and train_rows.mgf.

    Args:
        input_csv (_type_): _description_
        input_mgf (_type_): _description_
    """
    summary = pd.read_csv(input_csv)  
    
    spectra     = IndexedMGF(input_mgf, index_by_scans=False)
    
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
    
    # Select Test Points
    test_rows = summary.sample(n=num_test_points)
    print("Total number of test rows: ", len(test_rows), flush=True)
    train_rows = summary.loc[~ summary.scan.isin(test_rows.spectrum_id)]
    print("Total number of potential training rows: ", len(train_rows), flush=True)
    test_ids = test_rows['spectrum_id'].values
    train_ids = train_rows['spectrum_id'].values
    
    print("Writing Test Rows", flush=True)
    start_time = time()
    # Write Test Rows
    test_rows.to_csv('test_rows.csv', index=False)
    mgf.write(spectra[[str(spectrum_id) for spectrum_id in test_ids]],  output='test_rows.mgf')
    print(f"Done in {time() - start_time} seconds.", flush=True)
    
    print("Writing Train Rows", flush=True)
    start_time = time()
    # Write Train Rows
    train_rows.to_csv('train_rows.csv', index=False)
    mgf.write(spectra[[str(spectrum_id) for spectrum_id in train_ids]],  output='train_rows.mgf')
    print(f"Done in {time() - start_time} seconds.", flush=True)

def main():
    parser = argparse.ArgumentParser(description='Create Parallel Parameters.')
    parser.add_argument('--input_csv', type=str, help='Input CSV')
    parser.add_argument('--input_mgf', type=str, help='Input MGF')
    parser.add_argument('--num_test_points', help='Number of Test Points')
    args = parser.parse_args()
    
    select_test(args.input_csv, args.input_mgf, args.num_test_points)
    
if __name__ == '__main__':
    main()
    
    
def test_select_test():
    """This function tests the select_test function.
    """
    input_csv = '../data/test_data/test_summary.csv'
    input_mgf = '../data/test_data/test.mgf'
    num_test_points = 2
    
    select_test(input_csv, input_mgf, num_test_points)
    
    test_rows = pd.read_csv('test_rows.csv')
    train_rows = pd.read_csv('train_rows.csv')
    
    assert len(test_rows) == num_test_points
    assert len(train_rows) == 8
    
    os.remove('test_rows.csv')
    os.remove('test_rows.mgf')
    os.remove('train_rows.csv')
    os.remove('train_rows.mgf')