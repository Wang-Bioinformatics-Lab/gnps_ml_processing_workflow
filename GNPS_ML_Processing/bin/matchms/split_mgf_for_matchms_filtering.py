import argparse
import os

import matchms
import numpy as np
from time import time

def main():
    parser = argparse.ArgumentParser(description='Split MGF file into smaller chunks for parallel processing')
    parser.add_argument('--input_mgf_path', help='Input MGF file')
    parser.add_argument('--output_path', help='Folder to output spectra')
    parser.add_argument('--splits', help='Number of splits to make', default=1000, type=int)
    args = parser.parse_args()
    
    if not os.path.isdir(args.output_path):
        os.makedirs(args.output_path, exist_ok=True)
    
    # Load spectra
    print("Loading spectra...", flush=True)
    start_time = time()
    spectra = np.array(list(matchms.importing.load_from_mgf(args.input_mgf_path)))
    print("Loaded {} spectra in {:.2f} seconds".format(len(spectra), time() - start_time), flush=True)
    
    # Split spectra
    print(f"Splitting spectra into {args.splits} chunks...", flush=True)
    start_time = time()
    split_arrays = np.array_split(spectra, args.splits)
    print("Split {} spectra into {} chunks in {:.2f} seconds".format(len(spectra), args.splits, time() - start_time), flush=True)
    
    # Write arrays to mgf, note that matchms preserves metadata for scan nums
    print(f"Writing spectra to mgf in {args.output_path}...", flush=True)
    for i, split_array in enumerate(split_arrays):
        matchms.exporting.save_as_mgf(split_array.tolist(), args.output_path + "/split_" + str(i) + ".mgf")
    print("Wrote {} spectra to {} chunks in {:.2f} seconds".format(len(spectra), args.splits, time() - start_time), flush=True)
    
if __name__ == "__main__":
    main()