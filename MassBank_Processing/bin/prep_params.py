import argparse
import os
from glob import glob

import numpy as np

def main():
    parser = argparse.ArgumentParser(description='Create Parallel Parameters.')
    parser.add_argument('--glob', help="Glob to discover files.", required=True)
    parser.add_argument('--blacklist', help="Text file containing blacklisted file names.", required=False, default=None)
    parser.add_argument('-p', type=int, help='Parallelism', default=100)

    args = parser.parse_args()
    
    processors_numbers = args.p
    
    if not os.path.isdir('./params'):
        os.makedirs('./params',exist_ok=True)
        
    files_list = sorted(glob(args.glob))    # Ensures ordering is consistent across file system
    
    if args.blacklist is not None:
        # Open black list and see if it's a blacklisted file
        with open(args.blacklist, 'r') as f:
            blacklist = f.read().splitlines()
            # Check if the filename is a substring of any of the blacklisted files
            temp = []
            for filename in files_list:
                blacklisted = False
                for x in blacklist:
                    if os.path.basename(filename) in x.split('/', 1)[1]:
                        blacklisted = True
                        break
                if not blacklisted:
                    temp.append(filename)
                        
            print(f"Removed {len(files_list) - len(temp)} blacklisted files.")
            print(f"Remaining files: {len(temp)}")
            files_list = temp

   
    num_sections = min(max(1, processors_numbers), len(files_list))
    indices = np.array_split(np.arange(1,len(files_list)+1), num_sections)
    scan_start = [x[0] for x in indices]
    splits = np.array_split(files_list, num_sections)
    
    for section_idx in range(num_sections):
        # Save file name is params_splitNum_startScan.npy
        np.save(f'params/params_{section_idx}_{scan_start[section_idx]}.npy', splits[section_idx])
        
if __name__ == '__main__':
    main()
    