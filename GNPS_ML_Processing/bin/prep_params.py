import numpy as np
import argparse
import requests
import os
import json

def main():
    # We only want to generate these files quarter, so we'll check if it has already been done
    # if not os.path.isfile(final_csv_path):
    #     if not os.path.isfile(final_mgf_path):
    parser = argparse.ArgumentParser(description='Create Parallel Parameters.')
    # parser.add_argument('input_libraryname', default='ALL_GNPS')
    parser.add_argument('-s', '--structures_required', help="remove entries that don't include structures", action="store_true")
    parser.add_argument('-p', type=int, help='Parallelism', default=100)
    parser.add_argument('--all_GNPS_json_path', help='Path to all_GNPS.json', default="None")   # Using none in quotes for easy nextflow arguments
    args = parser.parse_args()

    # If a path to the JSON file is provided, use that instead of downloading from GNPS
    if args.all_GNPS_json_path != "None":
        all_spectra_list = json.load(open(args.all_GNPS_json_path, 'r'))
    else:
        gnps_url = "https://external.gnps2.org/gnpslibrary/{}.json".format("ALL_GNPS_NO_PROPOGATED")
        all_spectra_list = requests.get(gnps_url).json()

    if not os.path.isdir('./params'): os.makedirs('./params',exist_ok=True)

    if args.structures_required:
        org_len = len(all_spectra_list)
        all_spectra_list = [spectrum for spectrum in all_spectra_list if spectrum['Smiles'] != 'n/a' and spectrum['Smiles'] != 'n\/a']
        print("Found {} entries with structures out of {} structures: {:4.2f}%".format(len(all_spectra_list), org_len, len(all_spectra_list)/org_len*100))
    
    processors_numbers = args.p
    print("Recommended Parallelism:", max(1, int(len(all_spectra_list)/1000)))
    num_sections = max(1, processors_numbers)
    indices = np.array_split(np.arange(1,len(all_spectra_list)+1), num_sections)
    scan_start = [x[0] for x in indices]
    splits = np.array_split(all_spectra_list, num_sections)
    
    for section_idx in range(num_sections):
        # Save file name is params_splitNum_startScan.npy
        np.save('params/params_{}_{}.npy'.format(section_idx, scan_start[section_idx]), splits[section_idx])
        
if __name__ == '__main__':
    main()