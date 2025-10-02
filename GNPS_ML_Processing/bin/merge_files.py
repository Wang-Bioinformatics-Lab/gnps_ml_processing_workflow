import os
import argparse
from glob import glob
from pathlib import Path
import re

def merge(include_massbank=False, include_riken=False, include_mona=False, include_mgf_exports=False):
    """A simple script to combine csv and mgf files for our pipleline.
    """
    merged_csv_path = "ALL_GNPS_merged.csv"
    merged_mgf_path = "ALL_GNPS_merged.mgf"
    
    # The reason that we have these seperate flags is that we don't control the naming conventions of these files as tightly.
    # Therefore, we make sure they remain in the same order, based on their sources
    if include_massbank:
        massbank_csv = ['ALL_MassBank_merged.csv']
        massbank_mgf = ['ALL_MassBank_merged.mgf']
    else:
        massbank_csv = []
        massbank_mgf = []

    if include_riken:
        riken_csv = list(glob('./RIKEN-*.csv'))
        riken_mgf = list(glob('./RIKEN-*.mgf'))
        # Sort so we know the order is the same
        riken_csv.sort()
        riken_mgf.sort()
    else:
        riken_csv = []
        riken_mgf = []

    if include_mona:
        mona_csv = list(glob('./MONA-*.csv'))
        mona_mgf = list(glob('./MONA-*.mgf'))
        # Sort so we know the order is the same
        mona_csv.sort()
        mona_mgf.sort()
    else:
        mona_csv = []
        mona_mgf = []

    if include_mgf_exports:
        mgf_exports_csv = list(glob('MGF-Export-*.csv'))
        mgf_exports_mgf = list(glob('MGF-Export-*.mgf'))
        mgf_exports_csv.sort()
        mgf_exports_mgf.sort()
    else:
        mgf_exports_csv = []
        mgf_exports_mgf = []
    
    if not os.path.isfile(merged_csv_path):
        if not os.path.isfile(merged_mgf_path):
            file_pattern = re.compile(r'.*?(\d+).*?')
            def get_order(file,):
                match = file_pattern.match(Path(file).name)
                return int(match.groups()[0])

            sorted_csv_files = sorted(glob('./temp_*.csv'), key=get_order) + massbank_csv + riken_csv + mona_csv + mgf_exports_csv
            sorted_mgf_files = sorted(glob('./temp_*.mgf'), key=get_order) + massbank_mgf + riken_mgf + mona_mgf + mgf_exports_mgf

            print(f"Concatenating {len(sorted_csv_files)} csv files")
            print(f"Concatenating {len(sorted_mgf_files)} mgf files")

            # os.system("cat " + " ".join(sorted_csv_files) +"> " + merged_csv_path)
            # os.system("cat " + " ".join(sorted_mgf_files) +"> " + merged_mgf_path)
            for file in sorted_csv_files:
                os.system(f"cat {file} >> {merged_csv_path}")
                
            for file in sorted_mgf_files:
                os.system(f"cat {file} >> {merged_mgf_path}")

def main():
    parser = argparse.ArgumentParser(description='Merge csv and mgf files.')
    parser.add_argument('--include_massbank', help='Include MassBank files.', action='store_true')
    parser.add_argument('--include_riken', help='Include Riken files.', action='store_true')
    parser.add_argument('--include_mona', help='Include MONA files.', action='store_true')
    parser.add_argument('--include_mgf_exports', help='Include MGF exports.', action='store_true')
    args = parser.parse_args()

    merge(args.include_massbank, args.include_riken, args.include_mona, args.include_mgf_exports)

    
    
if __name__ == '__main__':
    main()