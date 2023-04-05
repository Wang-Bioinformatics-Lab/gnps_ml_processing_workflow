import os
import argparse
from glob import glob
from pathlib import Path
import re

def main():
    """A simple script to combine csv and mgf files for our pipleline.
    """
    merged_csv_path = "ALL_GNPS_merged.csv"
    merged_mgf_path = "ALL_GNPS_merged.mgf"
    
    if not os.path.isfile(merged_csv_path):
        if not os.path.isfile(merged_mgf_path):
            file_pattern = re.compile(r'.*?(\d+).*?')
            def get_order(file,):
                match = file_pattern.match(Path(file).name)
                return int(match.groups()[0])

            sorted_csv_files = sorted(glob('./temp_*.csv'), key=get_order)
            sorted_mgf_files = sorted(glob('./temp_*.mgf'), key=get_order)

            os.system("cat " + " ".join(sorted_csv_files) +"> " + merged_csv_path)
            os.system("cat " + " ".join(sorted_mgf_files) +"> " + merged_mgf_path)
            # os.system("rm " + " ".join(sorted_csv_files))
            # os.system("rm " + " ".join(sorted_mgf_files))

    
if __name__ == '__main__':
    main()