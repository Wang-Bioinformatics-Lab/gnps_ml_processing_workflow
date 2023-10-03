import os
import argparse
import csv
from glob import glob
from pathlib import Path
import re
def merge(add_headers=False):
    """A simple script to combine csv and mgf files for our pipleline.
    """
    merged_csv_path = "ALL_MassBank_merged.csv"
    merged_mgf_path = "ALL_MassBank_merged.mgf"
    
    if add_headers:
        column_names = ['scan', 'spectrum_id','collision_energy','retention_time','Adduct','Compound_Source','Compund_Name','Precursor_MZ','ExactMass','Charge','Ion_Mode','Smiles','INCHI','InChIKey_smiles','InChIKey_inchi','msModel','msManufacturer','msDetector','msMassAnalyzer','msIonisation','msDissociationMethod','GNPS_library_membership','GNPS_Inst']
    
    if not os.path.isfile(merged_csv_path):
        if not os.path.isfile(merged_mgf_path):
            # This will sort based on the name. It will not be in any particular order,
            # but it will enforce mgf/csv pairs to be in the same order
            file_pattern = re.compile(r'^output_(.*)\.\w+')
            def get_order(file,):
                try:
                    match = file_pattern.match(Path(file).name)
                    return match.groups()[0]
                except Exception as e:
                        print(f"Failing on {file}")
                        raise e
            
            sorted_csv_files = sorted(glob('./output_*.csv'), key=get_order)
            sorted_mgf_files = sorted(glob('./output_*.mgf'), key=get_order)
            
            if add_headers:
                with open('temp_column_names.csv', 'w') as f:
                    csv.DictWriter(f, column_names).writeheader()
                sorted_csv_files = ['temp_column_names.csv'] + sorted_csv_files

            os.system("cat " + " ".join(sorted_csv_files) +"> " + merged_csv_path)
            os.system("cat " + " ".join(sorted_mgf_files) +"> " + merged_mgf_path)
            # os.system("rm " + " ".join(sorted_csv_files))
            # os.system("rm " + " ".join(sorted_mgf_files))
            
            if add_headers:
                try:
                    if os.path.isfile('temp_column_names.csv'):
                        os.remove('temp_column_names.csv')
                except Exception as e:
                    pass

def main():
    parser = argparse.ArgumentParser(description='Merge files for MassBank processing')
    parser.add_argument('--add_headers', action='store_true', help='Add headers to the merged file')
    args = parser.parse_args()
    merge(args.add_headers)

    
if __name__ == '__main__':
    main()