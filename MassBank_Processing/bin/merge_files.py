import os
import argparse
import csv
from glob import glob
from pathlib import Path
import re
import tempfile
import pytest

def concatenate_files_in_batches(file_list, output_path, batch_size):
    for i in range(0, len(file_list), batch_size):
        batch = file_list[i:i+batch_size]
        os.system("cat " + " ".join(batch) + ">> " + output_path)

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

            # Merge files
            batch_size = 100
            concatenate_files_in_batches(sorted_csv_files, merged_csv_path, batch_size)
            concatenate_files_in_batches(sorted_mgf_files, merged_mgf_path, batch_size)
            
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
    
    
@pytest.fixture
def temp_output_file():
    # Create a temporary output file for testing
    with tempfile.NamedTemporaryFile(delete=False) as temp_file:
        yield temp_file.name

def test_concatenate_files_in_batches(temp_output_file):
    # List of files to concatenate
    file_list = ["file1.txt", "file2.txt", "file3.txt", "file4.txt"]

    # Create empty files for testing
    for file in file_list:
        with open(file, "w") as f:
            # Write the file name to the file
            f.write(file)

    # Concatenate files in batches of 2
    concatenate_files_in_batches(file_list, temp_output_file, 2)

    # Verify the content of the output file
    with open(temp_output_file, "r") as output_file:
        result = output_file.read()

    expected_result = "file1.txtfile2.txtfile3.txtfile4.txt"
    assert result == expected_result

    # Clean up temporary files
    for file in file_list:
        os.remove(file)
    os.remove(temp_output_file)