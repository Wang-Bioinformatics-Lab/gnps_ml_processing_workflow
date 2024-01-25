import pandas as pd
import argparse
import shutil

def main():
    parser = argparse.ArgumentParser(description='Spoof matchms caching')
    parser.add_argument('--input_csv_path', type=str, help='Path to input csv file')
    parser.add_argument('--cached_compound_name_annotation_path', type=str, help='Path to cached compound name annotation file')
    args = parser.parse_args()
    
    titles = pd.read_csv(args.input_csv_path)['spectrum_id']
    
    output_name = './compound_name_annotation.csv'
    
    shutil.copy(args.cached_compound_name_annotation_path, output_name)
    
    # Append to cached_compound_name_annotation_path
    with open(output_name, 'a') as f:
        with open(output_name, 'r') as f2:
            # Ensure we're on a newline
            if f2.read()[-1] != '\n':
                f.write('\n')
        
        # CSV headers: compound_name,smiles,inchi,inchikey,monoisotopic_mass
        for title in titles:
            if title.startswith("CCMSLIB") or title.startswith("MSBNK"):
                f.write(f"{title},,,,\n")
    
  
if __name__ == "__main__":
    main()