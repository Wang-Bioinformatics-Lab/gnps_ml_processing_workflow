import re
import pandas as pd
from pyteomics.mgf import IndexedMGF
import datetime
import os
import argparse

year = None
month = None

def Bruker_Fragmentation_Prediction(summary_path:str, mgf_path:str, output_path:str):
    """This function follows the cleaning in 3DMolMS applied to Bruker qtof instruments.

    Args:
        summary_path (str): _description_
        mgf_path (str): _description_

    """
    summary = pd.read_csv(summary_path)
    mgf = IndexedMGF(mgf_path, index_by_scans=True)
    
    allowed_atoms = ['C', 'H', 'O', 'N', 'F', 'S', 'Cl', 'P', 'B', 'Br', 'I']
    print("Starting Size:", len(summary))
    # Remove structure-less entries. select instrument = qTof and Adduct in ['M+H','M-H']
    reduced_df = summary[~ (summary['Smiles'].isna() ) & (summary['msMassAnalyzer'] == 'qtof') & ((summary['Adduct'] == 'M+H') | (summary['Adduct'] == 'M-H'))].copy(deep=True)
    print("Lost {} structures when requiring smiles, msMassAnalyzer == 'qtof', and M+H/M-H".format(len(summary) - len(reduced_df)))
    # Remove all entires with atoms not in ['C', 'H', 'O', 'N', 'F', 'S', 'Cl', 'P', 'B', 'Br', 'I']
    reduced_df['Smiles_letters_only'] = reduced_df['Smiles'].apply(lambda x: "".join(re.findall("[a-zA-Z]+", x)))
    reduced_df['Smiles_cleaned'] = reduced_df['Smiles_letters_only'].apply(lambda x: "".join(re.findall("^[" + "|".join(allowed_atoms) + "]+$", x)))
    reduced_df = reduced_df[reduced_df['Smiles_cleaned'] != ""]
    reduced_df.drop(['Smiles_letters_only','Smiles_cleaned'], inplace=True, axis=1)
    print(reduced_df.msManufacturer.value_counts())
    reduced_df = reduced_df[reduced_df.msManufacturer == 'Bruker Daltonics']  
    print("Ending Size:", len(reduced_df))
    
    reduced_df.to_csv(output_path + '/Bruker_Fragmentation_Prediction_{}_{}.csv'.format(month, year), index=False)
    
    assert all([int(mgf[i-1]['params']['scans']) == i+1 for i in reduced_df.scan])
    spectra = [mgf[i+1] for i in reduced_df.scan]
    mgf.write(spectra,output_path + '/Bruker_Fragmentation_Prediction_{}_{}.mgf'.format(month, year), file_mode='w')
    
    
def main():
    subsets = ['Bruker_Fragmentation_Prediction']
    parser = argparse.ArgumentParser(
                    prog = 'GNPS2 Subset Generator',
                    description = 'This program generates predetermined subsets splits from GNPS2.')
    parser.add_argument('subset', choices=subsets)
    parser.add_argument('-s', '--split', action='store_true')
    args = parser.parse_args()
    
    now = datetime.datetime.now()
    global year
    global month
    year = now.year
    month = now.month

    csv_path = "./GNPS_ml_exports/ALL_GNPS_cleaned_{}_{}.csv".format(month, year)
    mgf_path = "./GNPS_ml_exports/ALL_GNPS_cleaned_{}_{}.mgf".format(month, year)
    
    output_path = './nf_output'
    if not os.path.isdir(output_path): os.makedirs(output_path)
    
    if args.split:
        raise NotImplementedError
    
    if args.subset == 'Bruker_Fragmentation_Prediction':
        Bruker_Fragmentation_Prediction(csv_path, mgf_path, output_path)
        
            
if __name__ == '__main__':
    main()