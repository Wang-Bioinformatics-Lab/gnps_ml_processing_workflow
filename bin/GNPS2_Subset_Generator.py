import re
import pandas as pd
from pyteomics.mgf import IndexedMGF
import datetime
import os
import argparse

year = None
quarter = None

def sync_df_and_parquet(df:pd.DataFrame, parquet:pd.DataFrame):
    return parquet.loc[[x in df.spectrum_id for x in parquet.index]]

def Bruker_Fragmentation_Prediction(summary_path:str, parquet_path:str, output_path:str):
    """This function follows the cleaning in 3DMolMS applied to Bruker qtof instruments.

    Args:
        summary_path (str): _description_
        parquet_path (str): _description_

    """
    summary = pd.read_csv(summary_path)
    parquet_as_df = pd.read_parquet(parquet_path)
    
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
    
    reduced_df.to_csv(output_path + '/Bruker_Fragmentation_Prediction_{}_{}.csv'.format(quarter, year), index=False)
       
    parquet_as_df = sync_df_and_parquet(reduced_df, parquet_as_df)
    parquet_as_df.to_parquet(output_path + '/Bruker_Fragmentation_Prediction_{}_{}.mgf'.format(quarter, year))
    
def MH_MNA_Translation(summary_path:str, parquet_path:str, output_path:str):
    reduced_df = pd.read_csv(summary_path)
    parquet_as_df = pd.read_parquet(parquet_path)
    
    reduced_df = reduced_df.loc[(reduced_df.Adduct == 'M+H') | (reduced_df.Adduct == 'M+NA')]
    reduced_df.to_csv(output_path + '/Bruker_Fragmentation_Prediction_{}_{}.csv'.format(quarter, year), index=False)
       
    parquet_as_df = sync_df_and_parquet(reduced_df, parquet_as_df)
    parquet_as_df.to_parquet(output_path + '/Bruker_Fragmentation_Prediction_{}_{}.mgf'.format(quarter, year))
    
    
    raise NotImplementedError
    
def main():
    subsets = ['Bruker_Fragmentation_Prediction','MH_MNA_Translation']
    parser = argparse.ArgumentParser(
                    prog = 'GNPS2 Subset Generator',
                    description = 'This program generates predetermined subsets splits from GNPS2.')
    parser.add_argument('subset', choices=subsets)
    parser.add_argument('-s', '--split', action='store_true')
    args = parser.parse_args()
    
    now = datetime.datetime.now()
    global year
    global quarter
    year = now.year
    quarter = int(now.month/4) + 1

    csv_path     = "./GNPS_ml_exports/ALL_GNPS_cleaned_{}_{}.csv".format(quarter, year)
    parquet_path = "./GNPS_ml_exports/ALL_GNPS_cleaned_{}_{}.parquet".format(quarter, year)
    
    output_path = './nf_output'
    if not os.path.isdir(output_path): os.makedirs(output_path)
    
    if args.split:
        raise NotImplementedError
    
    if args.subset == 'Bruker_Fragmentation_Prediction':
        Bruker_Fragmentation_Prediction(csv_path, parquet_path, output_path)
    elif args.subset == 'MH_MNA_Translation':
        MH_MNA_Translation(csv_path, parquet_path, output_path)
        
            
if __name__ == '__main__':
    main()