import re
import pandas as pd
import argparse
import os
import vaex
import pickle

    
def Bruker_Fragmentation_Prediction(summary_path:str, parquet_path:str):
    """This function follows the cleaning in 3DMolMS applied to Bruker qtof instruments.

    Args:
        summary_path (str): _description_
        parquet_path (str): _description_

    """
    reduced_df = pd.read_csv(summary_path)

    allowed_atoms = ['C', 'H', 'O', 'N', 'F', 'S', 'Cl', 'P', 'B', 'Br', 'I']

    # Remove structure-less entries. select instrument = qTof and Adduct in ['M+H','M-H']
    reduced_df = reduced_df[~ (reduced_df['Smiles'].isna() ) & (reduced_df['msMassAnalyzer'] == 'qtof') & ((reduced_df['Adduct'] == 'M+H') | (reduced_df['Adduct'] == 'M-H'))].copy(deep=True)

    # Remove all entires with atoms not in ['C', 'H', 'O', 'N', 'F', 'S', 'Cl', 'P', 'B', 'Br', 'I']
    reduced_df['Smiles_letters_only'] = reduced_df['Smiles'].apply(lambda x: "".join(re.findall("[a-zA-Z]+", x)))
    reduced_df['Smiles_cleaned'] = reduced_df['Smiles_letters_only'].apply(lambda x: "".join(re.findall("^[" + "|".join(allowed_atoms) + "]+$", x)))
    reduced_df = reduced_df[reduced_df['Smiles_cleaned'] != ""]
    reduced_df.drop(['Smiles_letters_only','Smiles_cleaned'], inplace=True, axis=1)

    reduced_df = reduced_df[reduced_df.msManufacturer == 'Bruker Daltonics']  
    reduced_df.to_csv('./summary/Bruker_Fragmentation_Prediction.csv', index=False)
    
    id_list = list(reduced_df.spectrum_id )
    del reduced_df
    
    parquet_as_df = vaex.open(parquet_path)
    parquet_as_df = parquet_as_df[parquet_as_df.spectrum_id.isin(id_list)]
    parquet_as_df.export_parquet('./spectra/Bruker_Fragmentation_Prediction.parquet')
    
def MH_MNA_Translation(summary_path:str, parquet_path:str):

    reduced_df = pd.read_csv(summary_path)

    reduced_df = reduced_df.loc[reduced_df.msMassAnalyzer == 'orbitrap']
    reduced_df = reduced_df.loc[reduced_df.GNPS_Inst == 'orbitrap']
    reduced_df = reduced_df.loc[~reduced_df.Smiles.isna()]
    reduced_df = reduced_df.loc[(reduced_df.Adduct == 'M+H') | (reduced_df.Adduct == 'M+Na')]

    reduced_df.to_csv('./summary/MH_MNA_Translation.csv', index=False)
    
    def Generate_Pairs_List(summary:pd.DataFrame):
        output = []
        for _, row in summary[['spectrum_id', 'Smiles','Adduct']].iterrows():
            similar_ids = summary.loc[(row["Smiles"] == summary["Smiles"]) & (row["Adduct"] != summary["Adduct"]),"spectrum_id"].values
            output.append((row["spectrum_id"],list(similar_ids)))
        return output
    
    pairs_list = Generate_Pairs_List(reduced_df)
    with open("./util/MH_MNA_Translation_pairs.pkl", "wb") as fp:
        pickle.dump(pairs_list, fp)
        
    id_list = list(reduced_df.spectrum_id )
    del reduced_df
    del pairs_list
    
    parquet_as_df = vaex.open(parquet_path)
    parquet_as_df = parquet_as_df[parquet_as_df.spectrum_id.isin(id_list)]
    parquet_as_df.export_parquet('./spectra/MH_MNA_Translation.parquet')
    
def main():
    subsets = ['Bruker_Fragmentation_Prediction','MH_MNA_Translation','GNPS_default']
    parser = argparse.ArgumentParser(
                    prog = 'GNPS2 Subset Generator',
                    description = 'This program generates predetermined subsets splits from GNPS2.')
    parser.add_argument('subset', choices=subsets)
    parser.add_argument('-s', '--split', action='store_true')
    args = parser.parse_args()
    
    csv_path     = "ALL_GNPS_cleaned.csv"
    parquet_path = "ALL_GNPS_cleaned.parquet"
    
    if not os.path.isdir('./spectra'): os.makedirs('./spectra', exist_ok=True)
    if not os.path.isdir('./summary'): os.makedirs('./summary', exist_ok=True)
    if not os.path.isdir('./util'): os.makedirs('./util', exist_ok=True)
    
    
    if args.subset == 'Bruker_Fragmentation_Prediction':
        Bruker_Fragmentation_Prediction(csv_path, parquet_path)
    elif args.subset == 'MH_MNA_Translation':
        MH_MNA_Translation(csv_path, parquet_path)
    elif args.subset == 'GNPS_default':
        Bruker_Fragmentation_Prediction(csv_path, parquet_path)
        MH_MNA_Translation(csv_path, parquet_path)
        
            
if __name__ == '__main__':
    main()