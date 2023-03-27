import re
import pandas as pd
import argparse
import os
import vaex
import pickle
from utils import build_tanimoto_similarity_list_precomputed

    
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
    
def Orbitrap_Fragmentation_Prediction(summary_path:str, parquet_path:str):  
    """This function follows the cleaning in 3DMolMS applied to orbitrap instruments.

    Args:
        summary_path (str): _description_
        parquet_path (str): _description_

    """
    reduced_df = pd.read_csv(summary_path)

    allowed_atoms = ['C', 'H', 'O', 'N', 'F', 'S', 'Cl', 'P', 'B', 'Br', 'I']

    # Remove structure-less entries. select instrument = orbitrap, GNPS_Inst = orbitrap, adduct = M+H
    reduced_df = reduced_df[(~reduced_df['Smiles'].isna()) & (reduced_df['msMassAnalyzer'] == 'orbitrap') & (reduced_df['GNPS_Inst'] == 'orbitrap') & (reduced_df['Adduct'] == 'M+H') ]
    
    # Remove all entires with atoms not in ['C', 'H', 'O', 'N', 'F', 'S', 'Cl', 'P', 'B', 'Br', 'I']
    reduced_df['Smiles_letters_only'] = reduced_df['Smiles'].apply(lambda x: "".join(re.findall("[a-zA-Z]+", x)))
    reduced_df['Smiles_cleaned'] = reduced_df['Smiles_letters_only'].apply(lambda x: "".join(re.findall("^[" + "|".join(allowed_atoms) + "]+$", x)))
    reduced_df = reduced_df[reduced_df['Smiles_cleaned'] != ""]
    reduced_df.drop(['Smiles_letters_only','Smiles_cleaned'], inplace=True, axis=1)

    # reduced_df = reduced_df[reduced_df.msManufacturer == 'Bruker Daltonics']  
    reduced_df.to_csv('./summary/Orbitrap_Fragmentation_Prediction.csv', index=False)
    
    id_list = list(reduced_df.spectrum_id )
    del reduced_df
    
    parquet_as_df = vaex.open(parquet_path)
    parquet_as_df = parquet_as_df[parquet_as_df.spectrum_id.isin(id_list)]
    parquet_as_df.export_parquet('./spectra/Orbitrap_Fragmentation_Prediction.parquet')

def Thermo_Bruker_Translation(summary_path:str, parquet_path:str):
    """This function generates a csv, parquet file, and pairs list that contains the spectra of the same compounds from Thermo and Bruker instruments.
    More specifically, we include Q Exactive Instruments from Thermo Fisher and maXis impact instruments from Bruker Daltonics.
    In addition, we require that the adducts are identical. There are no conditions on the collision energy.

    Args:
        summary_path (str): _description_
        parquet_path (str): _description_
    """
    df = pd.read_csv(summary_path)
    df = df.loc[~df.Smiles.isna()]
    df = df.loc[df.Adduct.isin(['M+H','M-H'])]
    df.Smiles = df.Smiles.astype(str)
    thermo = df.loc[(df.msManufacturer=="Thermo") & (df.msModel == "Q Exactive"),['spectrum_id', 'Smiles', 'Adduct']]
    bruker = df.loc[(df.msManufacturer=="Bruker Daltonics") * (df.msModel == "maXis impact"),['spectrum_id', 'Smiles', 'Adduct']]
    spectrum_ids = bruker.merge(thermo, on=['Smiles','Adduct'], how='inner')
    
    # Generate a list of spectrum pairs
    pairs_list = []
    for id in spectrum_ids.spectrum_id_x.unique():
        pairs_list.append((id, list(spectrum_ids.loc[spectrum_ids.spectrum_id_x == id, "spectrum_id_y"].values)))
    with open("./util/Thermo_Bruker_Translation_pairs.pkl", "wb") as fp:
        pickle.dump(pairs_list, fp)
    
    spectrum_ids = pd.concat((spectrum_ids.spectrum_id_x, spectrum_ids.spectrum_id_y)).values
    
    # Select dataframe entried where spectrum_id is in spectrum_ids
    df = df[df.spectrum_id.isin(spectrum_ids)]
    
    df.to_csv('./summary/Thermo_Bruker_Translation.csv', index=False)
    
    id_list = list(df.spectrum_id )
    del df, bruker, thermo
    
    parquet_as_df = vaex.open(parquet_path)
    parquet_as_df = parquet_as_df[parquet_as_df.spectrum_id.isin(id_list)]
    parquet_as_df.export_parquet('./spectra/Thermo_Bruker_Translation.parquet')
    
def Structural_Modification(summary_path:str, parquet_path:str):
    # This dataset: Struct1 + Spec1 + struct2 -> spec2
    
    if not os.path.isdir('./spectra/Structural_Modification'):
        os.makedirs('./spectra/Structural_Modification', exist_ok=True)
    if not os.path.isdir('./summary/Structural_Modification'):
        os.makedirs('./summary/Structural_Modification', exist_ok=True)
    if not os.path.isdir('./util/Structural_Modification'):
        os.makedirs('./util/Structural_Modification', exist_ok=True)
    
    df = pd.read_csv(summary_path)
    df = df.loc[~df.Smiles.isna()]
    
    parquet_as_df = vaex.open(parquet_path)
        
    # Orbitrap Portion
    positive = df.loc[(df.Adduct == 'M+H') & (df.msMassAnalyzer == 'orbitrap') & (df.GNPS_Inst == 'orbitrap')]
    negative = df.loc[(df.Adduct == 'M-H') & (df.msMassAnalyzer == 'orbitrap') & (df.GNPS_Inst == 'orbitrap')]
    # Remove entries with duplicate SMILES
    positive = positive.drop_duplicates(subset=['Smiles'])
    negative = negative.drop_duplicates(subset=['Smiles'])
    # Calculate structural similarities
    positive_sim = build_tanimoto_similarity_list_precomputed(positive, similarity_threshold = 0.70)
    negative_sim = build_tanimoto_similarity_list_precomputed(negative, similarity_threshold = 0.70)
    # Remove entires not included in positive_sim
    positive_ids = pd.concat((positive_sim.spectrumid1,positive_sim.spectrumid2)).unique()
    negative_ids = pd.concat((negative_sim.spectrumid1,negative_sim.spectrumid2)).unique()
    positive = positive[positive.spectrum_id.isin(positive_ids)]
    negative = negative[negative.spectrum_id.isin(negative_ids)]
    # Save to csv
    positive.to_csv('./summary/Structural_Modification/orbitrap_positive.csv', index=False)
    negative.to_csv('./summary/Structural_Modification/orbitrap_negative.csv', index=False)
    # Save to parquet
    parquet_as_df[parquet_as_df.spectrum_id.isin(positive_ids)].export_parquet('./spectra/Structural_Modification/orbitrap_positive.parquet')
    parquet_as_df[parquet_as_df.spectrum_id.isin(negative_ids)].export_parquet('./spectra/Structural_Modification/orbitrap_negative.parquet')
    
    # Save Similarities to utils
    positive_sim.to_csv('./util/Structural_Modification/orbitrap_positive_sim.csv', index=False)
    negative_sim.to_csv('./util/Structural_Modification/orbitrap_negative_sim.csv', index=False)
    
    # qtof Portion
    positive = df.loc[(df.Adduct == 'M+H') & (df.msMassAnalyzer == 'qtof') & (df.GNPS_Inst == 'qtof')]
    negative = df.loc[(df.Adduct == 'M-H') & (df.msMassAnalyzer == 'qtof') & (df.GNPS_Inst == 'qtof')]
    # Remove entries with duplicate SMILES
    positive = positive.drop_duplicates(subset=['Smiles'])
    negative = negative.drop_duplicates(subset=['Smiles'])
    # Calculate structural similarities
    positive_sim = build_tanimoto_similarity_list_precomputed(positive, similarity_threshold = 0.70)
    negative_sim = build_tanimoto_similarity_list_precomputed(negative, similarity_threshold = 0.70)
    # Remove entires not included in positive_sim
    positive_ids = pd.concat((positive_sim.spectrumid1,positive_sim.spectrumid2)).unique()
    negative_ids = pd.concat((negative_sim.spectrumid1,negative_sim.spectrumid2)).unique()
    positive = positive[positive.spectrum_id.isin(positive_ids)]
    negative = negative[negative.spectrum_id.isin(negative_ids)]
    # Save to csv
    positive.to_csv('./summary/Structural_Modification/qtof_positive.csv', index=False)
    negative.to_csv('./summary/Structural_Modification/qtof_negative.csv', index=False)
    # Save to parquet
    parquet_as_df[parquet_as_df.spectrum_id.isin(positive_ids)].export_parquet('./spectra/Structural_Modification/qtof_positive.parquet')
    parquet_as_df[parquet_as_df.spectrum_id.isin(negative_ids)].export_parquet('./spectra/Structural_Modification/qtof_negative.parquet')
    
    # Save Similarities to utils
    positive_sim.to_csv('./util/Structural_Modification/qtof_positive_sim.csv', index=False)
    negative_sim.to_csv('./util/Structural_Modification/qtof_negative_sim.csv', index=False)
    

def main():
    subsets = ['Bruker_Fragmentation_Prediction','MH_MNA_Translation','Orbitrap_Fragmentation_Prediction',\
                'Thermo_Bruker_Translation','Structural_Modification','GNPS_default']
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
    elif args.subset == 'Orbitrap_Fragmentation_Prediction':
        Orbitrap_Fragmentation_Prediction(csv_path, parquet_path)
    elif args.subset == 'Thermo_Bruker_Translation':
        Thermo_Bruker_Translation(csv_path, parquet_path)
    elif args.subset == 'Structural_Modification':
        Structural_Modification(csv_path, parquet_path)
    elif args.subset == 'GNPS_default':
        Bruker_Fragmentation_Prediction(csv_path, parquet_path)
        MH_MNA_Translation(csv_path, parquet_path)
        Orbitrap_Fragmentation_Prediction(csv_path, parquet_path)
        Thermo_Bruker_Translation(csv_path, parquet_path)
        Structural_Modification(csv_path, parquet_path)
        
            
if __name__ == '__main__':
    main()