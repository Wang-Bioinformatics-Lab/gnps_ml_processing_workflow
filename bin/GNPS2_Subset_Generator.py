import re
import pandas as pd
import argparse
import os
import vaex
import pickle
from utils import build_tanimoto_similarity_list_precomputed, generate_fingerprints, harmonize_smiles_rdkit
import dask.dataframe as dd
    
def MH_MNA_Translation(summary_path:str, parquet_path:str):

    reduced_df = pd.read_csv(summary_path)

    # reduced_df = reduced_df.loc[reduced_df.msMassAnalyzer == 'orbitrap']
    # reduced_df = reduced_df.loc[reduced_df.GNPS_Inst == 'orbitrap']
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
    
def Fingerprint_Prediction(summary_path:str, parquet_path:str):
    reduced_df = pd.read_csv(summary_path)

    # reduced_df = reduced_df.loc[reduced_df.msMassAnalyzer == 'orbitrap']
    # reduced_df = reduced_df.loc[reduced_df.GNPS_Inst == 'orbitrap']
    reduced_df = reduced_df.loc[~reduced_df.Smiles.isna()]
    # reduced_df = reduced_df.loc[(reduced_df.Adduct == 'M+H')]
    reduced_df = generate_fingerprints(reduced_df)
    reduced_df.to_csv('./summary/Fingerprint_Prediction.csv', index=False)
    
    id_list = list(reduced_df.spectrum_id )
    del reduced_df
    
    parquet_as_df = vaex.open(parquet_path)
    parquet_as_df = parquet_as_df[parquet_as_df.spectrum_id.isin(id_list)]
    parquet_as_df.export_parquet('./spectra/Fingerprint_Prediction.parquet') 
    
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
    # df = df.loc[df.Adduct.isin(['M+H','M-H'])]
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
  
def Structural_Similarity_Prediction(summary_path:str, parquet_path:str):
    # This dataset: spec1 + spec2 -> struct similarity
    parquet_as_df = vaex.open(parquet_path)
            
    # df = df.loc[df.Adduct == 'M+H']
    # qtof = df.loc[(df.msMassAnalyzer == 'qtof') & (df.GNPS_Inst == 'qtof')]    
    # qtof = generate_fingerprints(qtof)
    # sim = build_tanimoto_similarity_list_precomputed(qtof, similarity_threshold=0.0)
    # Save to csv
    # qtof.to_csv('./summary/Structural_Similarity_Prediction.csv', index=False)
    # # Save to parquet
    # parquet_as_df[parquet_as_df.spectrum_id.isin(qtof.spectrum_id)].export_parquet('./spectra/Structural_Similarity_Prediction.parquet')
    
    df = pd.read_csv(summary_path, dtype=   {'Smiles':str,
                                            'msDetector':str,
                                            'msDissociationMethod':str,
                                            'msManufacturer':str,
                                            'msMassAnalyzer':str,
                                            'msModel':str})
    df.Smiles = df.Smiles.astype(str)
    df = df[(df.Smiles.notnull()) & (df.Smiles != 'nan')]
    df = generate_fingerprints(df)

    # Compute and Save Similarities
    build_tanimoto_similarity_list_precomputed(df, './util/Structural_Similarity_Prediction_Pairs.csv', similarity_threshold=0.0)
    
    if type(df) == dd.DataFrame:
        final_ids = df.spectrum_id.values.compute()
    else:
        final_ids = df.spectrum_id.values

    # Save to csv
    df.to_csv('./summary/Structural_Similarity_Prediction.csv', index=False)
    # Save to parquet
    # Make sure to call compute here to load the values into memory
    parquet_as_df[parquet_as_df.spectrum_id.isin(final_ids)].export_parquet('./spectra/Structural_Similarity_Prediction.parquet')

def Spectral_Similarity_Prediction(summary_path:str, parquet_path:str):
    # This dataset: struct1 + struct2 -> spectral similarity
    parquet_as_df = vaex.open(parquet_path)
            
    df = pd.read_csv(summary_path)
    df = df.loc[~df.Smiles.isna()]
    # df = df.loc[df.Adduct == 'M+H']
    # qtof = df.loc[(df.msMassAnalyzer == 'qtof') & (df.GNPS_Inst == 'qtof')]    
    
    # Save to csv
    # qtof.to_csv('./summary/Spectral_Similarity_Prediction.csv', index=False)
    df.to_csv('./summary/Spectral_Similarity_Prediction.csv', index=False)
    # Save to parquet
    parquet_as_df[parquet_as_df.spectrum_id.isin(df.spectrum_id)].export_parquet('./spectra/Spectral_Similarity_Prediction.parquet')
    
    # We will use the networking barebones workflow to generate the similarities (see workflow)
    
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
    # Generate Fingerprints
    positive = generate_fingerprints(positive)
    negative = generate_fingerprints(negative)
    # Compute and Save Similarities to utils
    positive_sim_path = './util/Structural_Modification_orbitrap_positive_sim.csv'
    negative_sim_path = './util/Structural_Modification_orbitrap_negative_sim.csv'
    build_tanimoto_similarity_list_precomputed(positive, positive_sim_path, similarity_threshold = 0.70)
    build_tanimoto_similarity_list_precomputed(negative, negative_sim_path, similarity_threshold = 0.70)
    # Remove entires not included in positive_sim
    positive_sim = vaex.open(positive_sim_path)
    negative_sim = vaex.open(negative_sim_path)
    positive_ids = pd.concat((positive_sim.spectrumid1.to_pandas_series().unique(),positive_sim.spectrumid2.to_pandas_series().unique())).unique()
    negative_ids = pd.concat((negative_sim.spectrumid1.to_pandas_series().unique(),negative_sim.spectrumid2.to_pandas_series().unique())).unique()
    positive = positive[positive.spectrum_id.isin(positive_ids)]
    negative = negative[negative.spectrum_id.isin(negative_ids)]
    # Save to csv
    positive.to_csv('./summary/Structural_Modification_orbitrap_positive.csv', index=False)
    negative.to_csv('./summary/Structural_Modification_orbitrap_negative.csv', index=False)
    # Save to parquet
    parquet_as_df[parquet_as_df.spectrum_id.isin(positive_ids)].export_parquet('./spectra/Structural_Modification_orbitrap_positive.parquet')
    parquet_as_df[parquet_as_df.spectrum_id.isin(negative_ids)].export_parquet('./spectra/Structural_Modification_orbitrap_negative.parquet')
    
    # qtof Portion
    positive = df.loc[(df.Adduct == 'M+H') & (df.msMassAnalyzer == 'qtof') & (df.GNPS_Inst == 'qtof')]
    negative = df.loc[(df.Adduct == 'M-H') & (df.msMassAnalyzer == 'qtof') & (df.GNPS_Inst == 'qtof')]
    # Remove entries with duplicate SMILES
    positive = positive.drop_duplicates(subset=['Smiles'])
    negative = negative.drop_duplicates(subset=['Smiles'])
    # Add fingerprints
    positive = generate_fingerprints(positive)
    negative = generate_fingerprints(negative)
    # Calculate and save structural similarities
    positive_sim_path = './util/Structural_Modification_qtof_positive_sim.csv'
    negative_sim_path = './util/Structural_Modification_qtof_negative_sim.csv'
    build_tanimoto_similarity_list_precomputed(positive, positive_sim_path, similarity_threshold = 0.70)
    build_tanimoto_similarity_list_precomputed(negative, negative_sim_path, similarity_threshold = 0.70)
    # Remove entires not included in positive_sim
    positive_sim = vaex.open(positive_sim_path)
    negative_sim = vaex.open(negative_sim_path)
    positive_ids = pd.concat((positive_sim.spectrumid1.to_pandas_series().unique(),positive_sim.spectrumid2.to_pandas_series().unique())).unique()
    negative_ids = pd.concat((negative_sim.spectrumid1.to_pandas_series().unique(),negative_sim.spectrumid2.to_pandas_series().unique())).unique()
    positive = positive[positive.spectrum_id.isin(positive_ids)]
    negative = negative[negative.spectrum_id.isin(negative_ids)]
    # Save to csv
    positive.to_csv('./summary/Structural_Modification_qtof_positive.csv', index=False)
    negative.to_csv('./summary/Structural_Modification_qtof_negative.csv', index=False)
    # Save to parquet
    parquet_as_df[parquet_as_df.spectrum_id.isin(positive_ids)].export_parquet('./spectra/Structural_Modification_qtof_positive.parquet')
    parquet_as_df[parquet_as_df.spectrum_id.isin(negative_ids)].export_parquet('./spectra/Structural_Modification_qtof_negative.parquet')
    

def main():
    subsets = ['Bruker_Fragmentation_Prediction','MH_MNA_Translation','Orbitrap_Fragmentation_Prediction',\
                'Thermo_Bruker_Translation','Structural_Modification','Structural_Similarity_Prediction',\
                'Spectral_Similarity_Prediction', \
                'GNPS_default']
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
    elif args.subset == 'Fingerprint_Prediction':
        Fingerprint_Prediction(csv_path, parquet_path)
    elif args.subset == 'Orbitrap_Fragmentation_Prediction':
        Orbitrap_Fragmentation_Prediction(csv_path, parquet_path)
    elif args.subset == 'Thermo_Bruker_Translation':
        Thermo_Bruker_Translation(csv_path, parquet_path)
    elif args.subset == 'Structural_Modification':
        Structural_Modification(csv_path, parquet_path)
    elif args.subset == 'Structural_Similarity_Prediction':
        Structural_Similarity_Prediction(csv_path, parquet_path)
    elif args.subset == 'Spectral_Similarity_Prediction':
        Spectral_Similarity_Prediction(csv_path, parquet_path)
    elif args.subset == 'GNPS_default':
        Bruker_Fragmentation_Prediction(csv_path, parquet_path)
        MH_MNA_Translation(csv_path, parquet_path)
        Fingerprint_Prediction(csv_path, parquet_path)
        Orbitrap_Fragmentation_Prediction(csv_path, parquet_path)
        Thermo_Bruker_Translation(csv_path, parquet_path)
        Structural_Modification(csv_path, parquet_path)
        Structural_Similarity_Prediction(csv_path, parquet_path)
        Spectral_Similarity_Prediction(csv_path, parquet_path)
        
            
if __name__ == '__main__':
    main()