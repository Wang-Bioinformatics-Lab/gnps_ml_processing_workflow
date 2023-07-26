import ast
import math
import os
import pickle
import numpy as np
import pandas as pd
from pyteomics.mgf import IndexedMGF
import re
from tqdm import tqdm
from utils import harmonize_smiles_rdkit
from rdkit import Chem

import sys

# Modify sys.path to include the parent directory
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 'formula_validation'))
sys.path.append(parent_dir)
print(sys.path)

# Import module2 using relative import
from Formula import Formula
from Adduct import Adduct
from IncorrectFormula import IncorrectFormula
from IncorrectAdduct import IncorrectAdduct

# Now you can use module2 in your code

# Restore sys.path to its original state if needed
sys.path.remove(parent_dir)



def basic_cleaning(summary):
    # scan
    summary.scan = summary.scan.astype('int')

    # smiles
    summary.Smiles = summary.Smiles.astype(str).apply(lambda x: x.strip() )
    summary.Smiles = summary.Smiles.apply(lambda x: 'nan' if x == '' or 'N/A' in x else x)
    
    # ionization
    summary.msIonisation = summary.msIonisation.astype(str)
    summary.loc[['esi' in x.lower() or 'electrospray' in x.lower() for x in summary.msIonisation], 'msIonisation'] = "ESI"
    summary.loc[['' == x or 'positive' == x.lower() or 'negative'  == x.lower() for x in summary.msIonisation], 'msIonisation'] = "nan"
    
    # Adduct translation from chemical names to adduct formulas -> 'M+TFA-H': '[M+C2HF3O2-H]-'
    # Adduct Table Credit: Yasin El Abiead
    adduct_mapping = pickle.load(open('./adduct_mapping.pkl', 'rb'))
    summary.Adduct = summary.Adduct.apply(lambda x: adduct_mapping.get(x.strip()))
    summary = summary[summary.Adduct.notna()]    

    # Conversion of numerical columns to numerical types to protect against contamination
    summary.Precursor_MZ = summary.Precursor_MZ.astype(float)
    summary.ExactMass = summary.ExactMass.astype(float)
    summary.Charge = summary.Charge.astype(int)

    # Ion Mode
    summary.Ion_Mode = summary.Ion_Mode.apply(lambda x: str(x).strip().lower())
    summary.loc[summary.Ion_Mode == 'positive-20ev','Ion_Mode'] = 'positive'

    # Manufacturer
    summary.msManufacturer = summary.msManufacturer.astype('str')
    summary.loc[['thermo' in x.lower() for x in summary.msManufacturer], 'msManufacturer'] = 'Thermo'

    # msMassAnalyzer
    def transform_analyzer(x:str):
        if 'quadrupole tof' == x:
            return 'qtof'
        if 'fourier transform ion cyclotron resonance mass spectrometer' == x:
            return 'ftms'
        if 'time-of-flight' == x:
            return 'tof'
        return x
    def merge_analyzer(x:list):
        if np.any(["orbitrap" in y for y in x]):
            return "orbitrap"
        if "quadrupole" in x and "tof" in x:
            return 'qtof'
        return x[0]
    
    summary.msMassAnalyzer = summary.msMassAnalyzer.astype('str')
    summary.msMassAnalyzer = summary.msMassAnalyzer.str.lower()
    
    # summary.msMassAnalyzer = summary.msMassAnalyzer.str.strip('[]').str.strip("'").str.split(',')
    mask = (summary.msMassAnalyzer.notna()) & (summary.msMassAnalyzer != 'nan')
    def literal_eval_helper(x):
        """
        An small helper function to handle weird entries in msMassAnalyzer
        """
        try:
            return ast.literal_eval(x)
        except Exception as e:
            print(e)
            print(f"Error in literal_eval_helper when trying to parse {x}")
            return ["nan"]
    
    summary.loc[mask,'msMassAnalyzer'] = summary.loc[mask,'msMassAnalyzer'].apply(lambda x: literal_eval_helper(x))
    summary.loc[mask,'msMassAnalyzer'] = summary.loc[mask,'msMassAnalyzer'].apply(lambda x: [transform_analyzer(y) for y in x])
    summary.loc[mask,'msMassAnalyzer'] = summary.loc[mask,'msMassAnalyzer'].apply(merge_analyzer)
    # summary.loc[mask].apply(lambda x: [y == 'quadrupole tof' for y in x].any())
    # summary.loc[mask,'msMassAnalyzer'] = 'qtof' 
    # summary.loc[summary.msMassAnalyzer == 'fourier transform ion cyclotron resonance mass spectrometer','msMassAnalyzer'] = 'ftms'

    # Detector
    '''A very specific set of files has MS:1000253 as the detector name which is used for mzML files, 
    however, it's written in mzXML files. We'll clean this up here.

    Conversion source: https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo
    '''

    summary.loc[['MS:1000253' == x for x in summary.msDetector], 'msDetector'] = 'EMT'
    summary.loc[summary.msDetector == "electron multiplier", 'msDetector'] = 'EMT'
    summary.loc[summary.msDetector == "microchannel plate detector", 'msDetector'] = 'MCP'
    
    
    return summary

def clean_smiles(summary):
    """This function will harmonize the tautomers for the smiles strings and remove invalid strings.

    Args:
        summary (DataFrame): The summary dataframe

    Returns:
        DataFrame: The modified summary dataframe
    """
    summary.Smiles = summary.Smiles.astype(str)
    summary.Smiles = summary.Smiles.apply(lambda x: harmonize_smiles_rdkit(x))
    return summary

def propogate_GNPS_Inst_field(summary):
    summary.GNPS_Inst = summary.GNPS_Inst.astype('str')
    summary.GNPS_Inst = summary.GNPS_Inst.map(lambda x: x.strip().lower())

    # Fragmentation Info (Done)
    summary.msDissociationMethod = summary.msDissociationMethod.astype(str)
    summary.loc[(pd.Series(["in source cid" == x for x in summary.GNPS_Inst]) & (summary.msDissociationMethod == 'nan')), 'msDissociationMethod'] = "is-cid"
   
    summary.loc[(pd.Series([("hid" in x) for x in summary.GNPS_Inst]) & (summary.msDissociationMethod == 'nan')), 'msDissociationMethod'] = "hid"    
    summary.loc[(pd.Series([("cid" in x and not "is-cid" in x) for x in summary.GNPS_Inst]) & (summary.msDissociationMethod == 'nan')), 'msDissociationMethod'] = "cid"    

    # Ionisation Info (Not Done)
    summary.loc[(pd.Series(["esi" in x for x in  summary.GNPS_Inst]) & (summary.msIonisation == 'nan')), 'msIonisation'] = 'ESI'
    summary.loc[(pd.Series(["apci" in x for x in  summary.GNPS_Inst]) & (summary.msIonisation == 'nan')), 'msIonisation'] = 'APCI'
    summary.loc[(pd.Series([("appi" in x and not "dappi" in x) for x in summary.GNPS_Inst]) & (summary.msIonisation == 'nan')), 'msIonisation'] = 'APPI'
    summary.loc[(pd.Series(["dappi" in x for x in summary.GNPS_Inst]) & (summary.msIonisation == 'nan')), 'msIonisation'] = 'DAPPI'

    # Mass Analyzer (Not Done)
    summary.loc[(pd.Series(["orbitrap" in x for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == 'nan')),"msMassAnalyzer"] = "orbitrap"
    summary.loc[(pd.Series([("quadrupole tof" in x or "qtof" in x or "q-tof" in x) and not "qq" in x for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == 'nan')),"msMassAnalyzer"] = "qtof"
    summary.loc[(pd.Series([("tof" in x) and not ("qq" in x or "qtof" in x or "q-tof" in x or "q tof" in x or "quadrupole tof" in x) for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == 'nan')),"msMassAnalyzer"] = "tof"

    # Manufacturer Info (Not Done)
    summary.loc[(pd.Series(["maxis" in x for x in summary.GNPS_Inst]) & (summary.msManufacturer == "nan")),"msManufacturer"] = "Bruker Daltonics"
    summary.loc[(pd.Series(["q exactive" in x or "q-exactive" in x for x in summary.GNPS_Inst]) & (summary.msManufacturer == "nan")),"msManufacturer"] = "Thermo"
    return summary

def propogate_msModel_field(summary):
    summary.msModel = summary.msModel.astype(str)

    summary.loc[["maxis" in x.lower() or "bruker daltonics" in x.lower() for x in summary.msModel],"msManufacturer"] = "Bruker Daltonics"
    summary.loc[["agilent" in x.lower() for x in summary.msModel],"msManufacturer"] = "Agilent"
    summary.loc[["waters" in x.lower() for x in summary.msModel],"msManufacturer"] = "Waters"
    summary.loc[["shimadzu" in x.lower() for x in summary.msModel],"msManufacturer"] = "Shimadzu"
    summary.loc[["q exactive" in x.lower() or "q-exactive" in x.lower() for x in summary.msModel],"msManufacturer"] = "Thermo"

    return summary


def sanity_checks(summary):
    assert len(summary[(summary.msManufacturer == 'Thermo') & (summary.msMassAnalyzer == 'qtof')]) == 0
    assert len(summary[(summary.msMassAnalyzer == 'orbitrap') & (summary.msManufacturer == 'Bruker Daltonics')]) == 0 
    assert len(summary[(summary.Adduct == 'None') & (summary.Adduct == 'nan') & (summary.Adduct.isna())]) == 0


def add_columns_formula_analysis(summary): 
    column_name_ppmBetweenExpAndThMass='ppmBetweenExpAndThMass'
    
    def helper(row):
        try:
            return Formula.formula_from_smiles(row['Smiles'], row['Adduct']).ppm_difference_with_exp_mass(row['Precursor_MZ'])
        except IncorrectFormula as incFor:
            return 'nan'
            
    summary[column_name_ppmBetweenExpAndThMass] = summary.apply(helper, axis=1)

def generate_parquet_df(input_mgf, spectrum_ids):
    """
    Details on output format:
    Columns will be [level_0, index, i, i_norm, mz, precmz]
    Index will be spectrum_id
    level_0 is the row index in file
    index is the row index in the spectra
    """
    output = []
    indexed_mgf = IndexedMGF(input_mgf,index_by_scans=True)
    level_0 = 0
    for m in tqdm(indexed_mgf):
        spectrum_id = m['params']['title']
        mz_array = m['m/z array']
        intensity_array = m['intensity array']
        precursor_mz = m['params']['pepmass']
        if spectrum_id in spectrum_ids.values: # Make sure that it didn't get removed during cleaning
            # charge = m['charge']
            for index, (mz, intensity) in enumerate(zip(mz_array, intensity_array)):
                output.append({'spectrum_id':spectrum_id, 'level_0': level_0, 'index':index, 'i':intensity, 'mz':mz, 'prec_mz':precursor_mz})
                level_0 += 1
                
    output = pd.DataFrame(output)
    output.set_index('spectrum_id')
    return output
            
                       

def postprocess_files(csv_path, mgf_path, output_csv_path, output_parquet_path):
    summary = pd.read_csv(csv_path)

    # Cleaning up files:
    summary = basic_cleaning(summary)
    
    # Clean smiles strings:
    summary = clean_smiles(summary)
    
    # Exploiting GNPS_Inst annotations:
    summary = propogate_GNPS_Inst_field(summary)

    # Exploiting Some of the info in msModel
    summary = propogate_msModel_field(summary)
    sanity_checks(summary)

    add_columns_formula_analysis(summary)
    
    #parquet_as_df = generate_parquet_df(mgf_path, summary.spectrum_id.astype('str'))
    #parquet_as_df.to_parquet(output_parquet_path)
    summary.to_csv(output_csv_path, index=False)

def main():
    csv_path = "/home/alberto/Downloads/sample_csv.csv"
    mgf_path = "/home/alberto/Downloads/sample_mgf.mgf"
    cleaned_csv_path = "ALL_GNPS_cleaned.csv"
    cleaned_parquet_path = "ALL_GNPS_cleaned.parquet"

    if not os.path.isfile(cleaned_csv_path):
        if not os.path.isfile(cleaned_parquet_path):
            postprocess_files(csv_path, mgf_path, cleaned_csv_path, cleaned_parquet_path)
            
if __name__ == '__main__':
    main()