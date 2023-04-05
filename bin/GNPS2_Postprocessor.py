import ast
import os
import numpy as np
import pandas as pd
from pyteomics import mgf
from pyteomics.mgf import IndexedMGF
from glob import glob
import re
from pathlib import Path 
import datetime
from rdkit.Chem import AllChem
from rdkit import Chem
from tqdm import tqdm

def sanity_checks(summary):
    assert len(summary[(summary.msManufacturer == 'Thermo') & (summary.msMassAnalyzer == 'qtof')]) == 0
    assert len(summary[(summary.msMassAnalyzer == 'orbitrap') & (summary.msManufacturer == 'Bruker Daltonics')]) == 0 

def basic_cleaning(summary):
    # scan
    summary.scan = summary.scan.astype('int')

    # smiles
    summary.Smiles = summary.Smiles.astype(str).apply(lambda x: x.strip() )
    summary.Smiles = summary.Smiles.apply(lambda x: 'nan' if x == '' or 'N/A' in x else x)
    
    # ionization
    summary.msIonisation = summary.msIonisation.astype(str)
    summary.loc[['esi' in x.lower() or 'electrospray' in x.lower() for x in summary.msIonisation], 'msIonisation'] = "ESI"

    # retention time
    assert pd.Series([x[-1] == 'S' for x in summary.retention_time[~ summary.retention_time.isna()]]).all() # Sanity check because right now everything is in seconds
    pattern = re.compile("([0-9]*[.])?[0-9]+")
    summary.loc[~ summary.retention_time.isna(),"retention_time"] = summary.loc[~ summary.retention_time.isna(),"retention_time"].apply(lambda x: re.search(pattern,x).group(0))

    # Adduct
    summary.Adduct = summary.Adduct.apply(lambda x: str(x).strip())
    # Strip starting and ending braces if no charge is specified
    summary.Adduct = summary.Adduct.map(lambda x: x[1:-1] if x.strip()[-1] == ']' and x[0] == '[' else x)
    # Strip starting and ending braces if no charge is specified
    summary.Adduct = summary.Adduct.map(lambda x: x[1:-2] if (x.strip()[-2:] == ']+' or x[-2:] == ']-') and x[0] == '[' else x)

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
    
    summary.msMassAnalyzer = summary.msMassAnalyzer.str.lower()
    # summary.msMassAnalyzer = summary.msMassAnalyzer.str.strip('[]').str.strip("'").str.split(',')
    mask = ~ summary.msMassAnalyzer.isna()
    summary.loc[mask,'msMassAnalyzer'] = summary.loc[mask,'msMassAnalyzer'].apply(lambda x: ast.literal_eval(x))
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

def propogate_GNPS_Inst_field(summary):
    summary.GNPS_Inst = summary.GNPS_Inst.astype('str')
    summary.GNPS_Inst = summary.GNPS_Inst.map(lambda x: x.strip().lower())

    # Fragmentation Info (Done)
    summary.msDissociationMethod = summary.msDissociationMethod.astype(str)
    summary.loc[pd.Series([("in source cid" == x) for x in summary.GNPS_Inst]) & summary.msDissociationMethod == 'nan', 'msDissociationMethod'] = "is-cid"    
    summary.loc[pd.Series([("hid" in x) for x in summary.GNPS_Inst]) & summary.msDissociationMethod == 'nan', 'msDissociationMethod'] = "hid"    
    summary.loc[pd.Series([("cid" in x and not "is-cid" in x) for x in summary.GNPS_Inst]) & summary.msDissociationMethod == 'nan', 'msDissociationMethod'] = "cid"    

    # Ionisation Info (Not Done)
    summary.loc[pd.Series(["esi" in x for x in  summary.GNPS_Inst]) & summary.msIonisation == 'nan', 'msIonisation'] = 'ESI'
    summary.loc[pd.Series(["apci" in x for x in  summary.GNPS_Inst]) & summary.msIonisation == 'nan', 'msIonisation'] = 'APCI'
    summary.loc[pd.Series([("appi" in x and not "dappi" in x) for x in summary.GNPS_Inst]) & summary.msIonisation == 'nan', 'msIonisation'] = 'APPI'
    summary.loc[pd.Series(["dappi" in x for x in summary.GNPS_Inst]) & summary.msIonisation == 'nan', 'msIonisation'] = 'DAPPI'

    # Mass Analyzer (Not Done)
    summary.loc[pd.Series(["orbitrap" in x for x in summary.GNPS_Inst]) & summary.msMassAnalyzer == 'nan',"msMassAnalyzer"] = "orbitrap"
    summary.loc[pd.Series([("quadrupole tof" in x or "qtof" in x or "q-tof" in x) and not "qq" in x for x in summary.GNPS_Inst]) & summary.msMassAnalyzer == 'nan',"msMassAnalyzer"] = "qtof"
    summary.loc[pd.Series([("tof" in x) and not ("qq" in x or "qtof" in x or "q-tof" in x or "q tof" in x or "quadrupole tof" in x) for x in summary.GNPS_Inst]) & summary.msMassAnalyzer == 'nan',"msMassAnalyzer"] = "tof"

    # Manufacturer Info (Not Done)
    summary.loc[pd.Series(["maxis" in x for x in summary.GNPS_Inst]) & summary.msManufacturer == "nan","msManufacturer"] = "Bruker Daltonics"
    summary.loc[pd.Series(["q exactive" in x or "q-exactive" in x for x in summary.GNPS_Inst]) & summary.msManufacturer == "nan","msManufacturer"] = "Thermo"
    return summary

def propogate_msModel_field(summary):
    summary.msModel = summary.msModel.astype(str)

    summary.loc[["maxis" in x.lower() or "bruker daltonics" in x.lower() for x in summary.msModel],"msManufacturer"] = "Bruker Daltonics"
    summary.loc[["agilent" in x.lower() for x in summary.msModel],"msManufacturer"] = "Agilent"
    summary.loc[["waters" in x.lower() for x in summary.msModel],"msManufacturer"] = "Waters"
    summary.loc[["shimadzu" in x.lower() for x in summary.msModel],"msManufacturer"] = "Shimadzu"
    summary.loc[["q exactive" in x.lower() or "q-exactive" in x.lower() for x in summary.msModel],"msManufacturer"] = "Thermo"

    return summary

def generate_fingerprints(summary):
    summary = summary.assign(mol=None)
    summary.loc[summary.Smiles != 'nan', 'mol'] = summary.loc[summary.Smiles != 'nan', 'Smiles'].apply(lambda x: Chem.MolFromSmiles(x, sanitize=True))
    smiles_parsing_success = [x is not None for x in summary.mol]
    summary.loc[smiles_parsing_success,'Morgan_2048_2'] = summary.loc[smiles_parsing_success,'mol'].apply( \
        lambda x: list(AllChem.GetMorganFingerprintAsBitVect(x,2,useChirality=False,nBits=2048)))
    summary.loc[smiles_parsing_success,'Morgan_4096_2'] = summary.loc[smiles_parsing_success,'mol'].apply( \
        lambda x: list(AllChem.GetMorganFingerprintAsBitVect(x,2,useChirality=False,nBits=4096)))
    summary.loc[smiles_parsing_success,'Morgan_2048_3'] = summary.loc[smiles_parsing_success,'mol'].apply( \
        lambda x: list(AllChem.GetMorganFingerprintAsBitVect(x,3,useChirality=False,nBits=2048)))
    summary.loc[smiles_parsing_success,'Morgan_4096_3'] = summary.loc[smiles_parsing_success,'mol'].apply( \
        lambda x: list(AllChem.GetMorganFingerprintAsBitVect(x,3,useChirality=False,nBits=4096)))
    # summary.loc[smiles_parsing_success,'RdKit_2048_5'] = summary.loc[smiles_parsing_success,'mol'].apply( \
    #     lambda x: Chem.RDKFingerprint(x,minPath=5,fpSize=2048))
    # summary.loc[smiles_parsing_success,'RdKit_4096_5'] = summary.loc[smiles_parsing_success,'mol'].apply( \
    #     lambda x: Chem.RDKFingerprint(x,minPath=5,fpSize=4096))
    summary.drop('mol', axis=1,inplace=True)
    return summary

def generate_parquet_file(input_mgf, spectrum_ids):
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
    
    # Exploiting GNPS_Inst annotations:
    summary = propogate_GNPS_Inst_field(summary)

    # Exploiting Some of the info in msModel
    summary = propogate_msModel_field(summary)
    sanity_checks(summary)
    
    # Add Fingerprints
    summary = generate_fingerprints(summary)

    parquet_as_df = generate_parquet_file(mgf_path, summary.spectrum_id.astype('str'))
    parquet_as_df.to_parquet(output_parquet_path)
    summary.to_csv(output_csv_path, index=False)

def main():
    csv_path = "ALL_GNPS_merged.csv"
    mgf_path = "ALL_GNPS_merged.mgf"
    cleaned_csv_path = "ALL_GNPS_cleaned.csv"
    cleaned_parquet_path = "ALL_GNPS_cleaned.parquet"

    if not os.path.isfile(cleaned_csv_path):
        if not os.path.isfile(cleaned_parquet_path):
            postprocess_files(csv_path, mgf_path, cleaned_csv_path, cleaned_parquet_path)
            
if __name__ == '__main__':
    main()