import ast
import math
import os
import pickle
import numpy as np
import pandas as pd
from pyteomics.mgf import IndexedMGF
import re
from tqdm import tqdm
from utils import harmonize_smiles_rdkit, INCHI_to_SMILES, synchronize_spectra, generate_parquet_df
from rdkit import Chem
from pandarallel import pandarallel
import time
import datetime
import argparse

PARALLEL_WORKERS = 32

import sys

from formula_validation.Formula import Formula
from formula_validation.Adduct import Adduct
from formula_validation.IncorrectFormula import IncorrectFormula
from formula_validation.IncorrectAdduct import IncorrectAdduct

# os.environ['JOBLIB_TEMP_FOLDER'] = '/tmp'

def basic_cleaning(summary):
    # scan
    summary.scan = summary.scan.astype('int')

    # smiles
    summary.Smiles = summary.Smiles.astype(str).parallel_apply(lambda x: x.strip() )
    summary.Smiles = summary.Smiles.parallel_apply(lambda x: 'nan' if (x == '') or ('N/A' in x) else x)
    
    # INCHI
    summary.INCHI = summary.INCHI.astype(str).parallel_apply(lambda x: x.strip() )
    summary.INCHI = summary.INCHI.parallel_apply(lambda x: 'nan' if (x == '') or ('N/A' in str(x)) else x)
    
    # In rare cases the user will use INCHI and not smiles, so we'll convert it to smiles
    mask = (summary.Smiles == 'nan') & (summary.INCHI != 'nan')
    summary.loc[mask, 'Smiles'] = summary.loc[mask, 'Smiles'].parallel_apply(INCHI_to_SMILES)
    
    # ionization
    summary.msIonisation = summary.msIonisation.astype(str)
    summary.loc[['esi' in x.lower() or 'electrospray' in x.lower() for x in summary.msIonisation], 'msIonisation'] = "ESI"
    summary.loc[['' == x or 'positive' == x.lower() or 'negative'  == x.lower() for x in summary.msIonisation], 'msIonisation'] = "nan"
    
    # Adduct translation from chemical names to adduct formulas -> 'M+TFA-H': '[M+C2HF3O2-H]-'
    # Adduct Table Credit: Yasin El Abiead
    adduct_mapping = pickle.load(open('./adduct_mapping.pkl', 'rb'))
    def helper(adduct):
        adduct = str(adduct).strip()
        mapped_adduct = None
        # Map the adduct if possible, if not, then we leave the orignal value
        if adduct not in adduct_mapping.values():
            # If the adduct is not in the values, attempt to map it:
            mapped_adduct = adduct_mapping.get(adduct)
            if mapped_adduct is not None:
                adduct = mapped_adduct
            
        # Add 1 if not None and last is
        if adduct is not None:
            if adduct[-2:] == "]+":
                adduct = adduct[:-2] + "]1+"
            elif adduct[-2:] == "]-":
                adduct = adduct[:-2] + "]1-"
            
        return adduct
    summary.Adduct = summary.Adduct.parallel_apply(helper)
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
    summary.Smiles = summary.Smiles.parallel_apply(lambda x: harmonize_smiles_rdkit(x, skip_tautomerization=True))
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
    # assert len(summary[(summary.Adduct == 'None') & (summary.Adduct == 'nan') & (summary.Adduct.isna())]) == 0 # Right now because of the adduct cleaning, there is a chance that we'll have nan adducts


def add_columns_formula_analysis(summary): 
    column_name_ppmBetweenExpAndThMass='ppmBetweenExpAndThMass'
    
    def helper(row):
        try:
            smiles = str(row['Smiles'])
            if row['Smiles'] != 'nan':
                formula = Formula.formula_from_smiles(smiles, row['Adduct'], no_api=False)  # Disabling API can improve speed 
                if formula is not None:
                    return float(formula.ppm_difference_with_exp_mass(row['Precursor_MZ']))
                else:
                    return np.nan
        except IncorrectFormula as incFor:
            return np.nan
        except IncorrectAdduct as incAdd:
            return np.nan
        except Exception as e:
            print(e, file=sys.stderr)
            return np.nan

    summary[column_name_ppmBetweenExpAndThMass] = summary.parallel_apply(helper, axis=1).astype(float)
    
def add_explained_intensity(summary, spectra):
    indexed_mgf = IndexedMGF(spectra,index_by_scans=True)
    
    column_name_ppmBetweenExpAndThMass='explainable_intensity'
    
    def helper(row):
        try:
            smiles = str(row['Smiles'])
            if smiles != 'nan':
                # Build a dictionary of mz, intensity
                this_spectra = indexed_mgf[row['scan']]
                if this_spectra['params']['title'] != row['spectrum_id']:
                    raise ValueError(f"Spectrum ID mismatch: {this_spectra['params']['title']} and {row['spectrum_id']}")
               
                mzs = this_spectra['m/z array']
                intensities = this_spectra['intensity array']
                fragments_mz_intensities = dict(zip(mzs, intensities))
                return Formula.formula_from_smiles(smiles, row['Adduct'], metadata={'ccms_id':row['spectrum_id']}).percentage_intensity_fragments_explained_by_formula(fragments_mz_intensities, ppm=50)
        except IncorrectFormula as incFor:
            return 'nan'
        except IncorrectAdduct as incAdd:
            return 'nan'
        except Exception as e:
            print(e, file=sys.stderr)
            return 'nan'
            
    filter = (summary['ppmBetweenExpAndThMass'].notna() & summary['ppmBetweenExpAndThMass']<=50)
    summary.loc[filter, column_name_ppmBetweenExpAndThMass] = summary.loc[filter].parallel_apply(helper, axis=1)
            
def postprocess_files(csv_path, mgf_path, output_csv_path, output_parquet_path, cleaned_mgf_path):
    pandarallel.initialize(progress_bar=False, nb_workers=PARALLEL_WORKERS, use_memory_fs = False)
    
    summary = pd.read_csv(csv_path)

    # Cleaning up files:
    print("Performing basic cleaning", flush=True)
    start = time.time()
    summary = basic_cleaning(summary)
    print("Done in {} seconds".format(datetime.timedelta(seconds=time.time() - start)), flush=True)
    
    # Clean smiles strings:
    print("Cleaning up smiles strings", flush=True)
    start = time.time()
    summary = clean_smiles(summary)
    print("Done in {} seconds".format(datetime.timedelta(seconds=time.time() - start)), flush=True)
    
    # Exploiting GNPS_Inst annotations:
    print("Attempting to propogate user instrument annotations", flush=True)
    start = time.time()
    summary = propogate_GNPS_Inst_field(summary)
    print("Done in {} seconds".format(datetime.timedelta(seconds=time.time() - start)), flush=True)

    # Exploiting Some of the info in msModel
    print("Attempting to propogate msModel field", flush=True)
    start = time.time()
    summary = propogate_msModel_field(summary)
    print("Done in {} seconds".format(datetime.timedelta(seconds=time.time() - start)), flush=True)

    # Calculating ppm error
    print("Calculating ppm error", flush=True)
    start = time.time()
    add_columns_formula_analysis(summary)
    print("Done in {} seconds".format(datetime.timedelta(seconds=time.time() - start)), flush=True)
    
    # Cleanup scan numbers, scan numbers will be reset in mgf by synchronize spectra
    summary.scan = np.arange(0, len(summary))
        
    # Cleanup MGF file. Must be done before explained intensity calculations in order to make sure spectra are in order
    print("Writing mgf file", flush=True)
    start = time.time()
    synchronize_spectra(mgf_path, cleaned_mgf_path, summary.spectrum_id.astype('str').values)
    print("Done in {} seconds".format(datetime.timedelta(seconds=time.time() - start)), flush=True)

    # # Because calculating explained intensity is slow, we'll save an output file just in case
    # print("Writing csv file", flush=True)
    # summary.to_csv(output_csv_path, index=False)

    # # Calculate explained intensity
    # print("Calculating explained intensity")
    # start = time.time()
    # add_explained_intensity(summary, mgf_path)
    # print("Done in {} seconds".format(datetime.timedelta(seconds=time.time() - start)))

    sanity_checks(summary)
    
    print("Writing output files...", flush=True)
    print("Writing parquet file", flush=True)
    parquet_as_df = generate_parquet_df(cleaned_mgf_path)
    parquet_as_df.to_parquet(output_parquet_path)
    print("Writing csv file", flush=True)
    summary.to_csv(output_csv_path, index=False)
    print("Postprocessing complete!", flush=True)

def main():
    parser = argparse.ArgumentParser(description='Postprocess GNPS files')
    parser.add_argument('--input_csv_path', type=str, default="ALL_GNPS_merged.csv", help='Path to the csv file')
    parser.add_argument('--input_mgf_path', type=str, default="ALL_GNPS_merged.mgf", help='Path to the mgf file')
    parser.add_argument('--output_csv_path', type=str, default="ALL_GNPS_cleaned.csv", help='Path to the output csv file')
    parser.add_argument('--output_parquet_path', type=str, default="ALL_GNPS_cleaned.parquet", help='Path to the output parquet file')
    parser.add_argument('--output_mgf_path', type=str, default="ALL_GNPS_cleaned.mgf", help='Path to the output mgf file')
    args= parser.parse_args()
    
    csv_path                = str(args.input_csv_path)
    mgf_path                = str(args.input_mgf_path)
    cleaned_csv_path        = str(args.output_csv_path)
    cleaned_parquet_path    = str(args.output_parquet_path)
    cleaned_mgf_path        = str(args.output_mgf_path)

    if not os.path.isfile(cleaned_csv_path):
        if not os.path.isfile(cleaned_parquet_path):
            postprocess_files(csv_path, mgf_path, cleaned_csv_path, cleaned_parquet_path, cleaned_mgf_path)
            
if __name__ == '__main__':
    main()