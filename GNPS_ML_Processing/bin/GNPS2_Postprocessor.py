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
    summary.Smiles = summary.Smiles.parallel_apply(lambda x: '' if ('N/A' in x) or ('nan' in x) else x)
    
    # check smiles validity
    summary.Smiles = summary.Smiles.apply(lambda x: x if Chem.MolFromSmiles(x) is not None else '')
    
    # INCHI
    summary.INCHI = summary.INCHI.astype(str).parallel_apply(lambda x: x.strip() )
    summary.INCHI = summary.INCHI.parallel_apply(lambda x: '' if ('N/A' in x) or ('nan' in x) else x)

    # ionization
    summary.msIonisation = summary.msIonisation.astype(str)
    summary.loc[['esi' in x.lower() or 'electrospray' in x.lower() for x in summary.msIonisation], 'msIonisation'] = "ESI"
    summary.loc[['' == x or 'positive' == x.lower() or 'negative'  == x.lower() for x in summary.msIonisation], 'msIonisation'] = "nan"
    
    # Adduct translation from chemical names to adduct formulas -> 'M+TFA-H': '[M+C2HF3O2-H]-'
    # Adduct Table Credit: Yasin El Abiead
    with open("./adduct_mapping.txt", "r") as data:
        adduct_mapping = ast.literal_eval(data.read())
    def helper(adduct):
        adduct = str(adduct).strip()
        adduct = adduct.replace(' ', '')
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
    
    # Get a mask of non-matching adducts
    pattern = r'\[(\d)*M([\+-].*?)\](\d)([\+-])'
    mask = summary.Adduct.apply(lambda x: re.match(pattern, x) is None)
    if sum(mask) > 0:
        print(f"Warning: {sum(mask)} entries have Adducts that are not in the expected format, these will be removed.")
        print(summary.Adduct.loc[mask].value_counts().head(10))
    summary = summary.loc[~mask]

    # Conversion of numerical columns to numerical types to protect against contamination
    summary.Precursor_MZ = summary.Precursor_MZ.astype(float)
    summary.ExactMass = summary.ExactMass.astype(float)
    summary.Charge = summary.Charge.astype(int)
    
    # Charge
    # Mask whether charge is equal to adduct charge
    adduct_charges = summary.Adduct.apply(lambda x: int(x[-1] + x.split(']')[-1][:-1]))
    mask = (summary.Charge != adduct_charges)
    if sum(mask) > 0:
        print(f"Warning: {sum(mask)} entries have Charge and Adduct Charge that are not equivalent, Adduct Charge will be prefered.")
        print(f"Of the {sum(mask)} entires, {sum(mask & (summary.Charge == 0))} have Charge of 0.")
    summary.loc[mask, 'Charge'] = adduct_charges[mask]

    # Collision Energy
    # Rather nicely, sometimes the collision energy is in the ion mode field, but we'll prefer the raw file data
    summary.Ion_Mode = summary.Ion_Mode.apply(lambda x: str(x).strip().lower())   
    mask = (summary.Ion_Mode == 'positive-20ev') & (summary.collision_energy.isna())
    if sum(mask) > 0:
        print(f"Imputing {sum(mask)} collision energies using the Ion_Mode field")
        summary.loc[mask, 'collision_energy'] = 20
    
    # Sometimes the collision energy is in the GNPS_inst field
    pattern = re.compile(r'(\d+)eV')
    extracted_eV = summary.GNPS_Inst.apply(lambda x: re.search(pattern, str(x)))
    extracted_eV = extracted_eV.apply(lambda x: x.group(1) if x is not None else None)
    mask = (extracted_eV.notna()) & (summary.collision_energy.isna())
    if sum(mask) > 0:
        print(f"Imputing {sum(mask)} collision energies using the GNPS_Inst field.")
        summary.loc[mask, 'collision_energy'] = extracted_eV[mask]

    # Ion Mode
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
    
    
    # Check the INCHI and SMILES are equivalent
    mask = (summary.Smiles != '') & (summary.INCHI != '')
    inchi_from_smiles = summary.loc[mask].apply(lambda x: Chem.inchi.MolToInchi(Chem.MolFromSmiles(x['Smiles'])), axis=1)
    equivalency_mask = (inchi_from_smiles != summary.loc[mask, 'INCHI'])
    if sum(equivalency_mask) > 0 :
        print(f"Warning: {sum(equivalency_mask)} entries have INCHI and SMILES that are not equivalent, SMILES values will be used to replace INCHI")
        summary.loc[mask, 'INCHI'].loc[equivalency_mask] = inchi_from_smiles.loc[equivalency_mask]
    
    # In rare cases the user will use INCHI and not smiles, so we'll convert it to smiles
    mask = (summary.Smiles == '') & (summary.INCHI != '')
    summary.loc[mask, 'Smiles'] = summary.loc[mask, 'Smiles'].apply(INCHI_to_SMILES)
    summary.Smiles = summary.Smiles.parallel_apply(lambda x: harmonize_smiles_rdkit(x, skip_tautomerization=True))
    # If no INCHI but we have SMILES, convert it to INCHI
    mask = (summary.Smiles != '') & (summary.INCHI == '')
    summary.loc[mask, 'INCHI'] = summary.loc[mask, 'Smiles'].apply(lambda x: Chem.inchi.MolToInchi(Chem.MolFromSmiles(x)))
    # Fill in INCHI key
    summary.InChIKey_smiles = summary.InChIKey_smiles.astype(str).parallel_apply(lambda x: x.strip() )
    summary.InChIKey_smiles = summary.InChIKey_smiles.parallel_apply(lambda x: '' if ('N/A' in str(x)) else x)
    # Generate all inchi keys
    all_keys = summary.loc[:, 'INCHI'].apply(lambda x: Chem.inchi.InchiToInchiKey(x) if x != '' else '')
    # Check existing keys
    incorrect_mask = (summary.InChIKey_smiles != all_keys)
    if sum(incorrect_mask) > 0 :
        print(f"Warning: {sum(incorrect_mask)} entries have INCHI and INCHIKey that are not equivalent, new INCHIKeys will be generated based on INCHI")
        summary.loc[incorrect_mask, 'InChIKey_smiles'] = summary.loc[incorrect_mask, 'INCHI'].apply(lambda x: Chem.inchi.InchiToInchiKey(x))
    
    return summary

def propogate_GNPS_Inst_field(summary):
    summary.GNPS_Inst = summary.GNPS_Inst.astype('str')
    summary.GNPS_Inst = summary.GNPS_Inst.map(lambda x: x.strip().lower())
    """
    Whenever we create our own series using a list comprehension and use '&' to combine it with something like 
    (summary.msDissociationMethod == 'nan'), we have to have a continuous, zero-indexed dataframe because the '&' 
    joins on the index 
    E.g.:
    >>> a = pd.Series([False,True,True], index=[1,2,3])
    >>> b = pd.Series([True,True,False], index=[2,3,4])
    >>> a & b
    1    False
    2     True
    3     True
    4    False
    dtype: bool
    # Note that there are four entires, not three
    
    A safer solution is to use numpy arrays, which will ignore the index and that's what we'll do here
    """
    summary.reset_index(inplace=True)

    # Fragmentation Info (Done)
    summary.msDissociationMethod = summary.msDissociationMethod.astype(str)
    summary.loc[(np.array(["in source cid" == x for x in summary.GNPS_Inst]) & (summary.msDissociationMethod == 'nan')), 'msDissociationMethod'] = "is-cid"
   
    summary.loc[(np.array([("hid" in x) for x in summary.GNPS_Inst]) & (summary.msDissociationMethod == 'nan')), 'msDissociationMethod'] = "hid"    
    summary.loc[(np.array([("hcd" in x) for x in summary.GNPS_Inst]) & (summary.msDissociationMethod == 'nan')), 'msDissociationMethod'] = "hcd"    
    summary.loc[(np.array([("cid" in x and not "is-cid" in x) for x in summary.GNPS_Inst]) & (summary.msDissociationMethod == 'nan')), 'msDissociationMethod'] = "cid"    

    # Ionisation Info (Not Done)
    summary.loc[(np.array(["esi" in x for x in  summary.GNPS_Inst]) & (summary.msIonisation == 'nan')), 'msIonisation'] = 'ESI'
    summary.loc[(np.array(["apci" in x for x in  summary.GNPS_Inst]) & (summary.msIonisation == 'nan')), 'msIonisation'] = 'APCI'
    summary.loc[(np.array([("appi" in x and not "dappi" in x) for x in summary.GNPS_Inst]) & (summary.msIonisation == 'nan')), 'msIonisation'] = 'APPI'
    summary.loc[(np.array(["dappi" in x for x in summary.GNPS_Inst]) & (summary.msIonisation == 'nan')), 'msIonisation'] = 'DAPPI'

    # Mass Analyzer (Not Done)
    summary.loc[(np.array(["orbitrap" in x for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == 'nan')),"msMassAnalyzer"] = "orbitrap"
    summary.loc[(np.array([("quadrupole tof" in x or "qtof" in x or "q-tof" in x) and not "qq" in x for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == 'nan')),"msMassAnalyzer"] = "qtof"
    summary.loc[(np.array([("tof" in x) and not ("qq" in x or "qtof" in x or "q-tof" in x or "q tof" in x or "quadrupole tof" in x) for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == 'nan')),"msMassAnalyzer"] = "tof"
    summary.loc[(np.array([("qft" in x) for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == 'nan')),"msMassAnalyzer"] = "quadurpole"
    summary.loc[(np.array([("ion trap" in x) or ('itms' in x) for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == 'nan')),"msMassAnalyzer"] = "ion trap"
    summary.loc[(np.array([("itft" in x) or ("fticr" in x) for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == 'nan')),"msMassAnalyzer"] = "ftms"
    summary.loc[(np.array([("lc-esi-q" == x) or ("lc-appi-qq" == x) or ("lcq" in x) or ("qqq" in x ) for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == 'nan')),"msMassAnalyzer"] = "quadurpole"
    summary.loc[(np.array([("impact hd" == x) for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == 'nan')),"msMassAnalyzer"] = "qtof"
  
    # Manufacturer Info (Not Done)
    summary.loc[(np.array([bool("maxis" in x) for x in summary.GNPS_Inst]) & (summary.msManufacturer == "nan")),"msManufacturer"] = "Bruker Daltonics"
    summary.loc[(np.array(["q exactive" in x or "q-exactive" in x for x in summary.GNPS_Inst]) & (summary.msManufacturer == "nan")),"msManufacturer"] = "Thermo"
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
    test_df = summary[(summary.msManufacturer == 'Thermo') & (summary.msMassAnalyzer == 'qtof')]
    if len(test_df) != 0:
        pd.options.display.max_columns = None
        print(test_df.head(10))
        raise ValueError("There are {} entries with Thermo and qtof".format(len(test_df)))
    test_df = summary[(summary.msMassAnalyzer == 'orbitrap') & (summary.msManufacturer == 'Bruker Daltonics')]
    if len(test_df) != 0:
        pd.options.display.max_columns = None
        print(test_df.head(10))
        raise ValueError("There are {} entries with Bruker Daltonics and orbitrap".format(len(test_df)))
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
    summary.scan = np.arange(1, len(summary)+ 1)
        
    sanity_checks(summary)

        
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