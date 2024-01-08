import ast
import math
import os
import pickle
import numpy as np
import pandas as pd
from pyteomics.mgf import IndexedMGF
import re
from tqdm import tqdm
from utils import harmonize_smiles_rdkit, neutralize_atoms, INCHI_to_SMILES, synchronize_spectra, generate_parquet_df
from rdkit import Chem
from rdkit.Chem import Descriptors
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

    # Precursor MZ
    # Drop all entries with precursor mz <= 1
    summary.Precursor_MZ = summary.Precursor_MZ.astype('float')
    summary = summary.loc[summary.Precursor_MZ > 1]

    # smiles
    summary.Smiles = summary.Smiles.astype(str).parallel_apply(lambda x: x.strip() )
    summary.Smiles = summary.Smiles.parallel_apply(lambda x: '' if ('N/A' in x) or ('nan' in x) else x)
    
    # check smiles validity
    summary.Smiles = summary.Smiles.apply(lambda x: x if Chem.MolFromSmiles(x) is not None else '')
    
    # INCHI
    summary.INCHI = summary.INCHI.astype(str).parallel_apply(lambda x: x.strip().replace('"', '') )
    summary.INCHI = summary.INCHI.parallel_apply(lambda x: '' if ('N/A' in x) or ('nan' in x) else x)
    
    # InChIKey_SMILES
    summary.InChIKey_smiles = summary.InChIKey_smiles.astype(str).parallel_apply(lambda x: x.strip() )
    summary.InChIKey_smiles = summary.InChIKey_smiles.parallel_apply(lambda x: '' if ('N/A' in str(x)) else x)
    
    # InChIKey_InChI
    summary.InChIKey_inchi = summary.InChIKey_inchi.astype(str).parallel_apply(lambda x: x.strip() )
    summary.InChIKey_inchi = summary.InChIKey_inchi.parallel_apply(lambda x: '' if ('N/A' in str(x)) else x)

    # ionization
    summary.msIonisation = summary.msIonisation.astype(str)
    summary.loc[['esi' in x.lower() or 'electrospray' in x.lower() for x in summary.msIonisation], 'msIonisation'] = "ESI"
    summary.loc[['' == x or 'positive' == x.lower() or 'negative'  == x.lower() for x in summary.msIonisation], 'msIonisation'] = ""
    summary.loc[['LC-APCI' in x.upper() for x in summary.msIonisation], 'msIonisation'] = "APCI"   
    
    # Compound Source
    summary.Compound_Source = summary.Compound_Source.astype(str)
    summary.Compound_Source = summary.Compound_Source.parallel_apply(lambda x: x.strip().lower())
    summary.loc[[('unknown' == x) or ('nan' == x) or ('other' == x) or ('lcms'==x) for x in summary.Compound_Source], 'Compound_Source'] = ''
    summary.loc[[('commercial standard' == x) or ('commercial' == x) or ('standard' == x) or ('cosmetic _raw meterial' == x) for x in summary.Compound_Source], 'Compound_Source'] = "commercial"
    summary.loc[[('prestwick' in x) or ('nih natural product library' == x)  or('nih pharmacologically active library' == x) for x in summary.Compound_Source], 'Compound_Source'] = "isolated"
    summary.loc[[('lysate' == x)  for x in summary.Compound_Source], 'Compound_Source'] = "crude"
    
    
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
            
            # This will drop adducts that are not mapped
            adduct = mapped_adduct
            
        # Add 1 if not None and last is
        if adduct is not None:
            if adduct[-2:] == "]+":
                adduct = adduct[:-2] + "]1+"
            elif adduct[-2:] == "]-":
                adduct = adduct[:-2] + "]1-"
            
        return adduct
        
    summary.loc[:, 'Adduct'] = summary.Adduct.parallel_apply(helper)
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
    summary.Charge = summary.Charge.fillna(0).astype(int)
    
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
        
    # Check if it's just in the collision energy field and needs cleaning
    pattern = re.compile(r'(\d+)')
    extracted_eV = summary.collision_energy.apply(lambda x: re.search(pattern, str(x)))
    summary.collision_energy = extracted_eV.apply(lambda x: x.group(1) if x is not None else None).astype(float)
    
    # Ion Mode
    summary.Ion_Mode = summary.Ion_Mode.parallel_apply(lambda x: '' if ('n/a' in x) or ('nan' in x) or ('unknown' in x) else x)
    summary.loc[summary.Ion_Mode == 'positive-20ev','Ion_Mode'] = 'positive'
    # We'll infer any missing ion modes from the charge
    mask = (summary.Ion_Mode == '') & (summary.Charge > 0)
    if sum(mask) > 0:
        print(f"Imputing {sum(mask)} positive ion modes using the Charge field.")
        summary.loc[mask, 'Ion_Mode'] = 'positive'
    mask = (summary.Ion_Mode == '') & (summary.Charge < 0)
    if sum(mask) > 0:
        print(f"Imputing {sum(mask)} negative ion modes using the Charge field.")
        summary.loc[mask, 'Ion_Mode'] = 'negative'
    # Any ion modes that disagree with the charge field will be replaced
    mask = (summary.Ion_Mode != 'negative') & (summary.Charge < 0)
    if sum(mask) > 0:
        print(f"Correcting {sum(mask)} positive ion modes using the Charge field.")
        summary.loc[mask, 'Ion_Mode'] = 'negative'
    mask = (summary.Ion_Mode != 'positive') & (summary.Charge > 0)
    if sum(mask) > 0:
        print(f"Correcting {sum(mask)} negative ion modes using the Charge field.")
        summary.loc[mask, 'Ion_Mode'] = 'positive'

    # Manufacturer
    summary.msManufacturer = summary.msManufacturer.astype('str')
    summary.loc[['thermo' in x.lower() for x in summary.msManufacturer], 'msManufacturer'] = 'Thermo'
    summary.msManufacturer = summary.msManufacturer.parallel_apply(lambda x: '' if ('n/a' in x) or ('nan' in x) else x)

    # msMassAnalyzer
    # This cleanup address lists of mass analyzers from mzML and mzXML files
    def transform_analyzer(x:str):
        x = x.lower()
        if 'quadrupole tof' == x:
            return 'qtof'
        if 'fourier transform ion cyclotron resonance mass spectrometer' == x:
            return 'ftms'
        if 'time-of-flight' == x or 'it-tof' == x:
            return 'tof'
        if 'ion trap' in x or 'itms' in x or 'lcq' in x:
            return 'ion trap'
        return x
    def merge_analyzer(x:list):
        if np.any(["orbitrap" in y for y in x]):
            return "orbitrap"
        if "quadrupole" in x and "tof" in x:
            return 'qtof'
        return x[0]
    
    summary.msMassAnalyzer = summary.msMassAnalyzer.astype('str')
    summary.msMassAnalyzer = summary.msMassAnalyzer.str.lower()
    summary.msMassAnalyzer = summary.msMassAnalyzer.parallel_apply(lambda x: '' if ('n/a' in x) or ('nan' in x) or ('unknown' in x) else x)
    
    # msModel
    summary.msModel = summary.msModel.astype('str')
    summary.msModel = summary.msModel.str.lower()
    summary.msModel = summary.msModel.parallel_apply(lambda x: '' if ('n/a' in x) or ('nan' in x) or ('unknown' in x) else x)
    
    # msDissociationMethod
    summary.msDissociationMethod = summary.msDissociationMethod.astype('str')
    summary.msDissociationMethod = summary.msDissociationMethod.str.lower()
    summary.loc[:,'msDissociationMethod'] = summary.msDissociationMethod.parallel_apply(lambda x: '' if ('n/a' in x) or ('nan' in x) or ('unknown' in x) else x)
    summary.loc[(np.array(['beam-type' in x for x in summary.msDissociationMethod])), 'msDissociationMethod'] = 'hcd'
    summary.loc[(np.array(['collision-induced' in x or 'low-energy cid' in x for x in summary.msDissociationMethod])), 'msDissociationMethod'] = 'cid'
    
    
    mask = (summary.msMassAnalyzer.notna()) & (summary.msMassAnalyzer != 'nan')
    def literal_eval_helper(x):
        """
        An small helper function to handle weird entries in msMassAnalyzer
        """
        if x == '':
            return ['']
        try:
            return ast.literal_eval(x)
        except Exception as e:
            print(e)
            print(f"Error in literal_eval_helper when trying to parse {x}")
            return ['']
    
    summary.loc[mask,'msMassAnalyzer'] = summary.loc[mask,'msMassAnalyzer'].apply(literal_eval_helper)
    summary.loc[mask,'msMassAnalyzer'] = summary.loc[mask,'msMassAnalyzer'].apply(lambda x: [transform_analyzer(y) for y in x])
    summary.loc[mask,'msMassAnalyzer'] = summary.loc[mask,'msMassAnalyzer'].apply(merge_analyzer)
    # summary.loc[mask].apply(lambda x: [y == 'quadrupole tof' for y in x].any())
    # summary.loc[mask,'msMassAnalyzer'] = 'qtof' 
    # summary.loc[summary.msMassAnalyzer == 'fourier transform ion cyclotron resonance mass spectrometer','msMassAnalyzer'] = 'ftms'

    # msIonisation
    summary.msIonisation = summary.msIonisation.astype('str')
    summary.msIonisation = summary.msIonisation.str.upper()
    summary.msIonisation = summary.msIonisation.parallel_apply(lambda x: '' if ('N/A' in x) or ('NAN' in x) else x) 

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
    
    # Check INCHI is parsable
    summary.loc[:, 'INCHI'] = summary.loc[:,'INCHI'].apply(lambda x: '' if Chem.inchi.MolFromInchi(x) is None else x if x != '' else '')
       
    # In rare cases the user will use INCHI and not smiles, so we'll convert it to smiles
    mask = (summary.Smiles == '') & (summary.INCHI != '')
    summary.loc[mask, 'Smiles'] = summary.loc[mask, 'INCHI'].apply(INCHI_to_SMILES)
    
    summary.Smiles = summary.Smiles.parallel_apply(lambda x: harmonize_smiles_rdkit(x, skip_tautomerization=True))
    
    # Check the INCHI and SMILES are equivalent
    mask = (summary.Smiles != '') & (summary.INCHI != '')
    inchi_from_smiles = summary.loc[mask].apply(lambda x: Chem.inchi.MolToInchi(Chem.MolFromSmiles(x['Smiles'])), axis=1)
    equivalency_mask = (inchi_from_smiles != summary.loc[mask, 'INCHI'])
    if sum(equivalency_mask) > 0 :
        print(f"Warning: {sum(equivalency_mask)} entries have INCHI and SMILES that are not equivalent, SMILES values will be used to replace INCHI")
        summary.loc[mask & equivalency_mask, 'INCHI'] = inchi_from_smiles.loc[equivalency_mask]
    
    # If no INCHI but we have SMILES, convert it to INCHI
    mask = (summary.Smiles != '') & (summary.INCHI == '')
    summary.loc[mask, 'INCHI'] = summary.loc[mask, 'Smiles'].apply(lambda x: Chem.inchi.MolToInchi(Chem.MolFromSmiles(x)))
    # Fill in INCHI key
    # Generate all inchi keys
    all_keys = summary.loc[:, 'INCHI'].apply(lambda x: Chem.inchi.InchiToInchiKey(x) if x != '' else '')
    # Check existing keys
    incorrect_mask = (summary.InChIKey_smiles != all_keys)
    if sum(incorrect_mask) > 0 :
        print(f"Warning: {sum(incorrect_mask)} entries have INCHI and INCHIKey that are not equivalent, new INCHIKeys will be generated based on INCHI")
        summary.loc[incorrect_mask, 'InChIKey_smiles'] = summary.loc[incorrect_mask, 'INCHI'].apply(lambda x: Chem.inchi.InchiToInchiKey(x))
       
    return summary

def validate_monoisotopic_masses(summary:pd.DataFrame):
    """This function takes a pandas dataframe with a 'Smiles' column and validates that the ExactMass column is correct.
        If the ExactMass column is incorrect, the function will attempt to correct it, printing a message indicating how many
        masses were incorrect.

    Args:
        summary (pd.DataFrame): A summary dataframe containing a 'Smiles' column, and a 'ExactMass' column
        
    Returns:
        pd.DataFrame: The modified summary dataframe
    """
    parsable_mask = summary.Smiles.apply(lambda x: Chem.MolFromSmiles(x) is not None if x != '' else False)
    # Anything uparsable should have an ExactMass of np.nan
    summary.loc[~parsable_mask, 'ExactMass'] = np.nan
    
    correct_masses = summary.loc[parsable_mask, 'Smiles'].parallel_apply(lambda x: Descriptors.ExactMolWt(Chem.MolFromSmiles(x)))
    correct_mask = summary.loc[parsable_mask, 'ExactMass'] != correct_masses
    
    if sum(correct_mask) > 0:
        print(f"Warning: {sum(correct_mask)} entries have ExactMasses that are not equivalent to the monoisotopic mass of the SMILES.")
        print(f"Of the incorrect masses, {sum(summary.loc[parsable_mask, 'ExactMass'].loc[correct_mask]==0)} entries have an ExactMass of 0")
        print("ExactMasses will be replaced with the monoisotopic mass of the SMILES")
        diff = summary.loc[parsable_mask & correct_mask, 'ExactMass'] - correct_masses.loc[correct_mask]
        # Save a csv of the difference in incorrect masses
        diff.to_csv('./incorrect_masses.csv')
        # Save a csv in the difference of incorrect masses that were not zero
        diff =  summary.loc[parsable_mask & (correct_mask & (summary.ExactMass !=0)), 'ExactMass'] - correct_masses.loc[correct_mask & (summary.ExactMass !=0)]
        diff.to_csv('./incorrect_masses_non_zero.csv')
        
        
        summary.loc[parsable_mask & correct_mask, 'ExactMass'] = correct_masses.loc[correct_mask]
    
    return summary

def check_M_H_adducts(summary):
    # Verify that there are no protonated or unprotonated smiles strings
    mask = (summary.Smiles != '')
    result = summary.loc[mask, 'Smiles'].parallel_apply(neutralize_atoms)   # Returns num_removed_charges, pos_and_neg, sum_of_charges, smiles
    
    num_removed_charges = result.apply(lambda x: x[0]).astype(int)
    net_charge = result.apply(lambda x: x[2]).astype(int)
    smiles = result.apply(lambda x: x[3]).astype(str)
    
    # Check for molecules with both positive and negative charges that can be removed with protonation/deprotonation
    # that sum to zero
    net_zero = result.apply(lambda x: x[1]) & (net_charge==0)
    
    # If these exist, they must have an adduct that is not [M], otherwise they could not be detected
    adduct_is_M = summary.loc[(mask & net_zero), 'Adduct'].apply(lambda x: bool('[M]' in str(x)))
    to_drop = (mask & net_zero & adduct_is_M)
    print(f"There are {sum(to_drop)} molecules with with a net charge of 0. And [M] Adduct. These will be dropped.")
    summary = summary.loc[~to_drop]
    
    # # Mask for any charges that have changed.
    # removed_charge_mask = num_removed_charges & (~ to_drop)
    # if sum(removed_charge_mask) > 0:
    #     print(f"Found {sum(removed_charge_mask)} structures with non-intrinsic charges, these will be updated. Adducts and charges will not be updated.")
    #     summary.loc[mask & removed_charge_mask, 'Smiles'] = smiles.loc[removed_charge_mask]
        
        
    #     # DEBUG
    #     summary.loc[mask & removed_charge_mask].to_csv('./structures_with_non_intrinsic_charges.csv')
    #     #### So there's a decision point here, we could update the charges on the adducts and in the charge column
    #     #### but that feels risky. I think the best way to go here is to trust that they got the overall charge + adduct right
    #     #### definitely double check this thought though
        
    #     # Recalculate the ExactMass. We don't care about the adduct mass because this will be used to check if there is no adduct
    #     summary.loc[mask & removed_charge_mask, 'ExactMass'] = summary.loc[mask & removed_charge_mask, 'Smiles'].parallel_apply(lambda x: Descriptors.ExactMolWt(Chem.MolFromSmiles(x)))
    
    # it's very hard to tell if the original molecule was charged or not, and whether the adduct is supposed to make up for it 
    # ex: CCCCCCCCC=CCCCCCCCC(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)O with [M+H]2+, 
    #     was it CCCCCCCCC=CCCCCCCCC(=O)OCC(O)COP(=O)(O)OCC[N+](C)(C)C with no adduct ([M]+)?
    #     or did it enter the machine with a positive charge?
    #     You can look at the adduct and see if there are hydrogens to make up for the charge on the element, but this doesn't feel very concrete.
    #     Here, we would be able to note that the instrument was in positive ion mode, and this is a negative charge, so it's likely that the molecule
    #     the molecule was not really charged.
    # In a future release we'll consider updating the adducts and charges based on the new smiles.
    
    # Mask all charges that have changed, if there is no remaining intrinsic charge, we'll drop the whole row because
    # to_drop_2 = (net_charge != 0) & (num_removed_charges > 0)
    # summary = summary.loc[~to_drop_2]

        
    # For mols with intrinsic charges, if the Precursor_MZ and ExactMass are approximately the same, we will change the adduct to [M(+/-)xe]x(+/-)
    intrinsic_charge_mask = (net_charge != 0) & (~ (to_drop))
    mass_diff = summary.loc[intrinsic_charge_mask, 'Precursor_MZ'] - summary.loc[intrinsic_charge_mask, 'ExactMass']
    mass_mask = mass_diff.abs() <= 0.01
    
    suspect_adduct_mask = (intrinsic_charge_mask & mass_mask)
    
    # Try rewriting all adducts
    corrected_adducts = net_charge.loc[suspect_adduct_mask].apply(lambda x: f'[M]{x}+' if x > 0 else f'[M]{x}-')
    # Check which ones are incorrect
    # Here, we '&' with mask, because we need to align the boolean indexers. This is just as good with anding with True for all rows
    incorrect_adduct_mask = (summary.loc[suspect_adduct_mask, 'Adduct'] != corrected_adducts) & mask   
    
    if sum(incorrect_adduct_mask) > 0:
        print(f"Found {sum(incorrect_adduct_mask)} structures with incorrect adducts, these will be changed to [M]+")
        summary.loc[incorrect_adduct_mask, 'Adduct'] = corrected_adducts.loc[incorrect_adduct_mask]
        
    # Drop everything we didn't update
    # i.e. keep everything that we corrected the adduct of, didn't have SMILES, or we didn't change the charge of
    summary = summary.loc[(incorrect_adduct_mask) | (~mask) | (num_removed_charges == 0)]
    
    return summary

def propagate_GNPS_Inst_field(summary):
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
    
    # Fragmentation Info (Done)
    summary.msDissociationMethod = summary.msDissociationMethod.astype(str)
    summary.loc[(np.array(["in source cid" == x for x in summary.GNPS_Inst]) & (summary.msDissociationMethod == '')), 'msDissociationMethod'] = "is-cid"
    summary.loc[(np.array([("hid" in x) for x in summary.GNPS_Inst]) & (summary.msDissociationMethod == '')), 'msDissociationMethod'] = "hid"
    summary.loc[(np.array([("hcd" in x) for x in summary.GNPS_Inst]) & (summary.msDissociationMethod == '')), 'msDissociationMethod'] = "hcd"
    summary.loc[(np.array([("cid" in x and not "is-cid" in x) for x in summary.GNPS_Inst]) & (summary.msDissociationMethod == '')), 'msDissociationMethod'] = "cid"

    # Ionisation Info (Not Done)
    summary.loc[(np.array(["FAB" in x for x in  summary.GNPS_Inst]) & (summary.msIonisation == '')), 'msIonisation'] = 'FAB'
    summary.loc[(np.array(["ESI" in x for x in  summary.GNPS_Inst]) & (summary.msIonisation == '')), 'msIonisation'] = 'ESI'
    summary.loc[(np.array(["APCI" in x for x in  summary.GNPS_Inst]) & (summary.msIonisation == '')), 'msIonisation'] = 'APCI'
    summary.loc[(np.array(["API" in x for x in  summary.GNPS_Inst]) & (summary.msIonisation == '')), 'msIonisation'] = 'API'
    summary.loc[(np.array([("APPI" in x and not "DAPPI" in x) for x in summary.GNPS_Inst]) & (summary.msIonisation == '')), 'msIonisation'] = 'APPI'
    summary.loc[(np.array(["DAPPI" in x for x in summary.GNPS_Inst]) & (summary.msIonisation == '')), 'msIonisation'] = 'DAPPI'

    # Mass Analyzer (Not Done)
    summary.loc[(np.array([("orbitrap" in x) or ("q-exactive" in x) for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == '')),"msMassAnalyzer"] = "orbitrap"
    summary.loc[(np.array(["maxis" in x for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == '')),"msMassAnalyzer"] = "qtof"
    summary.loc[(np.array([("quadrupole tof" in x or "qtof" in x or "q-tof" in x) and not "qq" in x for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == '')),"msMassAnalyzer"] = "qtof"
    summary.loc[(np.array([("tof" in x) and not ("qq" in x or "qtof" in x or "q-tof" in x or "q tof" in x or "quadrupole tof" in x) for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == '')),"msMassAnalyzer"] = "tof"
    summary.loc[(np.array([("itft" in x) or ("it-ft" in x) or ("qft" in x) or ("hybrid ft" in x) for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == '')),"msMassAnalyzer"] = "ftms"
    summary.loc[(np.array([("fticr" in x) or ("ft-icr" in x) for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == '')),"msMassAnalyzer"] = "fticr"
    summary.loc[(np.array([("lc-esi-q" == x) or ("lc-esi-qq" == x) or ("lc-appi-qq" == x) or ("qqq" in x ) or ("beqq" in x) for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == '')),"msMassAnalyzer"] = "quadrupole"
    summary.loc[(np.array([("impact hd" == x) for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == '')),"msMassAnalyzer"] = "qtof"
    summary.loc[(np.array([("ion trap" in x) or ('itms' in x) or ("lcq" in x) or ("qit" in x)  or ("lit" in x) or ("lc-esi-it" in x) for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == '')),"msMassAnalyzer"] = "ion trap"
    summary.loc[(np.array([("cid" in x) and ('lumos' in x) for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == '')),"msMassAnalyzer"] = "ion trap"
    summary.loc[(np.array([("hcd" in x) and ('lumos' in x) for x in summary.GNPS_Inst]) & (summary.msMassAnalyzer == '')),"msMassAnalyzer"] = "orbitrap"
  
    # Manufacturer Info (Not Done)
    summary.loc[(np.array([bool("maxis" in x) for x in summary.GNPS_Inst]) & (summary.msManufacturer == "")),"msManufacturer"] = "Bruker Daltonics"
    summary.loc[(np.array([("orbitrap" in x) for x in summary.GNPS_Inst]) & (summary.msManufacturer == "")),"msManufacturer"] = "Thermo"
    summary.loc[(np.array(["q exactive" in x or "q-exactive" in x for x in summary.GNPS_Inst]) & (summary.msManufacturer == "")),"msManufacturer"] = "Thermo"
    return summary

def propagate_msModel_field(summary):
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
        raise ValueError(f"There are {len(test_df)} entries with Thermo and qtof")

    test_df = summary[(summary.msMassAnalyzer == 'orbitrap') & ((summary.msManufacturer != 'Thermo') & (summary.msManufacturer != ''))]
    if len(test_df) != 0:
        pd.options.display.max_columns = None
        print(test_df.head(10))
        raise ValueError(f"There are {len(test_df)} entries with orbitrap that are not Thermo.")
    
    # assert len(summary[(summary.Adduct == 'None') & (summary.Adduct == 'nan') & (summary.Adduct.isna())]) == 0 # Right now because of the adduct cleaning, there is a chance that we'll have nan adducts


def add_columns_formula_analysis(summary): 
    column_name_ppmBetweenExpAndThMass='ppmBetweenExpAndThMass'
    
    def helper(row):
        try:
            smiles = str(row['Smiles'])
            if row['Smiles'] != '':
                # For molecules with intrinsic charges, the adduct is [M]+ or [M]-, so we'll remove it
                # This prevents the loss of charge from being double counted with the following +/- and on the molecule, 'M' itself
                adduct = row['Adduct']
                if adduct == '[M]+' or adduct == '[M]-':
                    adduct = '[M]'
                
                formula = Formula.formula_from_smiles(smiles, adduct, no_api=False)  # Disabling API can improve speed 
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
            if smiles != '':
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
            
def postprocess_files(csv_path, mgf_path, output_csv_path, output_parquet_path, cleaned_mgf_path, includes_massbank=False):
    pandarallel.initialize(progress_bar=False, nb_workers=PARALLEL_WORKERS, use_memory_fs = False)
    
    summary = pd.read_csv(csv_path)
    
    # If the merged files include massbank, drop the old massbank data
    if includes_massbank:
        summary = summary.loc[(summary.GNPS_library_membership != 'MASSBANK') & (summary.GNPS_library_membership != 'MASSBANKEU')]

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
       
    # Clean up monoistopic masses:
    print("Cleaning up monoisotopic masses", flush=True)
    start = time.time()
    summary = validate_monoisotopic_masses(summary)
    print("Done in {} seconds".format(datetime.timedelta(seconds=time.time() - start)), flush=True)
    
    # Check M+H adducts should not be [M]+:
    print("Checking for missing [M]+ Adducts")
    start = time.time()
    summary = check_M_H_adducts(summary)
    print("Done in {} seconds".format(datetime.timedelta(seconds=time.time() - start)), flush=True)
    
    # Exploiting GNPS_Inst annotations:
    print("Attempting to propagate user instrument annotations", flush=True)
    start = time.time()
    summary = propagate_GNPS_Inst_field(summary)
    print("Done in {} seconds".format(datetime.timedelta(seconds=time.time() - start)), flush=True)

    # Exploiting Some of the info in msModel
    print("Attempting to propagate msModel field", flush=True)
    start = time.time()
    summary = propagate_msModel_field(summary)
    print("Done in {} seconds".format(datetime.timedelta(seconds=time.time() - start)), flush=True)

    # Calculating ppm error
    print("Calculating ppm error", flush=True)
    start = time.time()
    add_columns_formula_analysis(summary)
    print("Done in {} seconds".format(datetime.timedelta(seconds=time.time() - start)), flush=True)
    
    # Cleanup scan numbers, scan numbers will be reset in mgf by synchronize spectra
    summary.scan = np.arange(1, len(summary)+ 1)
        
    sanity_checks(summary)
        
    print("Writing csv file", flush=True)
    summary = summary.drop(columns=['retention_time', 'msDetector', 'msModel', 'GNPS_Inst', 'InChIKey_inchi'])
    summary.to_csv(output_csv_path, index=False)
        
    # Cleanup MGF file. Must be done before explained intensity calculations in order to make sure spectra are in order
    print("Writing mgf file", flush=True)
    start = time.time()
    synchronize_spectra(mgf_path, cleaned_mgf_path, summary)
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
    parquet_as_df.to_parquet(output_parquet_path, index=False)
    print("Postprocessing complete!", flush=True)

def main():
    parser = argparse.ArgumentParser(description='Postprocess GNPS files')
    parser.add_argument('--input_csv_path', type=str, default="ALL_GNPS_merged.csv", help='Path to the csv file')
    parser.add_argument('--input_mgf_path', type=str, default="ALL_GNPS_merged.mgf", help='Path to the mgf file')
    parser.add_argument('--output_csv_path', type=str, default="ALL_GNPS_cleaned.csv", help='Path to the output csv file')
    parser.add_argument('--output_parquet_path', type=str, default="ALL_GNPS_cleaned.parquet", help='Path to the output parquet file')
    parser.add_argument('--output_mgf_path', type=str, default="ALL_GNPS_cleaned.mgf", help='Path to the output mgf file')
    parser.add_argument('--includes_massbank', action='store_true', help='Whether the merged files include reparsed massbank entries.')
    args= parser.parse_args()
    
    csv_path                = str(args.input_csv_path)
    mgf_path                = str(args.input_mgf_path)
    cleaned_csv_path        = str(args.output_csv_path)
    cleaned_parquet_path    = str(args.output_parquet_path)
    cleaned_mgf_path        = str(args.output_mgf_path)

    if not os.path.isfile(cleaned_csv_path):
        if not os.path.isfile(cleaned_parquet_path):
            postprocess_files(csv_path, mgf_path, cleaned_csv_path, cleaned_parquet_path, cleaned_mgf_path, args.includes_massbank)
            
if __name__ == '__main__':
    main()
    
    
def test_propagate_GNPS_Inst_field():
    # pip install -U pytest
    # Command to run: python -m pytest ./GNPS2_Postprocessor.py
    
    # This dataframe contains a manually annotated version of the unique GNPS_Inst values
    # We should expect that the logic in our code exactly matches the manual annotation
    test_df = pd.read_csv('../test_data/unique_GNPS_Inst_values.csv', dtype=str)
    # Make all NaN values ''
    test_df = test_df.fillna('')
       
    out = propagate_GNPS_Inst_field(test_df)
    correct_cols = [col for col in out.columns if col.startswith('correct_')]
    inconsistent_rows = pd.DataFrame(columns=out.columns)

    for index, row in out.iterrows():
        for col in correct_cols:
            partner_col = col.replace('correct_', '')
            if str(row[col]) != str(row[partner_col]):
                print(f"'{str(row[col])}', '{str(row[partner_col])}'")
                inconsistent_rows = inconsistent_rows._append(row, ignore_index=True)
                break
    
    if len(inconsistent_rows) > 0:
        print("The following rows are inconsistent:")
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)
        print(inconsistent_rows)
        inconsistent_rows.to_csv('./GNPS_inst_test_output.csv')        
        raise AssertionError("The GNPS_Inst field is not being propagated correctly")