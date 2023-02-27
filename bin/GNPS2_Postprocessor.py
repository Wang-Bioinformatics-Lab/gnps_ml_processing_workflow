import os
import pandas as pd
from pyteomics import mgf
from pyteomics.mgf import IndexedMGF
from glob import glob
import re
from pathlib import Path 
import datetime
from rdkit.Chem import AllChem
from rdkit import Chem

def sanity_checks(summary):
    assert len(summary[(summary.msManufacturer == 'Thermo') & (summary.msMassAnalyzer == 'qtof')]) == 0

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
    pattern = re.compile("[^\D]([0-9]*[.]{0,1}[0-9]*)")
    summary.loc[~ summary.retention_time.isna(),"retention_time"] = summary.loc[~ summary.retention_time.isna(),"retention_time"].apply(lambda x: re.search(pattern,x).group(0))

    # Adduct
    # Strip starting and ending braces if no charge is specified
    summary.Adduct = summary.Adduct.map(lambda x: x.strip()[1:-1] if x.strip()[-1] == ']' and x.strip()[0] == '[' else x.strip())

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
    summary.msMassAnalyzer = summary.astype(str).msMassAnalyzer.map(lambda x: x.lower())
    summary.loc[summary.msMassAnalyzer == 'quadrupole tof','msMassAnalyzer'] = 'qtof'
    summary.loc[summary.msMassAnalyzer == 'fourier transform ion cyclotron resonance mass spectrometer','msMassAnalyzer'] = 'ftms'

    # Detector
    '''A very specific set of files has MS:1000253 as the detector name which is used for mzML files, 
    however, it's written in mzXML files. We'll clean this up here.

    Conversion source: https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo
    '''

    summary.loc[['MS:1000253' == x for x in summary.msDetector], 'msDetector'] = 'electron multiplier tube'
    
    return summary

def propogate_GNPS_Inst_field(summary):
    summary.GNPS_Inst = summary.GNPS_Inst.astype('str')
    summary.GNPS_Inst = summary.GNPS_Inst.map(lambda x: x.strip().lower())

    # Fragmentation Info (Done)
    summary.msDissociationMethod = summary.msDissociationMethod.astype(str)
    summary.loc[[("in source cid" == x) for x in summary.GNPS_Inst] and summary.msDissociationMethod == 'nan', 'msDissociationMethod'] = "is-cid"    
    summary.loc[[("hid" in x) for x in summary.GNPS_Inst] and summary.msDissociationMethod == 'nan', 'msDissociationMethod'] = "hid"    
    summary.loc[[("cid" in x and not "is-cid" in x) for x in summary.GNPS_Inst] and summary.msDissociationMethod == 'nan', 'msDissociationMethod'] = "cid"    

    # Ionisation Info (Not Done)
    summary.loc[["esi" in x for x in  summary.GNPS_Inst] and summary.msIonisation == 'nan', 'msIonisation'] = 'ESI'
    summary.loc[["apci" in x for x in  summary.GNPS_Inst] and summary.msIonisation == 'nan', 'msIonisation'] = 'APCI'
    summary.loc[[("appi" in x and not "dappi" in x) for x in summary.GNPS_Inst] and summary.msIonisation == 'nan', 'msIonisation'] = 'APPI'
    summary.loc[["dappi" in x for x in summary.GNPS_Inst] and summary.msIonisation == 'nan', 'msIonisation'] = 'DAPPI'

    # Mass Analyzer (Not Done)
    summary.loc[["orbitrap" in x for x in summary.GNPS_Inst] and summary.msMassAnalyzer == 'nan',"msMassAnalyzer"] = "orbitrap"
    summary.loc[[("quadrupole tof" in x or "qtof" in x or "q-tof" in x) and not "qq" in x for x in summary.GNPS_Inst] and summary.msMassAnalyzer == 'nan',"msMassAnalyzer"] = "qtof"
    summary.loc[[("tof" in x) and not ("qq" in x or "qtof" in x or "q-tof" in x or "q tof" in x or "quadrupole tof" in x) for x in summary.GNPS_Inst] and summary.msMassAnalyzer == 'nan',"msMassAnalyzer"] = "tof"

    # Manufacturer Info (Not Done)
    summary.loc[["maxis" in x for x in summary.GNPS_Inst] and summary.msManufacturer == "nan","msManufacturer"] = "Bruker Daltonics"
    summary.loc[["q exactive" in x or "q-exactive" in x for x in summary.GNPS_Inst] and summary.msManufacturer == "nan","msManufacturer"] = "Thermo"
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
    summary.drop('mol', axis=1)
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
    indexed_mgf = IndexedMGF(input_mgf)
    level_0 = 0
    for m in indexed_mgf:
        spectrum_id = m['params']['title']
        if spectrum_id in spectrum_ids: # Make sure that it didn't get removed during cleaning
            mz_array = m['m/z array']
            intensity_array = m['intensity array']
            precursor_mz = m['params']['pepmass']
            # charge = m['charge']
            for index, (mz, intensity) in enumerate(zip(mz_array, intensity_array)):
                output.append({'spectrum_id':spectrum_id, 'level_0': level_0, 'index':index, 'i':intensity, 'mz':mz, 'prec_mz':precursor_mz})
                level_0 += 1
                
    output = pd.DataFrame(output)
    output.set_index('spectrum_id')
    return output
            
                       

def postprocess_files(csv_path, mgf_path, output_csv_path, output_parquet_path):
    if not os.path.isfile(csv_path):
        if not os.path.isfile(mgf_path):
            file_pattern = re.compile(r'.*?(\d+).*?')
            def get_order(file,):
                match = file_pattern.match(Path(file).name)
                return int(match.groups()[0])

            sorted_csv_files = sorted(glob('./temp/temp_*.csv'), key=get_order)
            sorted_mgf_files = sorted(glob('./temp/temp_*.mgf'), key=get_order)

            os.system("cat " + " ".join(sorted_csv_files) +"> " + csv_path)
            os.system("cat " + " ".join(sorted_mgf_files) +"> " + mgf_path)
    
    # TODO: Delete temp files

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

    parquet_as_df = generate_parquet_file(mgf_path, summary.spectrum_id)
    parquet_as_df.to_parquet(output_parquet_path)
    summary.to_csv(output_csv_path, index=False)

def main():
    now = datetime.datetime.now()
    year = now.year
    quarter = int(now.month/4) + 1

    csv_path = "./GNPS_ml_exports/ALL_GNPS_merged_{}_{}.csv".format(quarter, year)
    mgf_path = "./GNPS_ml_exports/ALL_GNPS_merged_{}_{}.mgf".format(quarter, year)
    cleaned_csv_path = "./GNPS_ml_exports/ALL_GNPS_cleaned_{}_{}.csv".format(quarter, year)
    cleaned_parquet_path = "./GNPS_ml_exports/ALL_GNPS_cleaned_{}_{}.parquet".format(quarter, year)
    if not os.path.isdir('./GNPS_ml_exports/'):
        os.makedirs('./GNPS_ml_exports/')


    if not os.path.isfile(cleaned_csv_path):
        if not os.path.isfile(cleaned_parquet_path):
            postprocess_files(csv_path, mgf_path, cleaned_csv_path, cleaned_parquet_path)
            
if __name__ == '__main__':
    main()