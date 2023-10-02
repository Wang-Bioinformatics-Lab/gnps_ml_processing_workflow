import argparse
import csv
import os
import re

from rdkit import Chem
from rdkit.Chem import inchi


def process_file(input_path, output_path = None):
    # Initialize an empty dictionary to store the field values
    data = {}
    
    reading_peaks = False
    output_peaks = []

    # Open and read the text file
    with open(input_path, 'r') as file:
        for line in file:

            if line.startswith('ACCESSION'):
                # spectrum_id
                data['ACCESSION'] = line.strip().split(' ',1)[1]
                continue
            
            if line.startswith('PK$PEAK:'):
                reading_peaks = True
                match = re.match(r'PK\$PEAK:\sm/z\sint.\srel.int.', line)
                if match is None:
                    raise ValueError(f"Unexpected ordering of PK$PEAK line {line} in file {input_path}")
                continue
            if ':' not in line or '$' not in line:
                # This is either a peak, or a line we don't need so we'll skip
                
                if reading_peaks:
                    # Peaks are in format m/z int. rel.int.
                    if '//' in line:
                        continue
                    peaks = line.strip().split(' ')
                    output_peaks.append((float(peaks[0]), float(peaks[1])))
                    continue
                continue
            
            if reading_peaks:
                reading_peaks = False
            
            # Split each line into key and value using the first ':' as the delimiter
            try:
                key, value = line.strip().split(': ', 1)
            except Exception as e:
                print(f"Failed to split line {line} in file {input_path}")
                raise e
            
            # Check if the key already exists in the dictionary, and if it does, append the value to the existing list
            if key in data:
                data[key].append(value)
            else:
                data[key] = [value]

    # Build double dict for the AC$MASS_SPECTROMETRY 
    if data.get('AC$MASS_SPECTROMETRY') is not None:
        data['AC$MASS_SPECTROMETRY'] = {x.split(' ', 1)[0]: x.split(' ', 1)[1] for x in data.get('AC$MASS_SPECTROMETRY')}
        try:
            if data['AC$MASS_SPECTROMETRY'].get('COLLISION_ENERGY') is not None:
                ce_as_str = data['AC$MASS_SPECTROMETRY'].get('COLLISION_ENERGY') 
                if ce_as_str is not None:
                    match = re.search(r'\d+\.\d+|\d+|\.\d+', ce_as_str)
                    if match is not None:
                        data['AC$MASS_SPECTROMETRY']['COLLISION_ENERGY'] = float(match.group())
        except Exception as e:
            print(f"Failed to parse collision energy for dict {data['AC$MASS_SPECTROMETRY']}")
            raise e
    else:
        # Make it an empty dict so .get() will always return None for us
        data['AC$MASS_SPECTROMETRY'] = {}
        
    if data['AC$MASS_SPECTROMETRY'].get('MS_TYPE') != 'MS2':
        # If it's not MS2, output nothing
        return
    
    # Build double dict for AC$CHROMATOGRAPHY
    if data.get('AC$CHROMATOGRAPHY') is not None:
        data['AC$CHROMATOGRAPHY'] = {x.split(' ', 1)[0]: x.split(' ', 1)[1] for x in data.get('AC$CHROMATOGRAPHY')}
    
    if data.get('MS$FOCUSED_ION') is not None:
        data['MS$FOCUSED_ION'] = {x.split(' ', 1)[0]: x.split(' ', 1)[1] for x in data.get('MS$FOCUSED_ION')}
        if data['MS$FOCUSED_ION'].get('PRECURSOR_M/Z') is not None:
            precursor_mz_as_str = data['MS$FOCUSED_ION'].get('PRECURSOR_M/Z')
            if precursor_mz_as_str is not None:
                match = re.search(r'\d+\.\d+|\d+|\.\d+', precursor_mz_as_str)
                if match is not None:
                    data['MS$FOCUSED_ION']['PRECURSOR_M/Z'] = float(match.group())
    else:
        # Make it an empty dict so .get() will always return None for us
        data['MS$FOCUSED_ION'] = {}
            
    if data.get('CH$EXACT_MASS') is not None:
        data['CH$EXACT_MASS'] = float(data['CH$EXACT_MASS'][0])
    else:
        data['CH$EXACT_MASS'] = None
        
    # Verify SMILES is rdkit parsable:
    if data.get('CH$SMILES', ['']) is not None:
        smiles = data.get('CH$SMILES', [''])[0]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            data['CH$SMILES'] = None
        else:
            data['CH$SMILES'] = smiles
        
    # Build dict of dicts for identifiers
    if data.get('CH$LINK') is not None:
        data['CH$LINK'] = {x.split(' ', 1)[0]: x.split(' ', 1)[1] for x in data.get('CH$LINK')}
        if data['CH$LINK'].get('INCHIKEY') is not None and data.get('CH$SMILES') is not None:
            # Check INCHIkey is valid
            inchikey = data['CH$LINK'].get('INCHIKEY')
            mol = Chem.MolFromSmiles(data['CH$SMILES'])
            generated_key = Chem.inchi.MolToInchiKey(mol)
            if Chem.inchi.MolToInchiKey(mol) != inchikey:
                data['CH$LINK']['INCHIKEY'] = generated_key
    else:
        # Make it an empty dict so .get() will always return None for us
        data['CH$LINK'] = {}
        
    # Check IUPAC is INCHI
    if data.get('CH$IUPAC') is not None:
        if data['CH$IUPAC'][0].startswith('InChI='):
            data['CH$IUPAC'] = data['CH$IUPAC'][0]
        else:
            data['CH$IUPAC'] = None
        # Verify INCHI is rdkit parsable:
        if data.get('CH$IUPAC') is not None:
            mol = Chem.inchi.MolFromInchi(data['CH$IUPAC'])
            if mol is None:
                data['CH$IUPAC'] = None
                
    # Compound name
    if data.get('CH$NAME') is not None:
        data['CH$NAME'] = data['CH$NAME'][0]
        
    # Instrument
    if data.get('AC$INSTRUMENT') is not None:
        data['AC$INSTRUMENT'] = data['AC$INSTRUMENT'][0]
        
    if data.get('AC$INSTRUMENT_TYPE') is not None:
        data['AC$INSTRUMENT_TYPE'] = data['AC$INSTRUMENT_TYPE'][0]
    

    summary_dict = {
        'scan': -1,                                                             # Assign later
        'spectrum_id': data.get('ACCESSION'),
        'collision_energy': data['AC$MASS_SPECTROMETRY'].get('COLLISION_ENERGY'),
        'retention_time': None,                                                 # Here for legacy reasons. Unclear what the unit is and dropped later
        'Adduct': data['MS$FOCUSED_ION'].get('PRECURSOR_TYPE'),
        'Compound_Source': None,
        'Compund_Name': data.get('CH$NAME'),
        'Precursor_MZ': data['MS$FOCUSED_ION'].get('PRECURSOR_M/Z'),
        'ExactMass': data.get('CH$EXACT_MASS'),
        'Charge': None,                                                         # It seems that no charge is stored in the file
        'Ion_Mode': data['AC$MASS_SPECTROMETRY'].get('ION_MODE'),
        'Smiles': data.get('CH$SMILES'),
        'INCHI': data.get('CH$IUPAC'),
        'InChIKey_smiles': data['CH$LINK'].get('INCHIKEY'),
        'InChIKey_inchi': None,
        'msModel': data.get('AC$INSTRUMENT'),
        'msManufacturer': None,                                                 # Will have to be parsed later
        'msDetector': None,                                                     # Will have to be parsed later,  deprecated
        'msMassAnalyzer': None,                                                 # Will have to be parsed later
        'msIonisation': data['AC$MASS_SPECTROMETRY'].get('IONIZATION'),
        'msDissociationMethod': data['AC$MASS_SPECTROMETRY'].get('FRAGMENTATION_MODE'),
        'GNPS_library_membership': "MassBank_ML_Export",
        'GNPS_Inst': data.get('AC$INSTRUMENT_TYPE')                  # It's not from gnps, but we'll keep the column name the same and parse it later
    }
    
    if output_path is None:
        output_path = f'output_{summary_dict["spectrum_id"]}'
    if output_path is not None:
        dir_path = os.path.dirname(output_path)
        if dir_path != '':
            if not os.path.isdir(dir_path):
                os.makedirs(dir_path)
        output_path = os.path.join(dir_path, os.path.basename(output_path))
    
    # Write csv file
    with open(output_path + '.csv', 'w') as output_csv:
        w = csv.DictWriter(output_csv, ['scan', 'spectrum_id','collision_energy','retention_time','Adduct','Compound_Source','Compund_Name','Precursor_MZ','ExactMass','Charge','Ion_Mode','Smiles','INCHI','InChIKey_smiles','InChIKey_inchi','msModel','msManufacturer','msDetector','msMassAnalyzer','msIonisation','msDissociationMethod','GNPS_library_membership','GNPS_Inst'])
        w.writerow(summary_dict)
    # Write mgf file
    with open(output_path + '.mgf', 'w') as output_mgf:
        output_mgf.write("BEGIN IONS\n")
        output_mgf.write("PEPMASS={}\n".format(summary_dict["Precursor_MZ"]))
        output_mgf.write("CHARGE={}\n".format(summary_dict["Charge"]))
        output_mgf.write("MSLEVEL={}\n".format(2))
        output_mgf.write("TITLE="+str(summary_dict.get("spectrum_id")+"\n"))
        output_mgf.write("SCANS={}\n".format(summary_dict['scan']))
        for peak in output_peaks:
            output_mgf.write("{} {}\n".format(peak[0], peak[1]))

        output_mgf.write("END IONS\n")

def main():
    parser = argparse.ArgumentParser(description='Process MassBank file')
    parser.add_argument('--input', help='Input MassBank file')
    parser.add_argument('--blacklist', help='Blacklist file')
    parser.add_argument('--output', help='Output file directory and basename', default=None)
    args = parser.parse_args()
    
    # Open black list and see if it's a blacklisted file
    with open(args.blacklist, 'r') as f:
        blacklist = f.read().splitlines()
        # Check if the filename is a substring of any of the blacklisted files
        if any(args.input in x.split('/', 1) for x in blacklist):
            return
    
    process_file(args.input, output_path=args.output)

if __name__ == "__main__":
    main()


def test_data_processing():
    # Tests can be run with  python3 -m pytest ./bin/MassBank_processing.py 
    import pandas as pd
    path = './data/tests/MSBNK-AAFC-AC000870.txt'
    process_file(path, output_path = './test_outputs/test.csv')
    
    # Compare csvs
    df = pd.read_csv('./test_outputs/test.csv')
    label_df = pd.read_csv('./data/tests/MSBNK-AAFC-AC000870_output.csv')
    
    # Compare mgfs
    with open('./test_outputs/test.mgf', 'r') as f:
        mgf = f.read()
        with open('./data/tests/MSBNK-AAFC-AC000870_output.mgf', 'r') as f:
            label_mgf = f.read()
            assert mgf == label_mgf
    assert df.equals(label_df)