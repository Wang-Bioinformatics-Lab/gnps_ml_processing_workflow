import argparse
import csv
from pathlib import Path

from tqdm import tqdm
from matchms.importing import load_from_msp

# Import headers from main processor
from GNPS2_Processor import csv_headers

def parse_msp(msp_path:Path, summary_output_path:Path, spectra_output_path:Path)->None:
    """
    Parse an MSP file and generate summary and spectra output files.

    Args:
        msp_path (Path): The path to the MSP file to be parsed.
        summary_output_path (Path): The path to the summary output file.
        spectra_output_path (Path): The path to the spectra output file.

    Raises:
        ValueError: If the script does not support alternative ms levels. Only MS/MS spectra should be supplied.

    Returns:
        None
    """
    spectra_generator = load_from_msp(msp_path)
    
    summary_output_path.parent.mkdir(parents=True, exist_ok=True)
    spectra_output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(summary_output_path, 'w', encoding='utf-8') as output_summary:
        output_summary_writer = csv.DictWriter(output_summary, csv_headers)
        
        with open(spectra_output_path, 'w', encoding='utf-8') as output_mgf:
            
            for idx, spectrum in enumerate(tqdm(spectra_generator)):
                summary_dict = {}
                
                if spectrum.metadata.get('ms_level') is not None:
                    if int(spectrum.metadata.get('ms_level')) != 2:
                        raise ValueError("This script does not support alternative ms levels. Only MS/MS spectra should be supplied.")
                    
                summary_dict['scan']                    = idx+1
                summary_dict['spectrum_id']             = f"{msp_path.stem}_scan_{idx+1}"
                summary_dict['collision_energy']        = spectrum.metadata.get('collision_energy')
                summary_dict['retention_time']          = spectrum.metadata.get('retention_time')
                summary_dict['Adduct']                  = spectrum.metadata.get('adduct')
                if summary_dict['Adduct'] is None:
                    # Fall back to precursortype
                    summary_dict['Adduct'] = spectrum.metadata.get('precursortype')
                if summary_dict['Adduct'] is None:
                    print(f"Warning: Failed to get adduct for scan {idx+1} in {msp_path.stem}")
                    print(f"Metadata: {spectrum.metadata}")
                summary_dict['Compound_Source']         = 'MSP Import'
                summary_dict['Compound_Name']           = spectrum.metadata.get('compound_name')
                summary_dict['Precursor_MZ']            = spectrum.metadata.get('precursor_mz')
                summary_dict['ExactMass']               = None # Postprocessing will get this
                summary_dict['Charge']                  = None   # Postprocessing will get this
                summary_dict['Ion_Mode']                = spectrum.metadata.get('ionmode')
                summary_dict['Smiles']                  = spectrum.metadata.get('smiles')
                summary_dict['INCHI']                   = spectrum.metadata.get('inchikey')
                summary_dict['InChIKey_smiles']         = spectrum.metadata.get('inchikey')
                summary_dict['InChIKey_inchi']          = None
                summary_dict['msModel']                 = spectrum.metadata.get('instrument')   # This seems to contain more specific information (incl column)
                summary_dict['msManufacturer']          = None  # Defer parsing these from msModel, GNPS_Inst until postprocessing
                summary_dict['msDetector']              = None
                summary_dict['msMassAnalyzer']          = None
                summary_dict['msIonisation']            = None
                summary_dict['msDissociationMethod']    = None
                summary_dict['GNPS_library_membership'] = msp_path.stem
                summary_dict['GNPS_Inst']               = spectrum.metadata.get('instrument_type')
                
                output_summary_writer.writerow(summary_dict)
                
                # Not implemented elsewhere in the pipeline
                if False:
                    ccs = spectrum.metadata.get('ccs')
                    summary_dict['ccs']                     = None if ccs is None or ccs == -1 else ccs
                    summary_dict['compound_class']          = spectrum.metadata.get('compound_class')
            
                output_mgf.write("BEGIN IONS\n")
                output_mgf.write(f"PEPMASS={summary_dict['Precursor_MZ']}\n")
                # output_mgf.write(f"CHARGE={summary_dict['Charge']}\n")
                # This is a bit dicey, we raise the error above if the ms_level is known but we will rely on the fact that only ms2 data is supplied
                output_mgf.write(f"MSLEVEL={2}\n")  
                output_mgf.write(f"TITLE={summary_dict['spectrum_id']}\n")
                output_mgf.write(f"SCANS={idx+1}\n")

                peaks = list(spectrum.peaks)
                for peak in peaks:
                    output_mgf.write("{} {}\n".format(peak[0], peak[1]))

                output_mgf.write("END IONS\n")

def main():
    parser = argparse.ArgumentParser(description='MSP Ingest')
    parser.add_argument("--msp_path", help="Path to the MSP file")
    parser.add_argument("--summary_output_path", help="Path to output the summary")
    parser.add_argument("--spectra_output_path", help="Path to output the spectra")
    args = parser.parse_args()
    
    msp_path            = Path(args.msp_path)
    summary_output_path = Path(args.summary_output_path)
    spectra_output_path = Path(args.spectra_output_path)
    
    parse_msp(msp_path, summary_output_path, spectra_output_path)

if __name__ == "__main__":
    main()