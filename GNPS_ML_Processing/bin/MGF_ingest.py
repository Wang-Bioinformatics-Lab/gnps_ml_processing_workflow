import argparse
import csv
from pathlib import Path

from tqdm import tqdm
from matchms.importing import load_from_mgf

# Import headers from main processor
from GNPS2_Processor import csv_headers

from pyteomics.mgf import MGF as py_mgf

def parse_mgf(mgf_path:Path, summary_output_path:Path, spectra_output_path:Path)->None:
    """
    Parse an MGF file and generate summary and spectra output files.

    Args:
        mgf_path (Path): The path to the MGF file to be parsed.
        summary_output_path (Path): The path to the summary output file.
        spectra_output_path (Path): The path to the spectra output file.

    Raises:
        ValueError: If the script does not support alternative ms levels. Only MS/MS spectra should be supplied.

    Returns:
        None
    """   
    summary_output_path.parent.mkdir(parents=True, exist_ok=True)
    spectra_output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with py_mgf(str(mgf_path)) as mgf_file:
        with open(summary_output_path, 'w', encoding='utf-8') as output_summary:
            output_summary_writer = csv.DictWriter(output_summary, csv_headers)
            with open(spectra_output_path, 'w', encoding='utf-8') as output_mgf:
                
                for idx, spectrum in enumerate(tqdm(mgf_file)):
                    summary_dict = {}
                    
                    if spectrum['params'].get('ms_level') is not None:
                        if int(spectrum.metadata.get('ms_level')) != 2:
                            raise ValueError("This script does not support alternative ms levels. Only MS/MS spectra should be supplied.")
                        
                    summary_dict['scan']                    = idx+1
                    summary_dict['spectrum_id']             = spectrum['params'].get('spectrumid')
                    if summary_dict['spectrum_id'] is None:
                        summary_dict['spectrum_id']         = spectrum['params'].get('spectrum_id')
                    if summary_dict['spectrum_id'] is None:
                        summary_dict['spectrum_id']         = f"{mgf_path.stem}_scan_{idx+1}"
                    summary_dict['collision_energy']        = spectrum['params'].get('collision_energy')
                    summary_dict['retention_time']          = spectrum['params'].get('retention_time')
                    summary_dict['Adduct']                  = spectrum['params'].get('adduct')
                    if summary_dict['Adduct'] is None:
                        # Fall back to precursortype
                        summary_dict['Adduct'] = spectrum['params'].get('precursortype')
                    if summary_dict['Adduct'] is None:
                        print(f"Warning: Failed to get adduct for scan {idx+1} in {mgf_path.stem}")
                        print(f"Metadata: {spectrum['params']}")
                    summary_dict['Compound_Source']         = 'MGF Import'
                    summary_dict['Compound_Name']           = spectrum['params'].get('compound_name')
                    if summary_dict['Compound_Name'] is None:
                        summary_dict['Compound_Name']       = spectrum['params'].get('name')
                    summary_dict['Precursor_MZ']            = spectrum['params'].get('pepmass')
                    if summary_dict['Precursor_MZ'] is None:
                        summary_dict['Precursor_MZ']            = spectrum['params'].get('precursor_mz')
                    if isinstance(summary_dict['Precursor_MZ'], (list, tuple)):
                        summary_dict['Precursor_MZ'] = summary_dict['Precursor_MZ'][0]
                    summary_dict['ExactMass']               = None # Postprocessing will get this
                    summary_dict['Charge']                  = None   # Postprocessing will get this
                    summary_dict['Ion_Mode']                = spectrum['params'].get('ionmode')
                    if summary_dict['Ion_Mode'] is None:
                        summary_dict['Ion_Mode']            = spectrum['params'].get('ion_mode')
                    summary_dict['Smiles']                  = spectrum['params'].get('smiles')
                    summary_dict['INCHI']                   = spectrum['params'].get('inchikey')
                    summary_dict['InChIKey_smiles']         = spectrum['params'].get('inchikey')
                    summary_dict['InChIKey_inchi']          = None
                    summary_dict['msModel']                 = spectrum['params'].get('instrument')   # This seems to contain more specific information (incl column)
                    summary_dict['msManufacturer']          = None  # Defer parsing these from msModel, GNPS_Inst until postprocessing
                    summary_dict['msDetector']              = None
                    summary_dict['msMassAnalyzer']          = None
                    summary_dict['msIonisation']            = spectrum['params'].get('ion_source')
                    if summary_dict['msIonisation'] is None:
                        summary_dict['msIonisation']        = spectrum['params'].get('ionization')
                    if summary_dict['msIonisation'] is None:
                        summary_dict['msIonisation']        = spectrum['params'].get('ionisation')
                    if summary_dict['msIonisation'] is None:
                        summary_dict['msIonisation']        = spectrum['params'].get('ionsource')
                    summary_dict['msDissociationMethod']    = None
                    summary_dict['GNPS_library_membership'] = mgf_path.stem
                    summary_dict['GNPS_Inst']               = spectrum['params'].get('instrument_type')
                    if summary_dict['GNPS_Inst'] is None:
                        summary_dict['GNPS_Inst']           = spectrum['params'].get('instrument')
                    
                    output_summary_writer.writerow(summary_dict)
                    
                    # Not implemented elsewhere in the pipeline
                    if False:
                        ccs = spectrum['params'].get('ccs')
                        summary_dict['ccs']                     = None if ccs is None or ccs == -1 else ccs
                        summary_dict['compound_class']          = spectrum['params'].get('compound_class')
                
                    output_mgf.write("BEGIN IONS\n")
                    output_mgf.write(f"PEPMASS={summary_dict['Precursor_MZ']}\n")
                    # output_mgf.write(f"CHARGE={summary_dict['Charge']}\n")
                    # This is a bit dicey, we raise the error above if the ms_level is known but we will rely on the fact that only ms2 data is supplied
                    output_mgf.write(f"MSLEVEL={2}\n")  
                    output_mgf.write(f"TITLE={summary_dict['spectrum_id']}\n")
                    output_mgf.write(f"SCANS={idx+1}\n")

                    for mz, i in zip(spectrum['m/z array'], spectrum['intensity array']):
                        output_mgf.write(f"{mz} {i}\n")

                    output_mgf.write("END IONS\n")

def main():
    parser = argparse.ArgumentParser(description='MGF Ingest')
    parser.add_argument("--mgf_path", help="Path to the MGF file")
    parser.add_argument("--summary_output_path", help="Path to output the summary")
    parser.add_argument("--spectra_output_path", help="Path to output the spectra")
    args = parser.parse_args()
    
    mgf_path            = Path(args.mgf_path)
    summary_output_path = Path(args.summary_output_path)
    spectra_output_path = Path(args.spectra_output_path)
    
    parse_mgf(mgf_path, summary_output_path, spectra_output_path)

if __name__ == "__main__":
    main()