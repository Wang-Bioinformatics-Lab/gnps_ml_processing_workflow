import argparse
import pathlib
import pyteomics.mgf
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed  
from glob import glob
import vaex
import re
from utils import synchronize_spectra

def generate_mgf(parquet_file, output_path):
    name = parquet_file[:-8]
    # df = pd.read_parquet(parquet_file)
    df = vaex.open(parquet_file)
    # to_mgf = []
    spectrum_ids = df.spectrum_id.unique()
    # for i, spectrum_id in enumerate(spectrum_ids):
    #     # spectra = df.loc[df.spectrum_id == spectrum_id]
    #     spectra = df[df.spectrum_id == spectrum_id]
    #     to_mgf.append({'m/z array':spectra.mz.values, 'intensity array':spectra.i.values, 'params':{'SCAN': i, 'TITLE':spectrum_id,'prec_mz':spectra.prec_mz.values[0][0]}})
    # pyteomics.mgf.write(to_mgf, output=name+'.mgf',write_charges=False, write_ions=False)
    
    
    # output_mgf = open(name+'.mgf', "w")
    output_mgf = open(output_path, "w")
    for i, spectrum_id in enumerate(spectrum_ids):
        spectra = df[df.spectrum_id == spectrum_id]
        
    
        output_mgf.write("BEGIN IONS\n")
        output_mgf.write("PEPMASS={}\n".format(spectra.prec_mz.values[0][0]))
        output_mgf.write("TITLE={}\n".format(spectrum_id))
        output_mgf.write("SCANS={}\n".format(i))

        peaks = zip(spectra.mz.values, spectra.i.values)
        for peak in peaks:
            output_mgf.write("{} {}\n".format(peak[0], peak[1]))

        output_mgf.write("END IONS\n")
        
    output_mgf.close()
    
def main():
    parser = argparse.ArgumentParser(description='Generate MGF Files.')
    # parser.add_argument('-p', type=int, help='number or processors to use', default=10)
    # parser.add_argument('-path', type=str, help='Path to locate parquet files', default='./*.parquet')
    parser.add_argument('-input_path', type=str, help='Path to locate parquet files')
    # parser.add_argument('-output_path', type=str, help='Path to output MGF files')
    args = parser.parse_args()
    
    # files = glob(args.path)
    # parquet_outputs = [x for x in files if x != './ALL_GNPS_cleaned.parquet']   # We don't want to generate an MGF for the entire dataset
    
    # Parallel(n_jobs=args.p)(delayed(generate_mgf)(parquet_file) for parquet_file in tqdm(parquet_outputs))
    
    # Little hack to get nextflow to preserve the directories for the outputs
    # slash_indices = [x.start() for x in re.finditer(r"/",args.input_path)]
    # if args.input_path != './ALL_GNPS_cleaned.parquet':
    #     generate_mgf(args.input_path, args.input_path[slash_indices[2]::-7]+"mgf")

    # Get extension of args.input_path
    ext = pathlib.Path(args.input_path).suffix
    if ext == ".parquet":
        if 'ALL_GNPS_cleaned.parquet' not in args.input_path:
            generate_mgf(args.input_path, args.input_path.split('/')[-1].split('.')[0]+".mgf")
    else:
        raise IOError(f"Expected a parquet file bug got {ext}")
        
if __name__ == '__main__':
    main()