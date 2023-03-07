import argparse
import pyteomics.mgf
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed  
from glob import glob

def generate_mgf(parquet_file):
    name = parquet_file[:-8]
    df = pd.read_parquet(parquet_file)
    to_mgf = []
    spectrum_ids = df.spectrum_id.unique()
    for i, spectrum_id in enumerate(spectrum_ids):
        spectra = df.loc[df.spectrum_id == spectrum_id]
        to_mgf.append({'m/z array':list(spectra.mz), 'intensity array':list(spectra.i), 'params':{'SCAN': i, 'TITLE':spectrum_id,'prec_mz':spectra.prec_mz.iloc[0][0]}})
    pyteomics.mgf.write(to_mgf, output=name+'.mgf',write_charges=False, write_ions=False)
    
def main():
    parser = argparse.ArgumentParser(description='Generate MGF Files.')
    parser.add_argument('-p', type=int, help='number or processors to use', default=10)
    args = parser.parse_args()
    
    files = glob('./*.parquet')
    parquet_outputs = [x for x in files if x != './ALL_GNPS_cleaned.parquet']
    
    Parallel(n_jobs=args.p)(delayed(generate_mgf)(parquet_file) for parquet_file in tqdm(parquet_outputs))
        
if __name__ == '__main__':
    main()