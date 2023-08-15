import sys
import os
import simile.simile as sml
import numpy as np
import argparse
import pandas as pd
from tqdm import tqdm

sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "./GNPS_sharedcode"))
import ming_spectrum_library

def calculate_simile(spectrum1_dict, spectrum2_dict, alignment_params={}):
    s1 = spectrum1_dict
    s2 = spectrum2_dict

    # import json
    # s1 = json.loads(s1)
    mz1 = [float(x[0]) for x in s1['peaks']]
    i1 = [float(x[1]) for x in s1['peaks']]
    mzi1= np.asarray([mz1,i1])

    # s2 = json.loads(s2)
    mz2 = [float(x[0]) for x in s2['peaks']]
    i2 = [float(x[1]) for x in s2['peaks']]
    mzi2 = np.asarray([mz2,i2])

    # Generate pair-specific substitution matrix
    # S = sml.similarity_matrix(s1[0], s2[0], tolerance=0.01)
    S = sml.similarity_matrix(mzi1[0], mzi2[0], tolerance=float(alignment_params.get("peak_tolerance", 0.1)))

    # Align and score using upper-right quadrant of substitution matrix
    simile_score, simile_alignment1 = sml.pairwise_align(S[:mzi1.shape[1], mzi1.shape[1]:])

    # Calculate significance of the alignment
    pval = sml.significance_test(S, mzi1[0], mzi2[0], max_log_iter=2)

    scores = {}
    scores["score"] = simile_score
    scores["pval"] = pval
    scores["matched_peaks"] = len(simile_alignment1)
    
    return scores

def divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n): 
        yield l[i:i + n]

def main():
    parser = argparse.ArgumentParser(description='Reformat blink')
    parser.add_argument('input_mgf', help='input mgf')
    parser.add_argument('output_results', help='output results tsv')
    parser.add_argument('--nodenumber', type=int, default=0, help='current node number')
    parser.add_argument('--nodetotal', type=int, default=1, help='total node number')
    
    args = parser.parse_args()

    spectrum_collection = ming_spectrum_library.SpectrumCollection(args.input_mgf)
    spectrum_collection.load_from_file()

    alignment_params = {}
    alignment_params["peak_tolerance"] = 0.01

    output_results_list = []

    nodetotal = args.nodetotal
    nodenumber = args.nodenumber

    if len(spectrum_collection.spectrum_list) < nodetotal:
        nodetotal = 1

    chunk_size = int(len(spectrum_collection.spectrum_list) / nodetotal)
    if nodenumber == nodetotal - 1:
        sliced_spectrum_list = spectrum_collection.spectrum_list[nodenumber * chunk_size:]
    else:
        sliced_spectrum_list = spectrum_collection.spectrum_list[nodenumber * chunk_size:nodenumber * chunk_size + chunk_size]

    for spectrum1 in tqdm(sliced_spectrum_list):
        if len(spectrum1.peaks) < 2:
            continue

        for spectrum2 in spectrum_collection.spectrum_list:
            if len(spectrum2.peaks) < 2:
                continue

            if spectrum1.scan > spectrum2.scan:
                spectrum1_dict = {}
                spectrum2_dict = {}

                spectrum1_dict["peaks"] = spectrum1.peaks
                spectrum2_dict["peaks"] = spectrum2.peaks

                try:
                    scores = calculate_simile(spectrum1_dict, spectrum2_dict, alignment_params=alignment_params)
                except KeyboardInterrupt:
                    raise
                except:
                    print("Error")
                    continue

                print(scores, spectrum1.scan, spectrum2.scan)

                result_dict = {}
                result_dict["scan1"] = spectrum1.scan
                result_dict["scan2"] = spectrum2.scan
                result_dict["score"] = scores['score']
                result_dict["pval"] = scores['pval']
                result_dict["matched_peaks"] = scores['matched_peaks']

                output_results_list.append(result_dict)

    df = pd.DataFrame(output_results_list)

    print(df)

    if len(df) > 0:
        # Renaming columns
        df["CLUSTERID1"] = df["scan1"]
        df["CLUSTERID2"] = df["scan2"]

        df.to_csv(args.output_results, sep="\t", index=False)

if __name__ == "__main__":
    main()
