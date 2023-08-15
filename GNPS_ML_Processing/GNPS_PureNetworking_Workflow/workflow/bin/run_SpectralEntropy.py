import sys
import os
import numpy as np
import argparse
import pandas as pd
from tqdm import tqdm
from spectral_entropy.spectral_similarity import multiple_similarity

#path_to_remove = "/" + os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))).replace('\\', '/') + "/GNPS_sharedcode"
#sys.path.insert(0, path_to_remove)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.realpath(__file__)), "./GNPS_sharedcode"))
import ming_spectrum_library


def calculate_entropy(spectrum1_dict, spectrum2_dict, alignment_params={}):
    mzi1 = np.array(spectrum1_dict['peaks']).T.astype(float)
    mzi2 = np.array(spectrum2_dict['peaks']).T.astype(float)

    diff = np.abs(mzi1[0][:, None] - mzi2[0][None, :])
    mask = diff <= float(alignment_params.get("ms2_tolerance", 0.1))
    matched_peaks = np.count_nonzero(mask.any(axis=1))

    scores = {}
    scores["entropy_score"] = float(0)
    scores["unweighted_entropy_score"] = float(0)
    scores["matched_peaks"] = int(0)

    if alignment_params.get("min_matched_peaks", 6) <= matched_peaks:
        similarity_entropy = multiple_similarity(spectrum_query=mzi1,
                                                 spectrum_library=mzi2,
                                                 methods=['entropy', 'unweighted_entropy'],
                                                 need_normalize_result=True,
                                                 ms2_da=float(alignment_params.get("peak_tolerance", 0.1)))


        scores["entropy_score"] = float(similarity_entropy['entropy'])
        scores["unweighted_entropy_score"] = float(similarity_entropy['unweighted_entropy'])
        scores["matched_peaks"] = int(matched_peaks)

    return scores


def divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i:i + n]


def main():
    parser = argparse.ArgumentParser(description='Reformat entropy')
    parser.add_argument('input_mgf', help='input mgf')
    parser.add_argument('output_results', help='output results tsv')
    parser.add_argument('--nodenumber', type=int, default=0, help='current node number')
    parser.add_argument('--nodetotal', type=int, default=1, help='total node number')
    parser.add_argument('--min_cosine', default=0.5, type=float, help='min_cosine')
    parser.add_argument('--min_matched_peaks', default=0.5, type=float, help='min_matched_peaks')
    parser.add_argument('--ms2_tolerance', default=0.5, type=float, help='ms2_tolerance')

    args = parser.parse_args()

    spectrum_collection = ming_spectrum_library.SpectrumCollection(args.input_mgf)
    spectrum_collection.load_from_file()


    output_results_list = []

    nodetotal = args.nodetotal
    nodenumber = args.nodenumber

    if len(spectrum_collection.spectrum_list) < nodetotal:
        nodetotal = 1

    chunk_size = int(len(spectrum_collection.spectrum_list) / nodetotal)
    if nodenumber == nodetotal - 1:
        sliced_spectrum_list = spectrum_collection.spectrum_list[nodenumber * chunk_size:]
    else:
        sliced_spectrum_list = spectrum_collection.spectrum_list[
                               nodenumber * chunk_size:nodenumber * chunk_size + chunk_size]

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

                    ms2_tolerance = args.ms2_tolerance
                    min_matched_peaks = args.min_matched_peaks

                    scores = calculate_entropy(spectrum1_dict, spectrum2_dict, alignment_params={'ms2_tolerance': ms2_tolerance,
                                                                                                 'min_matched_peaks': min_matched_peaks})

                    min_cosine = args.min_cosine


                    if (scores['matched_peaks'] < min_matched_peaks) or (scores['entropy_score'] <= min_cosine):
                        continue
                except KeyboardInterrupt:
                    raise
                except Exception as e:
                    print(f"Error: {str(e)}")
                    continue

                #print(scores, spectrum1.scan, spectrum2.scan)


                result_dict = {}
                result_dict["scan1"] = spectrum1.scan
                result_dict["scan2"] = spectrum2.scan
                result_dict["entropy_score"] = scores['entropy_score']
                result_dict["unweighted_entropy_score"] = scores['unweighted_entropy_score']
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

