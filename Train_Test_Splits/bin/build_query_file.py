import os
import uuid
import json
import pandas as pd
from pathlib import Path
import subprocess
import tempfile
import argparse
import ast
from pyteomics.mgf import IndexedMGF

# Debugging
import time

# Helping to keep caching things in memory, specifically for library correspondence
from expiringdict import ExpiringDict
memory_cache = ExpiringDict(max_len=100, max_age_seconds=84600)

METADATA_PREFIX = 'Index Metadata:'
MASSIVE_ID = 'Precursor'
GNPS_ID = 'GNPSLibraryAccession'
NO_ID = 'Annotation'

def search_library(input_mgf,
                    library_name,
                    analog_search=False,
                    lower_delta=130,
                    upper_delta=200,
                    pm_tolerance=0.2,
                    fragment_tolerance=0.2,
                    cosine_threshold=0.7):

    TIMEOUT_SEC = 600 # Originally 180
    TIMEOUT_MESSAGE = "Query timed out in {} seconds.".format(TIMEOUT_SEC)

    starttime = time.time()

    is_index = 'index' in library_name
    
    # temp_dir = Path(tempfile.gettempdir())
    temp_dir = Path('./')

    random_string = str(uuid.uuid4()).replace("-", "")
    temp_query_mgf = temp_dir.joinpath("query_{}.mgf".format(random_string))
    temp_results_tsv = temp_dir.joinpath("{}.tsv".format(random_string))
    temp_results_json_folder = temp_dir.joinpath("{}".format(random_string))
    temp_results_json_folder.mkdir()

    results_summary_json = temp_dir.joinpath("{}_summary.json".format(random_string))

    spectra = IndexedMGF(input_mgf, index_by_scans=True)

    with open(temp_query_mgf, "w") as o:
        for _, spectrum in enumerate(spectra):
            scan = spectrum['params']['scans']
            charge = int(spectrum['params']['charge'][0])
            
            if charge == 0:
                charge = 1
                
            o.write("BEGIN IONS\n")
            o.write("SEQ=*..*\n")
            o.write(f"SCANS={scan}\n")
            o.write(f"CHARGE={charge}\n")
            o.write("MSLEVEL=2\n")
            o.write(f"PEPMASS={spectrum['params']['pepmass'][0]}\n")
            for peak in zip(spectrum['m/z array'], spectrum['intensity array']):
                o.write(f"{peak[0]} {peak[1]}\n")
            o.write("END IONS\n")
    
    # Return to command line
    print(temp_query_mgf, cosine_threshold)
    return

    library_index_path = Path("./libraries").joinpath(library_name).joinpath("build_index")
    # library_spectrum_correspondence_path = os.path.join("./libraries", library_name, "library_id_correspondence.json")
    # library_spectrum_correspondence_feather = os.path.join("./libraries", library_name, "library_id_correspondence.feather")

    cmd = [
        Path('./bin/main_execmodule'), 'ExecIndex', Path('./bin/generic_params'),
        '-autoINPUT_INDEX', library_index_path,
        '-autoOUTPUT_RESULTS', temp_results_tsv,
        '-autoINPUT_SPECS', temp_query_mgf,
        '-autoPM_TOLERANCE', str(float(pm_tolerance)),
        '-autoFRAG_TOLERANCE', str(float(fragment_tolerance)),
        '-autoDELTA_MZ_ABOVE', str(int(lower_delta) if analog_search else 0),
        '-autoDELTA_MZ_BELOW', str(int(upper_delta) if analog_search else 0),
        '-autoTHETA', str(float(cosine_threshold)),
        '-autoSPECTRUM_DISK_ACCESS_POLICY', 'DISK', 
        '-autoINDEX_DISK_ACCESS_POLICY', 'DISK',
        '-autoOUTPUT_FILTERED_QUERIES_JSON_FOLDER', temp_results_json_folder,
        '-autoDELTA_F', '0.4' if is_index else '0.2',
        '-autoVALIDATE', '0'
    ]

    # The following commands make things fast but also not fully working. 
    #-autoSPECTRUM_DISK_ACCESS_POLICY DISK \
    #-autoINDEX_DISK_ACCESS_POLICY DISK \

    # timing_list = []
    
    print("1", time.time() - starttime)
    # timing_list.append({"desc" : "prep", "timing": time.time() - starttime})
    starttime = time.time()
    
    import sys
    cmd_string = ' '.join([str(s) for s in cmd])
    print(cmd_string, file=sys.stderr, flush=True)

    log_txt = cmd_string + "\n---\n"

    try:
        proc = subprocess.run([s for s in cmd], 
                                text=True,
                                timeout=TIMEOUT_SEC,
                                capture_output=True
                            )
        log_txt += proc.stderr
        log_txt += proc.stdout
    except subprocess.TimeoutExpired as e:
        return pd.DataFrame(), [], pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), TIMEOUT_MESSAGE + "\n---\n" + log_txt + e.stderr + e.stdout
    
    try:
        print("About to read output")
        results_df = pd.read_csv(temp_results_tsv, sep="\t")
    except:
        return pd.DataFrame(), [], pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), log_txt

    if len(results_df) == 0:
        return results_df, [], pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), log_txt


    print("*************", temp_results_tsv)
    # Let's remove the results from the file system
    # try:
    #     os.remove(temp_results_tsv)
    # except:
    #     pass

    # Writing small summary for tracking purposes
    summary_dict = {}
    summary_dict["number_results"] = len(results_df)

    # write json
    with open(results_summary_json, 'w') as outfile:
        json.dump(summary_dict, outfile)

    # Processing the results

    # print("X2", time.time() - starttime)
    # timing_list.append({"desc" : "Search", "timing": time.time() - starttime})
    # starttime = time.time()

    # results_df = results_df.drop("Query Filename", axis=1)
    # results_df = results_df.drop("Query Scan", axis=1)

    # print("X4", time.time() - starttime)
    # timing_list.append({"desc" : "Parse Query", "timing": time.time() - starttime})
    # starttime = time.time()

    results_df["libraryname"] = results_df["Index Filename"].apply(lambda x: os.path.basename(x))
    results_df["uniquekey"] = results_df["libraryname"] + ":" + results_df["Index Scan"].astype(str)

    # print("X5", time.time() - starttime)
    # timing_list.append({"desc" : "Parse Query 2", "timing": time.time() - starttime})
    # starttime = time.time()

    # print("X6", time.time() - starttime)
    # timing_list.append({"desc" : "Get DB Info", "timing": time.time() - starttime})
    # starttime = time.time()
    
        # This is raw data
        # beware!! (TODO?) assumes MSV identifier
    merged_df = results_df.round(2)

    all_metadata = set([k.replace(METADATA_PREFIX,'') for k in merged_df if METADATA_PREFIX in k])
    
    identification_column = None
    
    for potential_identification in [MASSIVE_ID,GNPS_ID]:
        if potential_identification in all_metadata:
            identification_column = potential_identification
            all_metadata.remove(identification_column)
            break

    merged_df = merged_df.rename(lambda x: x.replace(METADATA_PREFIX,''), axis='columns')

    merged_df["Status"] = merged_df[identification_column].apply(lambda x: 'NoID' if pd.isna(x) else 'HasID') if identification_column else 'NoID'

    if GNPS_ID in merged_df:
        merged_df["Dataset"] = 'GNPS Library'
        merged_df["USI"] = merged_df[GNPS_ID].apply(lambda x: "mzspec:GNPS:GNPS-LIBRARY:accession:{}".format(x))
    else:
        merged_df["Dataset"] = merged_df["Index Filename"].apply(lambda x: x[:12] if 'MSV' in x else x.split('/')[1])
        merged_df["USI"] = merged_df["Index Filename"].apply(lambda x: "mzspec:{}:{}:scan:".format(*((x[:12], Path(x[13:]).stem) if 'MSV' in x else ('MassIVE','TASK-{}/{}'.format('-'.join(x.split('/')[1:-1]),x.split('/')[-1])))))
        merged_df["USI"] = merged_df["USI"] + merged_df["Index Scan"].astype(str)

    merged_df["Charge"] = merged_df["Index Charge"]
    merged_df["Filtered Input Spectrum Path"] = str(next(temp_results_json_folder.glob('*')))
    merged_df["Unit Delta Mass"] = merged_df["Delta Mass"].apply(lambda x: int(round(x)))
    merged_df = merged_df[["Delta Mass"] + ([identification_column] if identification_column else []) + ["USI", "Charge", "Cosine", "Matching Peaks", "Unit Delta Mass", "Dataset", "Status"] + list(all_metadata) + ["Query Filename", "Query Scan", "Index UnitPM", "Index IdxInUnitPM", "Filtered Input Spectrum Path"]]
    merged_df = merged_df.drop_duplicates()
    merged_df = merged_df.sort_values("Cosine", ascending=False)

    print("X7", time.time() - starttime)
    timing_list.append({"desc" : "Demangle", "timing": time.time() - starttime})
    starttime = time.time()

    grouped_by_unit_delta_mass = merged_df.groupby("Unit Delta Mass").size().reset_index(name='Frequency')
    grouped_by_dataset = merged_df.groupby("Dataset").size().reset_index(name='Frequency')
    grouped_by_dataset_unit_delta = merged_df.groupby(["Dataset","Unit Delta Mass"]).size().reset_index(name='Frequency')
    
    grouped_by_unit_delta_mass = grouped_by_unit_delta_mass[ ['Frequency'] + [ col for col in grouped_by_unit_delta_mass.columns if col != 'Frequency' ] ]
    grouped_by_dataset = grouped_by_dataset[ ['Frequency'] + [ col for col in grouped_by_dataset.columns if col != 'Frequency' ] ]
    grouped_by_dataset_unit_delta = grouped_by_dataset_unit_delta[ ['Frequency'] + [ col for col in grouped_by_dataset_unit_delta.columns if col != 'Frequency' ] ]

    grouped_by_unit_delta_mass = grouped_by_unit_delta_mass.sort_values("Frequency", ascending=False)
    grouped_by_dataset = grouped_by_dataset.sort_values("Frequency", ascending=False)
    grouped_by_dataset_unit_delta = grouped_by_dataset_unit_delta.sort_values("Frequency", ascending=False)

    metadata_groupings = []

    for metadata in [identification_column] + list(all_metadata):
    
        just_ids = pd.DataFrame(columns = merged_df.columns)
        if metadata:
            just_ids = merged_df[merged_df[metadata] != '']

        if len(just_ids.index) > 0:
            frequency_per_id = just_ids.groupby([metadata]).size().reset_index(name='Frequency')
            grouped_by_id = just_ids.groupby([metadata]).first().reset_index()
            grouped_by_id = grouped_by_id[[metadata, "Cosine", "Delta Mass", "USI", "Charge"]]
            grouped_by_id = grouped_by_id.rename(columns=lambda s: "Best {}".format(s) if s != metadata else s)
            grouped_by_id.insert(0,"Frequency",frequency_per_id["Frequency"])
            grouped_by_id = grouped_by_id.sort_values("Frequency", ascending=False)
            metadata_groupings.append((metadata,grouped_by_id))


    print("X8", time.time() - starttime)
    timing_list.append({"desc" : "Group", "timing": time.time() - starttime})

    # Enriching with Dataset Information
    try:
        grouped_by_dataset = _encrich_dataset_results(grouped_by_dataset)
    except:
        pass

    print("X9", time.time() - starttime)
    timing_list.append({"desc" : "Enrich datasets", "timing": time.time() - starttime})
    
    return merged_df, timing_list, grouped_by_unit_delta_mass, grouped_by_dataset, grouped_by_dataset_unit_delta, metadata_groupings, log_txt

def main():
    parser = argparse.ArgumentParser(description='Build Query File')
    parser.add_argument('--input_mgf', help='Input MGF File')
    parser.add_argument('--spectral_threshold', help='Spectral Threshold', type=str)
    args = parser.parse_args()
    
    input_mgf = args.input_mgf
    spectral_threshold = ast.literal_eval(args.spectral_threshold)
    
    # Nextflow workaround
    if len(spectral_threshold) != 1:
        raise ValueError(f'Spectral Threshold must be a list of length 1 but got {spectral_threshold}')
    spectral_threshold=spectral_threshold[0]
    
    search_library(input_mgf,
                    "input_library",               # ./libraries/input_library
                    analog_search=False,
                    lower_delta=130,
                    upper_delta=200,
                    pm_tolerance=0.2,
                    fragment_tolerance=0.2,
                    cosine_threshold=spectral_threshold)
    
if __name__ == "__main__":
    main()