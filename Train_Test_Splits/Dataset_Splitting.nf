#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Workflow Boiler Plate
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

TOOL_FOLDER_LS = "$moduleDir/bin"
FAST_SEARCH_LIBRARY_BIN = "$TOOL_FOLDER_LS/GNPS_FastSearch_Library/bin"

params.input_csv = '/home/user/SourceCode/GNPS_ML_Processing_Workflow/GNPS_ML_Processing/nf_output/ALL_GNPS_cleaned.csv'
params.input_mgf = '/home/user/SourceCode/GNPS_ML_Processing_Workflow/GNPS_ML_Processing/nf_output/ALL_GNPS_cleaned.mgf'

params.test_set_num = 0.10  // An integer (or float) representing the number (or percentage of) data points to use as a test set

params.lowest_spectral_threshold = '0.7' 
params.lowest_structural_threshold = '0.7' 

// Library Search Parameters
params.library_index_path = 
params.temp_results_tsv   = 
params.pm_tolerance       = 0.2
params.fragment_tolerance = 0.2
params.lower_delta        = 130  // To perform analog search, set lower_delta = upper_delta =0
params.upper_delta        = 200

/*
 * Generates a random subset of params.test_set_num  data points to use as a test set
 * These data points should be consistent across splitting methodology
 */
process generate_test_set {
  conda "$TOOL_FOLDER_LS/conda_env.yml"

  input:
  path csv_file
  path mgf_file

  output:
  path 'test_rows.csv',  emit: test_rows_csv
  path 'test_rows.mgf',  emit: test_rows_mgf
  path 'train_rows.csv', emit: train_rows_csv
  path 'train_rows.mgf', emit: train_rows_mgf

  """
  python3 $TOOL_FOLDER_LS/calc_test.py \
          --input_csv "$csv_file" \
          --input_mgf "$mgf_file" \
          --num_test_points $params.test_set_num
  """
}

// Build a library for the test set to search against
process build_library {
  conda "$TOOL_FOLDER_LS/conda_env.yml"

  input:
  path mgf_file

  output:
  path 'libraries', emit: libraries

  """
  ret_value=(\$(python3 $TOOL_FOLDER_LS/build_library.py --input_mgf "$mgf_file" --output_folder "./libraries"))
  build_input_folder=\${ret_value[0]}
  build_index_folder=\${ret_value[1]}
  spectrum_level_metadata_path=\${ret_value[2]}
  file_level_metadata_path=\${ret_value[3]}

  echo "Build Input Folder: " \$build_input_folder
  echo "Build Index Folder: " \$build_index_folder
  echo "Spectrum Level Metadata Path: " \$spectrum_level_metadata_path
  echo "File Level Metadata Path: " \$file_level_metadata_path

  $FAST_SEARCH_LIBRARY_BIN/main_execmodule ExecIndex $FAST_SEARCH_LIBRARY_BIN/generic_params \
  -autoINDEX_TYPE mxc \
  -autoINPUT_SPECS  \$build_input_folder \
  -autoOUTPUT_INDEX  \$build_index_folder \
  -autoSNR 0 \
  -autoMIN_MZ_PEAK 0 \
  -autoWINDOW_FILTER_RANK 6 \
  -autoWINDOW_FILTER_SIZE 50 \
  -autoMIN_PEAKS_AFTER_PROCESSING 1 \

  $FAST_SEARCH_LIBRARY_BIN/main_execmodule ExecIndex $FAST_SEARCH_LIBRARY_BIN/generic_params \
  -autoINPUT_INDEX  \$build_index_folder \
  -autoINPUT_METADATA_ANNOTATIONS \$spectrum_level_metadata_path \
  -autoINPUT_METADATA_FILENAMES \$file_level_metadata_path
  """
}

// Perform search against the library for similarities
process spectral_similarity_split {
  conda "$TOOL_FOLDER_LS/conda_env.yml"

  input:
  path library_dir
  path comparison_mgf

  """
  ret_value=(\$(python3 $TOOL_FOLDER_LS/build_query_file.py --input_mgf "$comparison_mgf")
  temp_query_mgf=\${ret_value[0]}

  echo "Query MGF File Path: " \$temp_query_mgf
  echo "Minimum Cosine Threshold: " \$params.lowest_spectral_threshold

 $FAST_SEARCH_LIBRARY_BIN/main_execmodule ExecIndex $FAST_SEARCH_LIBRARY_BIN/generic_params \
  -autoINPUT_INDEX ${library_dir}/build_index \
  -autoOUTPUT_RESULTS output.csv  \
  -autoINPUT_SPECS \$temp_query_mgf \
  -autoPM_TOLERANCE $params.pm_tolerance  \
  -autoFRAG_TOLERANCE $params.fragment_tolerance  \
  -autoDELTA_MZ_ABOVE $params.lower_delta\
  -autoDELTA_MZ_BELOW $params.upper_delta \
  -autoTHETA \$params.lowest_spectral_threshold \
  -autoSPECTRUM_DISK_ACCESS_POLICY DISK   \
  -autoINDEX_DISK_ACCESS_POLICY DISK  \
  -autoOUTPUT_FILTERED_QUERIES_JSON_FOLDER temp_results_json_folder \
  -autoDELTA_F 0.4 \
  -autoVALIDATE 0
  """
}

workflow {
  csv_file = Channel.fromPath(params.input_csv)
  mgf_file = Channel.fromPath(params.input_mgf)

  generate_test_set(csv_file, mgf_file)
  build_library(generate_test_set.out.test_rows_mgf)

  spectral_similarity_split(build_library.out.libraries, generate_test_set.out.train_rows_mgf)
}