#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Workflow Boiler Plate
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

TOOL_FOLDER_LS = "$moduleDir/bin"
FAST_SEARCH_LIBRARY_BIN = "$TOOL_FOLDER_LS/GNPS_FastSearch_Library/bin"

// params.input_csv = '/home/user/SourceCode/GNPS_ML_Processing_Workflow/GNPS_ML_Processing/nf_output/summary/Structural_Similarity_Prediction.csv'
// params.input_mgf = '/home/user/SourceCode/GNPS_ML_Processing_Workflow/GNPS_ML_Processing/nf_output/spectra/Structural_Similarity_Prediction.mgf'
params.input_csv = "../GNPS_ML_Processing/nf_output/ML_ready_subset_positive/selected_summary.csv"
params.input_mgf = "../GNPS_ML_Processing/nf_output/ML_ready_subset_positive/selected_spectra.mgf"

// Which task subset to use
params.subset = "Structural_Similarity_Prediction"
params.split_type = "structure_smart"  // 'structure_smart', 'random' (random spectra), or 'structure' (random inchi14)

params.test_set_num = 2300  // An integer (or float) representing the number (or percentage of) data points to use as a test set

params.lowest_spectral_threshold = '0.6' 
params.lowest_structural_threshold = '0.2' 
params.thresholds = "0.2 0.99"
params.structural_similarity_fingerprint = "Morgan_2048_3"
params.ion_mode = "positive"  // 'positive' or 'negative'

// Library Search Parameters
params.pm_tolerance       = 0.2
params.fragment_tolerance = 0.2
params.lower_delta        = 130  // To perform analog search, set lower_delta = upper_delta =0
params.upper_delta        = 200

/*
 * This process selects the subset of data that is usable for ML (e.g., more than one peak, not excessively noisy, etc.)
 * It also splits the data into the selected ion mode, and normalizes intensities.
 */
// process select_data_for_ml {
//   conda "$TOOL_FOLDER_LS/conda_env.yml"

//   input:
//   path csv_file
//   path mgf_file

//   output:
//   path 'selected_summary.csv',  emit: selected_summary
//   path 'selected_spectra.mgf',  emit: selected_spectra

//   """
//   python3 $TOOL_FOLDER_LS/select_data.py \
//           --input_csv "$csv_file" \
//           --input_mgf "$mgf_file" \
//           --ion_mode $params.ion_mode
//   """
// }

process generate_subset {
  conda "$TOOL_FOLDER_LS/conda_env.yml"
 
  publishDir "./nf_output", mode: 'copy'

  cache true

  input:
  path cleaned_csv
  path cleaned_mgf

  output:
  path "summary/*", emit: output_summary
  path "spectra/*.parquet", emit: output_parquet, optional: true
  path "spectra/*.mgf", emit: output_mgf
  path "json_outputs/*.json", emit: output_json, optional: true
  path "util/*", optional: true

  """
  python3 $TOOL_FOLDER_LS/subset_generator.py "$params.subset"
  """
}

/*
 * Generates a random subset of params.test_set_num  data points to use as a test set
 * These data points should be consistent across splitting methodology
 */
process generate_test_set {
  conda "$TOOL_FOLDER_LS/conda_env.yml"
  publishDir "./nf_output/${params.subset}/${params.split_type}", mode: 'copy'

  cache false

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
          --num_test_points $params.test_set_num \
          --sampling_strategy "$params.split_type" \
          --threshold 0.5 \
          --debug
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
process spectral_similarity_calculation {
  conda "$TOOL_FOLDER_LS/conda_env.yml"

  input:
  path library_dir
  path comparison_mgf

  output:
  path 'output.tsv', emit: spectral_similarities

  """
  ret_value=(\$(python3 $TOOL_FOLDER_LS/build_query_file.py --input_mgf "$comparison_mgf"))
  temp_query_mgf=\${ret_value[0]}

  echo "Query MGF File Path: " \$temp_query_mgf
  echo "Minimum Cosine Threshold: " $params.lowest_spectral_threshold

 $FAST_SEARCH_LIBRARY_BIN/main_execmodule ExecIndex $FAST_SEARCH_LIBRARY_BIN/generic_params \
  -autoINPUT_INDEX ${library_dir}/build_index \
  -autoOUTPUT_RESULTS output.tsv  \
  -autoINPUT_SPECS \$temp_query_mgf \
  -autoPM_TOLERANCE $params.pm_tolerance  \
  -autoFRAG_TOLERANCE $params.fragment_tolerance  \
  -autoDELTA_MZ_ABOVE $params.lower_delta \
  -autoDELTA_MZ_BELOW $params.upper_delta \
  -autoTHETA $params.lowest_spectral_threshold \
  -autoSPECTRUM_DISK_ACCESS_POLICY DISK   \
  -autoINDEX_DISK_ACCESS_POLICY DISK  \
  -autoOUTPUT_FILTERED_QUERIES_JSON_FOLDER temp_results_json_folder \
  -autoDELTA_F 0.4 \
  -autoVALIDATE 0
  """
}

process structural_similarity_calculation {
  conda "$TOOL_FOLDER_LS/conda_env.yml"
  publishDir './nf_output', mode: 'copy'

  cache true

  input:
  path train_rows_csv
  path test_rows_csv

  output:
  path 'output.csv', emit: structural_similarities

  """
  python3 $TOOL_FOLDER_LS/structural_similarity.py \
          --train_csv "$train_rows_csv" \
          --test_csv "$test_rows_csv" \
          --similarity_threshold $params.lowest_structural_threshold \
          --fingerprint "$params.structural_similarity_fingerprint"
  """
}

process split_data {
  publishDir './nf_output', mode: 'copy'
  conda "$TOOL_FOLDER_LS/conda_env.yml"

  cache false

  // Gather inputs avoiding name collisions
  input:
  path train_rows_csv, stageAs: 'train_rows.csv'
  path test_rows_csv, stageAs: 'test_rows.csv'
  path train_rows_mgf, stageAs: 'train_rows.mgf'
  path test_rows_mgf, stageAs: 'test_rows.mgf'
  path spectral_similarities
  path structural_similarities

  output:
  path 'train_rows_*.csv'
  path 'train_rows_*.mgf'
  path 'test_rows_*.csv'
  path 'test_rows_*.mgf'
  
  """
  python3 $TOOL_FOLDER_LS/split_data.py \
          --input_train_csv $train_rows_csv \
          --input_train_mgf $train_rows_mgf \
          --input_test_csv $test_rows_csv \
          --input_test_mgf $test_rows_mgf \
          --spectral_similarities "$spectral_similarities" \
          --strucutral_similarities "$structural_similarities" \
          --similarity_thresholds $params.thresholds  # No quotes around this param is necessary
  """
}

// Outputs the results (avoids seding data to output if proceses are included in another script)
process output_handler {
  publishDir "./nf_output", mode: 'copy'

  input:
  path train_rows_csv
  path test_rows_csv
  path train_rows_mgf
  path test_rows_mgf

  output: 
  path "*", includeInputs: true

  """
  echo "Outputting data to nf_output"
  """
}

workflow {
  csv_file = Channel.fromPath(params.input_csv)
  mgf_file = Channel.fromPath(params.input_mgf)

  // select_data_for_ml(csv_file, mgf_file)
  // generate_test_set(select_data_for_ml.out.selected_summary, select_data_for_ml.out.selected_spectra)


  generate_subset(csv_file, mgf_file)
  generate_test_set(generate_subset.out.output_summary, generate_subset.out.output_mgf)
  // build_library(generate_test_set.out.test_rows_mgf)

  // spectral_similarity_calculation(build_library.out.libraries, generate_test_set.out.train_rows_mgf)

  // structural_similarity_calculation(generate_test_set.out.train_rows_csv, 
                                    // generate_test_set.out.test_rows_csv)

  // split_data(generate_test_set.out.train_rows_csv, 
  //           generate_test_set.out.test_rows_csv, 
  //           generate_test_set.out.train_rows_mgf, 
  //           generate_test_set.out.test_rows_mgf,
  //           spectral_similarity_calculation.out.spectral_similarities,
  //           structural_similarity_calculation.out.structural_similarities)

  // output_handler(split_data.out.train_rows_csv_spectral, 
  //               split_data.out.test_rows_csv_spectral, 
  //               split_data.out.train_rows_mgf_spectral, 
  //               split_data.out.test_rows_mgf_spectral)
}