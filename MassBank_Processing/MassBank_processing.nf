#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.massbank_url = "https://github.com/MassBank/MassBank-data.git"

// Workflow Boiler Plate
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

TOOL_FOLDER_MB = "$moduleDir/bin"

// Use git to pull the latest MassBank Library
process fetch_data_massbank {
    output: 
    path "MassBank-data/MSSJ/*.txt", emit: massbank_data  // for testing
    // path "MassBank-data/*/*.txt", emit: massbank_data
    path "MassBank-data/legacy.blacklist", emit: blacklist

    """
    git clone $params.massbank_url
    """
}

// Process all data into a unified csv, mgf format
process export_massbank {
    conda "$TOOL_FOLDER_MB/conda_env.yml"

    input: 
    each input_file
    path blacklist

    // Optional because MS1 scans will not be read, and blacklist files are skipped
    output:
    path 'output_*', optional: true

    """
    python3 $TOOL_FOLDER_MB/MassBank_processing.py \
    --input "$input_file" \
    --blacklist "$blacklist"
    """
}

// Merges all the export_massbanks together
process merge_export_massbank {
  conda "$TOOL_FOLDER_MB/conda_env.yml"

  publishDir "./MassBank_export_massbank", mode: 'copy'

  input:
  path temp_files

  output:
  path 'ALL_MassBank_merged.*', emit: merged_files

  """
  python3 $TOOL_FOLDER_MB/merge_files.py 
  """
}

workflow {
  fetch_data_massbank()
  temp_files = export_massbank(fetch_data_massbank.out.massbank_data, fetch_data_massbank.out.blacklist)
  (merged_mgf, merged_csv) = merge_export_massbank(temp_files.collect())
}