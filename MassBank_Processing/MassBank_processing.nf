#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.massbank_url = "https://github.com/MassBank/MassBank-data.git"

// Workflow Boiler Plate
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

TOOL_FOLDER = "$baseDir/bin"

// Use git to pull the latest MassBank Library
process fetch_data {
    output: 
    // path "MassBank-data/MSSJ/*.txt", emit: massbank_data  // for testing
    path "MassBank-data/*/*.txt", emit: massbank_data
    path "MassBank-data/legacy.blacklist", emit: blacklist

    """
    git clone $params.massbank_url
    """
}

// Process all data into a unified csv, mgf format
process export {
    conda "$TOOL_FOLDER/conda_env.yml"

    cache true

    input: 
    each input_file
    path blacklist

    // Optional because MS1 scans will not be read, and blacklist files are skipped
    output:
    path 'output_*', optional: true

    """
    python3 $TOOL_FOLDER/MassBank_processing.py \
    --input "$input_file" \
    --blacklist "$blacklist"
    """
}

// Merges all the exports together
process merge_export {
  conda "$TOOL_FOLDER/conda_env.yml"

  publishDir "./MassBank_Export", mode: 'copy'

  cache false
  input:
  path temp_files

  output:
  path './ALL_MassBank_merged.mgf'
  path './ALL_MassBank_merged.csv'

  """
  python3 $TOOL_FOLDER/merge_files.py 
  """
}

// Allows calls from other workflows
process export_massbank {
  fetch_data()
  temp_files = export(fetch_data.out.massbank_data, fetch_data.out.blacklist)
  (merged_mgf, merged_csv) = merge_export(temp_files.collect())
}

workflow {
  fetch_data()
  temp_files = export(fetch_data.out.massbank_data, fetch_data.out.blacklist)
  (merged_mgf, merged_csv) = merge_export(temp_files.collect())
}