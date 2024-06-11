#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.massbank_url = "https://github.com/MassBank/MassBank-data.git"

// Workflow Boiler Plate
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

TOOL_FOLDER_MB = "$moduleDir/bin"

/*
 * Number of bins to bin MassBank data into.
 * Note that this is limited by the number of file handles that your system
 * can tolerate. If you are running into issues, try reducing this number.
 * If you run into load balancing issues, increase this number.
*/
params.spectra_parallelism = 5000

// Use git to pull the latest MassBank Library
process fetch_data_massbank {
    output: 
    path "MassBank-data/", emit: massbank_data
    path "MassBank-data/legacy.blacklist", emit: blacklist

    """
    git clone $params.massbank_url
    """
}

// Bin the massbank data into params.spectra_parallelism bins
process prep_params_massbank {
  conda "$TOOL_FOLDER_MB/conda_env.yml"

  input:
  path blacklist
  path massbank_data

  output:
  path 'params/params_*.npy', emit: params

  """
  python3 $TOOL_FOLDER_MB/prep_params.py \
            --glob "MassBank-data/*/*.txt" \
            --blacklist "$blacklist" \
            -p "$params.spectra_parallelism"

  # python3 $TOOL_FOLDER_MB/prep_params.py \
  #         --glob "MassBank-data/MSSJ/*.txt" \
  #         --blacklist "$blacklist" \
  #         -p "$params.spectra_parallelism"
  """
}

// Process all data into a unified csv, mgf format
process export_massbank {
    conda "$TOOL_FOLDER_MB/conda_env.yml"

    cache true

    input: 
    path masbank_data
    each params

    // Optional because MS1 scans will not be read, and blacklist files are skipped
    output:
    path 'output_*', optional: true

    """
    python3 "$TOOL_FOLDER_MB/MassBank_processing.py" \
    --params "$params"
    """
}

// Merges all the export_massbanks together
process merge_export_massbank {
  conda "$TOOL_FOLDER_MB/conda_env.yml"

  input:
  path temp_files

  output:
  path 'ALL_MassBank_merged.*', emit: merged_files

  """
  python3 $TOOL_FOLDER_MB/merge_files.py 
  """
}

// Output the data into nf_output
process output_files {
  publishDir "./nf_output", mode: 'copy'

  input:
  path merged_files

  output:
  path 'ALL_MassBank_merged.*', emit: merged_files, includeInputs: true

  """
  echo "Outputting data to nf_output"
  """
}

workflow {
  fetch_data_massbank()
  prep_params_massbank(fetch_data_massbank.out.blacklist, fetch_data_massbank.out.massbank_data)
  temp_files = export_massbank(fetch_data_massbank.out.massbank_data, prep_params_massbank.out.params)
  merged_files = merge_export_massbank(temp_files.collect())
  output_files(merged_files)
}