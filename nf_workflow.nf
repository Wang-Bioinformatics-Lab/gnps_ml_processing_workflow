#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// params.subset = "Bruker_Fragmentation_Prediction"
params.subset = "MH_MNA_Translation"

params.parallelism = 10
// Workflow Boiler Plate
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

TOOL_FOLDER = "$baseDir/bin"
GNPS_EXPORTS_FOLDER = "$baseDir/GNPS_ml_exports"

process export {
    conda "$TOOL_FOLDER/conda_env.yml"

    output:
    path './temp', emit: temp_path

    """
    python3 $TOOL_FOLDER/GNPS2_Processor.py -p "$params.parallelism"  
    """
}

process postprocess {
  conda "$TOOL_FOLDER/conda_env.yml"

  input:
  path temp_path

  output:
  path "ALL_GNPS_cleaned.csv", emit: cleaned_csv
  path "ALL_GNPS_cleaned.parquet", emit: cleaned_parquet

  """
  python3 $TOOL_FOLDER/GNPS2_Postprocessor.py 
  """
}

process generate_subset {
  publishDir "./nf_output", mode: 'copy'

  conda "$TOOL_FOLDER/conda_env.yml"

  input:
  path cleaned_csv
  path cleaned_parquet

  output:
  path "ALL_GNPS_cleaned.csv"
  path "ALL_GNPS_cleaned.parquet"
  path "summary*.csv"
  path "spectra*.parquet"

  """
  python3 $TOOL_FOLDER/GNPS2_Subset_Generator.py "$params.subset"
  """
}

workflow {
  export()
  postprocess(export.out.temp_path)
  generate_subset(postprocess.out.cleaned_csv, postprocess.out.cleaned_parquet)    
}