#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// params.subset = "Bruker_Fragmentation_Prediction"
params.subset = "MH_MNA_Translation"
params.split  = true

params.spectra_parallelism = 10
// Workflow Boiler Plate
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

TOOL_FOLDER = "$baseDir/bin"
GNPS_EXPORTS_FOLDER = "$baseDir/GNPS_ml_exports"

params.parallelism = 24

process export {
    conda "$TOOL_FOLDER/conda_env.yml"

    output:
    path './temp', emit: temp_path

    """
    python3 $TOOL_FOLDER/GNPS2_Processor.py -p "$params.spectra_parallelism"  
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
  path "summary/*.csv"
  path "spectra/*.parquet", emit: output_parquet

  """
  python3 $TOOL_FOLDER/GNPS2_Subset_Generator.py "$params.subset"
  """
}

process generate_mgf {
  conda "$TOOL_FOLDER/conda_env.yml"

  input:
  path output_parquet

  output:
  path "*.mgf", emit: output_mgf

  """
  python3 $TOOL_FOLDER/GNPS2_MGF_Generator.py -p "$params.parallelism"
  """
}

process calculate_similarities {
  input:
  path mgf

  output:
  path "similarity_calculations/*", emit: spectral_similarities

  """
  alias python=python3
  nextflow run $baseDir/GNPS_PureNetworking_Workflow/workflow/workflow.nf \
                --inputspectra $mgf \
                --parallelism $params.parallelism \
                --publishdir "similarity_calculations/${mgf.baseName}"
  """
}

process split_subsets {
  conda "$TOOL_FOLDER/conda_env.yml"
  publishDir "./nf_output", mode: 'copy'

  cache false

  input:
  path spectral_similarities

  // output: 
  // path "subsets/*"

  // exec:
  // println "$spectral_similarities"

  // Want to split CSV, parquet
  // Requires CSV, parquet, similarities
  """
  python3 $TOOL_FOLDER/GNPS2_Subset_Split.py  ${spectral_similarities}/merged_pairs.tsv \
                                              $baseDir/nf_output/spectra/${spectral_similarities}.parquet \
                                              $baseDir/nf_output/summary/${spectral_similarities}.csv
  """
}

workflow {
  export()
  postprocess(export.out.temp_path)
  generate_subset(postprocess.out.cleaned_csv, postprocess.out.cleaned_parquet)    
  if ("$params.split") {
    generate_mgf(generate_subset.out.output_parquet)
    generate_mgf.out.output_mgf  | calculate_similarities 
    calculate_similarities.out.spectral_similarities | split_subsets
  }
}