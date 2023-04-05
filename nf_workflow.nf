#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// params.subset = "Bruker_Fragmentation_Prediction"
// params.subset = "MH_MNA_Translation"
params.subset = "GNPS_default"
params.split  = false
use_default_path = true

params.spectra_parallelism = 100
// Workflow Boiler Plate
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

TOOL_FOLDER = "$baseDir/bin"
GNPS_EXPORTS_FOLDER = "$baseDir/GNPS_ml_exports"

params.parallelism = 12

process prep_params {
  conda "$TOOL_FOLDER/conda_env.yml"

  output:
  path 'params/params_*.npy', emit: export_params

  """
  python3 $TOOL_FOLDER/prep_params.py -p "$params.spectra_parallelism"  
  """
}

process export {
    conda "$TOOL_FOLDER/conda_env.yml"

    input: 
    each input_file

    output:
    path 'temp/*', emit: temp_files

    """
    python3 $TOOL_FOLDER/GNPS2_Processor.py -f $input_file
    """
}

process merge_export {

  input:
  path temp_files

  output:
  path './ALL_GNPS_merged.mgf', emit: merged_csv
  path './ALL_GNPS_merged.csv', emit: merged_mgf

  """
  python3 $TOOL_FOLDER/merge_files.py 
  """

}

process postprocess {
  conda "$TOOL_FOLDER/conda_env.yml"

  cache true

  input:
  path merged_csv
  path merged_parquet

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
  path "summary/*"
  path "spectra/*", emit: output_parquet
  path "util/*", optional: true

  """
  python3 $TOOL_FOLDER/GNPS2_Subset_Generator.py "$params.subset"
  """
}

process generate_mgf {
  conda "$TOOL_FOLDER/conda_env.yml"

  input:
  path dummy

  output:
  path "*.mgf", emit: output_mgf, optional: true

  // The default path will process all files while the other path will process only the spectral_similarity_prediction subset
  if (use_default_path) {
    """
    python3 $TOOL_FOLDER/GNPS2_MGF_Generator.py -p "$params.parallelism"
    """
  } else {
    """
    python3 $TOOL_FOLDER/GNPS2_MGF_Generator.py -p "$params.parallelism" -path "./Spectral_Similarity_Prediction.parquet"
    """
  }
  """
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

process publish_similarities_for_prediction {
  publishDir "./nf_output/", mode: 'copy'

  input:
  path inputFile

  output:
  path "./util/Spectral_Similarity_Prediction_Similarities.tsv"

  """
  cp ${inputFile} ./util/Spectral_Similarity_Prediction_Similarities.tsv
  """
}

workflow {
  prep_params()
  export(prep_params.out.export_params)
  merge_export(export.out.temp_files.collect())
  postprocess(merge_export.out.merged_csv, merge_export.out.merged_mgf)
  generate_subset(postprocess.out.cleaned_csv, postprocess.out.cleaned_parquet)    
  // For the spectral similarity prediction task, we need to calculate all pairs similarity in the training set
  if (params.subset == "GNPS_default" || params.subset == "Spectral_Similarity_Prediction") {
    use_default_path = false
    generate_mgf(generate_subset.out.output_parquet)
    generate_mgf.out.output_mgf  | calculate_similarities 
    publish_similarities_for_prediction(calculate_similarities.out.spectral_similarities)
  }

  if (params.split) {
    generate_mgf(generate_subset.out.output_parquet)
    generate_mgf.out.output_mgf  | calculate_similarities 
    calculate_similarities.out.spectral_similarities | split_subsets
  }
}