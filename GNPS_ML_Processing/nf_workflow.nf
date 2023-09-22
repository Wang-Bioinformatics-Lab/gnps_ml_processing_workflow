#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// params.subset = "Orbitrap_Fragmentation_Prediction"
params.split  = false

params.subset = "Structural_Similarity_Prediction"
// params.subset = "MH_MNA_Translation"
// params.subset = "GNPS_default"

use_default_path = true

params.spectra_parallelism = 100

params.path_to_provenance = "/home/user/LabData/GNPS_Library_Provenance/"

// Workflow Boiler Plate
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

TOOL_FOLDER = "$baseDir/bin"
params.parallelism = 12
params.pure_networking_parallelism = 5000
params.pure_networking_forks = 32

process prep_params {
  conda "$TOOL_FOLDER/conda_env.yml"

  cache 'lenient'

  output:
  path 'params/params_*.npy'

  """
  python3 $TOOL_FOLDER/prep_params.py -p "$params.spectra_parallelism"  
  """
}

process export {
    conda "$TOOL_FOLDER/conda_env.yml"

    maxForks 16

    input: 
    each input_file

    output:
    path 'temp/*'

    """
    python3 $TOOL_FOLDER/GNPS2_Processor.py -f "$input_file" --path_to_provenance "$params.path_to_provenance"
    """
}

process merge_export {
  //conda "$TOOL_FOLDER/gnps_ml_processing_env2/"
  conda "$TOOL_FOLDER/conda_env.yml"

  input:
  path temp_files

  output:
  path './ALL_GNPS_merged.mgf'
  path './ALL_GNPS_merged.csv'

  """
  python3 $TOOL_FOLDER/merge_files.py 
  """

}

process postprocess {
  //conda "$TOOL_FOLDER/gnps_ml_processing_env2/"
  conda "$TOOL_FOLDER/conda_env.yml"

  publishDir "./nf_output", mode: 'copy'

  cache true

  input:
  path merged_csv
  path merged_parquet
  path adduct_mapping

  output:
  path "ALL_GNPS_cleaned.csv", emit: cleaned_csv
  path "ALL_GNPS_cleaned.parquet", emit: cleaned_parquet
  path "ALL_GNPS_cleaned.mgf", emit: cleaned_mgf

  script:
  """
  Rscript -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org")'
  Rscript -e 'BiocManager::install("Rdisop")'
  Rscript -e "devtools::install_github('mjhelf/MassTools')"
  """

  """
  python3 $TOOL_FOLDER/GNPS2_Postprocessor.py 
  """
}

process export_full_json {
  //conda "$TOOL_FOLDER/gnps_ml_processing_env2/"
  conda "$TOOL_FOLDER/conda_env.yml"

  cache true

  input:
  path cleaned_csv
  path cleaned_mgf

  output:
  // path "*.json", emit: cleaned_json
  path "json_outputs/*.json"
  val 1, emit: dummy

  """
  python3 $TOOL_FOLDER/GNPS2_JSON_Export.py --input_mgf_path ${cleaned_mgf}\
                                            --input_csv_path ${cleaned_csv}\
                                            --output_path "./json_outputs" \
                                            --progress \
                                            --per_library
  """
}

process generate_subset {
  conda "$TOOL_FOLDER/gnps_ml_processing_env2/"
 
  publishDir "./nf_output", mode: 'copy'

  cache true

  input:
  path cleaned_csv
  path cleaned_parquet
  path cleaned_mgf
  val json_dummy          // Dummy input to force this process to run after the export_full_json process

  output:
  path "summary/*"
  path "spectra/*.parquet", emit: output_parquet
  path "spectra/*.mgf", emit: output_mgf
  path "json_outputs/*.json", emit: output_json, optional: true
  path "util/*", optional: true

  """
  python3 $TOOL_FOLDER/GNPS2_Subset_Generator.py "$params.subset"
  """
}

process generate_mgf {
  conda "$TOOL_FOLDER/gnps_ml_processing_env2/"

  input:
  each parquet_file

  output:
  path "*.mgf", emit: output_mgf//, optional: true

  // The default path will process all files while the other path will process only the spectral_similarity_prediction subset
  //old if 
  //python3 $TOOL_FOLDER/GNPS2_MGF_Generator.py -p "$params.parallelism"
  //old else
  // // python3 $TOOL_FOLDER/GNPS2_MGF_Generator.py -p "$params.parallelism" -path "./Spectral_Similarity_Prediction.parquet"
  
  // if (use_default_path) {
  //   """ 
  //   python3 $TOOL_FOLDER/GNPS2_MGF_Generator.py -input_path $parquet_file -output_path ${name}.mgf"
  //   """
  // } else {
  //   """
  //   python3 $TOOL_FOLDER/GNPS2_MGF_Generator.py -input_path "./Spectral_Similarity_Prediction.parquet" -output_path "./Spectral_Similarity_Prediction.mgf"
  //   """
  // }
  """
  python3 $TOOL_FOLDER/GNPS2_MGF_Generator.py -input_path "$parquet_file"
  """
}

process calculate_similarities_pure_networking {
  // Currently, this process is not used and it is scheduled for deletion
  publishDir "./nf_output", mode: 'copy'
  
  input:
  each mgf

  output:
  path "similarity_calculations/*", emit: spectral_similarities

  """
  nextflow run $TOOL_FOLDER/GNPS_PureNetworking_Workflow/workflow/workflow.nf \
                --inputspectra $mgf \
                --parallelism $params.pure_networking_parallelism \
                --maxforks $params.pure_networking_forks \
                --publishdir "similarity_calculations/${mgf.baseName}"
  """
}

process calculate_similarities {
  // Similarities using fasst search have not been implemented yet
  //conda "$TOOL_FOLDER/gnps_ml_processing_env2/"
  conda "$TOOL_FOLDER/conda_env.yml"

  publishDir "./nf_output", mode: 'copy'
  
  input:
  each mgf

  output:
  path "similarity_calculations/*", emit: spectral_similarities

  """
  exit()
  nextflow run $TOOL_FOLDER/GNPS_PureNetworking_Workflow/workflow/workflow.nf \
                --inputspectra $mgf \
                --parallelism $params.pure_networking_parallelism \
                --maxforks $params.pure_networking_forks \
                --publishdir "similarity_calculations/${mgf.baseName}"
  """
}

process split_subsets {
  conda "$TOOL_FOLDER/gnps_ml_processing_env2/"

  publishDir "./nf_output", mode: 'copy'

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
  // I'm not sure if there's a better way to skip the files we don't need, for now we'll just ignore the errors
  publishDir "./nf_output/", mode: 'copy'
  errorStrategy 'ignore'

  input:
  path inputFiles

  output:
  path "./util/Spectral_Similarity_Prediction_Similarities.tsv", optional: true

  """
  [ -f "Spectral_Similarity_Prediction/merged_pairs.tsv" ] && mkdir util && cp Spectral_Similarity_Prediction/merged_pairs.tsv ./util/Spectral_Similarity_Prediction_Similarities.tsv
  """
}

workflow {
  export_params = prep_params()
  temp_files = export(export_params)
  (merged_mgf, merged_csv) = merge_export(temp_files.collect())

  // A python dictionary that maps GNPS adducts to a unified set of adducts used in the GNPS2 workflow
  adduct_mapping_ch = channel.fromPath("$TOOL_FOLDER/adduct_mapping.pkl")
  
  postprocess(merged_csv, merged_mgf, adduct_mapping_ch)
  export_full_json(postprocess.out.cleaned_csv, postprocess.out.cleaned_mgf)
  generate_subset(postprocess.out.cleaned_csv, postprocess.out.cleaned_parquet, postprocess.out.cleaned_mgf, export_full_json.out.dummy)    

  calculate_similarities(generate_subset.out.output_mgf)

  // For the spectral similarity prediction task, we need to calculate all pairs similarity in the training set
  if (params.subset == "GNPS_default" || params.subset == "Spectral_Similarity_Prediction") {
    use_default_path = false
    //// generate_subset.out.output_parquet.collect()  // Make sure this finishes first
    //// channel.fromPath(".nf_output/Spectral_Similarity_Prediction.parquet") | generate_mgf
    //generate_mgf(generate_subset.out.output_parquet)
    //calculate_similarities(generate_mgf.out.output_mgf)

    publish_similarities_for_prediction(calculate_similarities.out.spectral_similarities)

  }
  // } else if (params.subset == "Spectral_Similarity_Prediction") {
  //   generate_subset.out.output_parquet.collect()  // Make sure this finishes first
  //   channel.fromPath(".nf_output/Spectral_Similarity_Prediction.parquet") | generate_mgf
  //   // generate_mgf(generate_subset.out.output_parquet)
  //   generate_mgf.out.output_mgf  | calculate_similarities 
  //   publish_similarities_for_prediction(calculate_similarities.out.spectral_similarities)
  // }

  if (params.split) {
    calculate_similarities.out.spectral_similarities | split_subsets
  }
}