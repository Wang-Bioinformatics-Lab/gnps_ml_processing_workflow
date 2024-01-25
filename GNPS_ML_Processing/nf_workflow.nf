#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// params.subset = "Orbitrap_Fragmentation_Prediction"
params.split  = false

params.subset = "Structural_Similarity_Prediction"
// params.subset = "MH_MNA_Translation"
// params.subset = "GNPS_default"

params.output_dir = "./nf_output"
params.GNPS_json_path = "None"

use_default_path = true

params.spectra_parallelism = 5000

params.path_to_provenance = "/home/user/LabData/GNPS_Library_Provenance/"

// If true, will download and reparse massbank, additionally removing all massbank entires from GNPS
params.include_massbank = true

// Workflow Boiler Plate
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

TOOL_FOLDER = "$baseDir/bin"
params.parallelism = 12
params.pure_networking_parallelism = 5000
params.pure_networking_forks = 32


// Include MassBank parser in case params.include_massbank is true
include { fetch_data_massbank; prep_params_massbank; export_massbank; merge_export_massbank } from '../MassBank_Processing/MassBank_processing.nf'

// This environment won't build if called by nextflow, but works fine here
process environment_creation {
  
  output: 
  path "dummy.text", emit: dummy
  
  """
  mamba env update --file $TOOL_FOLDER/conda_env.yml --prefix $TOOL_FOLDER/gnps_ml_processing_env2/

  touch dummy.text
  """
}

// Splitting out the GNPS Libraries into smaller chunks
process prep_params {
  conda "$TOOL_FOLDER/gnps_ml_processing_env2"

  cache 'lenient'

  input:
  path dummy

  output:
  path 'params/params_*.npy'

  """
  python3 $TOOL_FOLDER/prep_params.py \
  -p "$params.spectra_parallelism" \
  --all_GNPS_json_path $params.GNPS_json_path
  """
}

// Pull additional data from Provenance File
process export {
    conda "$TOOL_FOLDER/gnps_ml_processing_env2"

    maxForks 8
    errorStrategy 'retry'

    input: 
    each input_file

    output:
    path 'temp/*'

    """
    python3 $TOOL_FOLDER/GNPS2_Processor.py \
    -f "$input_file" \
    --path_to_provenance "$params.path_to_provenance"
    """
}

// Merges all the exports together
process merge_export {
  //conda "$TOOL_FOLDER/gnps_ml_processing_env2/"
  conda "$TOOL_FOLDER/gnps_ml_processing_env2"

  input:
  path temp_files

  output:
  path './ALL_GNPS_merged.mgf'
  path './ALL_GNPS_merged.csv'

  script:
  if (params.include_massbank)
    """
    python3 $TOOL_FOLDER/merge_files.py --include_massbank 
    """
  else
    """
    python3 $TOOL_FOLDER/merge_files.py 
    """
}

// Cleaning work - unifying the Controlled Vocabulary
process postprocess {
  //conda "$TOOL_FOLDER/gnps_ml_processing_env2/"
  conda "$TOOL_FOLDER/gnps_ml_processing_env2"

  publishDir "$params.output_dir", mode: 'copy'

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
  if (params.include_massbank)
    """
    python3 $TOOL_FOLDER/GNPS2_Postprocessor.py --includes_massbank
    """
  else
    """
    python3 $TOOL_FOLDER/GNPS2_Postprocessor.py
    """
}

// // A serial piece of code to cache the PubChem names for MatchMS Filtering
// // Doing so avoids excessive parallel API calls
// process cache_pubchem_names_for_matchms {
//   publishDir "./bin/matchms", mode: 'copy'
//   conda "$TOOL_FOLDER/gnps_ml_processing_matchms.yml"

//   maxForks 1

//   cache true

//   input:
//   path cleaned_mgf

//   output:
//   path "pubchem_names.csv", includeInputs: true, emit: pubchem_names

//   """
//   python3 $TOOL_FOLDER/matchms/cache_pubchem_names_for_matchms.py --input_mgf_path ${cleaned_mgf} \
//                                                                   --cached_compound_name_annotation_path "$TOOL_FOLDER/matchms/pubchem_names.csv"
//   """
// }

// Splits the output mgf into smaller chunks to parallelize MatchMS Filtering
process split_mgf_for_matchms_filtering {
  conda "$TOOL_FOLDER/gnps_ml_processing_matchms.yml"

  cache true

  input:
  path cleaned_mgf

  output: 
  path "mgf_chunks/*", emit: mgf_chunks

  // The below script has an optional argument --splitsize, which defaults to 1000. This can be changed to increase or decrease the number of mgf files
  """
  python3 $TOOL_FOLDER/matchms/split_mgf_for_matchms_filtering.py --input_mgf_path ${cleaned_mgf} \
                                                                  --output_path "./mgf_chunks"
  """
}

/* This process reads the csv file, and appends all CCMS IDs to the pubchem name searching cache file to signficantly
reduce the number of api calls */
process spoof_matchms_caching {
  // Curerntly deprecated, see TODO in workflow
  conda "$TOOL_FOLDER/gnps_ml_processing_matchms.yml"
  publishDir "$TOOL_FOLDER/matchms", mode: 'copy', pattern: "compound_name_annotation.csv", saveAs: { filename -> "pubchem_names.csv" } // The script will create this file and copy it back

  cache 'lenient'

  input:
  path cleaned_csv

  output:
  path "compound_name_annotation.csv"
  path "dummy.txt", emit: cache_dummy

  """
  python3 $TOOL_FOLDER/matchms/spoof_matchms_caching.py --input_csv_path ${cleaned_csv} \
                                                        --cached_compound_name_annotation_path "$TOOL_FOLDER/matchms/pubchem_names.csv"
  
  touch dummy.txt
  """
}

// Incoperate MatchMS Filtering into the Pipeline
process matchms_filtering {
  conda "$TOOL_FOLDER/gnps_ml_processing_matchms.yml"

  publishDir "$params.output_dir/matchms_output", mode: 'copy'
  publishDir "$TOOL_FOLDER/matchms", mode: 'copy', pattern: "compound_name_annotation.csv", saveAs: { filename -> "pubchem_names.csv" } // The script will create this file and copy it back

  cache false

  input:
  each cleaned_mgf_chunk
  // file cache_dummy // See TODO: in workflow
  
  output:
  path "matchms_output/*"
  path "compound_name_annotation.csv"

  """
  python3 $TOOL_FOLDER/matchms/matchms_cleaning.py  --input_mgf_path ${cleaned_mgf_chunk}\
                                                    --cached_compound_name_annotation_path "$TOOL_FOLDER/matchms/pubchem_names.csv" \
                                                    --output_path "./matchms_output/" \
  """
}

// Exports the output in JSON format
process export_full_json {
  conda "$TOOL_FOLDER/gnps_ml_processing_env2"

  publishDir "$params.output_dir", mode: 'copy'

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
                                            --progress
  """
}

process generate_subset {
  conda "$TOOL_FOLDER/gnps_ml_processing_env2/"
 
  publishDir "$params.output_dir", mode: 'copy'

  cache false

  input:
  path cleaned_csv
  path cleaned_parquet
  path cleaned_mgf
  // val json_dummy          // Dummy input to force this process to run after the export_full_json process

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
  publishDir "$params.output_dir", mode: 'copy'
  
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
  conda "$TOOL_FOLDER/gnps_ml_processing_env2"

  publishDir "$params.output_dir", mode: 'copy'
  
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

  publishDir "$params.output_dir", mode: 'copy'

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
  environment_creation()
  export_params = prep_params(environment_creation.out.dummy)
  temp_files = export(export_params)

  if (params.include_massbank) {
    fetch_data_massbank()
    prep_params_massbank(fetch_data_massbank.out.blacklist, fetch_data_massbank.out.massbank_data)
    massbank_files = export_massbank(fetch_data_massbank.out.massbank_data, prep_params_massbank.out.params)
    merge_export_massbank(massbank_files.collect())
    temp_files = temp_files.concat(merge_export_massbank.out.merged_files)
  }

  (merged_mgf, merged_csv) = merge_export(temp_files.collect())

  // A python dictionary that maps GNPS adducts to a unified set of adducts used in the GNPS2 workflow
  adduct_mapping_ch = channel.fromPath("$TOOL_FOLDER/adduct_mapping.txt")

  postprocess(merged_csv, merged_mgf, adduct_mapping_ch)
  // cache_pubchem_names_for_matchms(postprocess.out.cleaned_mgf)
  // split_mgf_for_matchms_filtering(postprocess.out.cleaned_mgf)

  /******* 
  * Since every time the cache is used matchms reopens the csv, it's actually faster to just 
  * let the API calls fail, fixing this is TODO
  * spoof_matchms_caching(postprocess.out.cleaned_csv) 
  * matchms_filtering(postprocess.out.cleaned_mgf, spoof_matchms_caching.out.cache_dummy)  
  *******/
  matchms_filtering(postprocess.out.cleaned_mgf)



  /*********** FROM HERE DOWN IS THE ML SPLITS ************/
  if (false) {
    // export_full_json(postprocess.out.cleaned_csv, postprocess.out.cleaned_mgf)
  generate_subset(postprocess.out.cleaned_csv, postprocess.out.cleaned_parquet, postprocess.out.cleaned_mgf)    //, export_full_json.out.dummy
  // if (false) {
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
}