#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.split  = false
params.matchms_pipeline = true //false
params.export_json = true //false

// Whether to subset the data for ML and perform basic cleaning
params.select_clean_subset = true
params.ion_mode = "positive"  // Which ion mode to select during cleaning
// Which subset ot select, if select_clean_subset is true
// "GNPS_default" for all, other options include "Orbitrap_Fragmentation_Prediction", "MH_MNA_Translation"
params.subset = "Structural_Similarity_Prediction" // None
// params.subset = "Structural_Similarity_Prediction" 

params.output_dir = "./nf_output"
params.GNPS_json_path = "None"

use_default_path = true

// How many chunks to split input files into for data extraction
params.spectra_parallelism = 10000

params.path_to_provenance = "/home/user/LabData/data/GNPS_Library_Provenance/GNPS_Library_Provenance/"

// API Cache path
params.api_cache = ""

params.path_to_nist = ""

// If true, will download and reparse massbank, additionally removing all massbank entires from GNPS
params.include_massbank = true
params.include_riken = false
params.include_mona = true

// Workflow Boiler Plate
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

TOOL_FOLDER = "$baseDir/bin"
params.parallelism = 12
params.pure_networking_parallelism = 5000
params.pure_networking_forks = 32

// Include an option to change the conda path
params.conda_path = "$TOOL_FOLDER/gnps_ml_processing_env2/"
params.matchms_conda_path = "$TOOL_FOLDER/gnps_ml_processing_matchms_env/"

// Include MassBank parser in case params.include_massbank is true
include { fetch_data_massbank; prep_params_massbank; export_massbank; merge_export_massbank } from '../MassBank_Processing/MassBank_processing.nf'

// This environment won't build if called by nextflow, but works fine here
process environment_creation {
  
  output: 
  path "dummy.text", emit: dummy
  
  """
  # If it doesn't exist, create the conda environment
  if [ ! -d "$params.conda_path" ]; then
    mamba env create --prefix $params.conda_path --file $TOOL_FOLDER/conda_env.yml
  else
    echo "Conda environment already exists at $params.conda_path"
  fi

  if [ ! -d "$params.matchms_conda_path" ]; then
    mamba env create --prefix $params.matchms_conda_path --file $TOOL_FOLDER/gnps_ml_processing_matchms.yml
  else
    echo "MatchMS conda environment already exists at $params.matchms_conda_path"
  fi
  # Update the conda environment with the latest packages
  mamba env update --file $TOOL_FOLDER/conda_env.yml --prefix $params.conda_path 
  mamba env update --file $TOOL_FOLDER/gnps_ml_processing_matchms.yml --prefix $params.matchms_conda_path

  touch dummy.text
  """
}

// Splitting out the GNPS Libraries into smaller chunks
process prep_params {
  conda "$params.conda_path"

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
    conda "$params.conda_path"

    maxForks 1
    errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
    maxRetries 5

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

// Incorporate Spectra From Riken Libaries
process export_riken {
  conda "$params.conda_path"

  maxForks 8
  // errorStrategy { sleep(Math.pow(2, task.attempt) * 200 as long); return 'retry' }
  // maxRetries 5

  cache true

  input:
  each download_url

  output:
  path 'temp/*'

  """
  wget "$download_url"

  filename=\$(basename $download_url)
  filename_no_ext="\${filename%.*}"

  mkdir -p "./temp/"

  python3 $TOOL_FOLDER/MSP_ingest.py \
          --msp_path "\$filename" \
          --summary_output_path "./temp/RIKEN-\$filename_no_ext.csv" \
          --spectra_output_path "./temp/RIKEN-\$filename_no_ext.mgf"
  """
}

process export_mona {
  conda "$params.conda_path"

  maxForks 8

  cache true

  input:
  each download_url
  val dummy

  output:
  path 'temp/*'

  """
  wget -qO- "$download_url" | funzip > "MONA_ML_Export.msp"

  mkdir -p "./temp/"

  python3 $TOOL_FOLDER/MSP_ingest.py \
          --msp_path "MONA_ML_Export.msp" \
          --summary_output_path "./temp/MONA-Export.csv" \
          --spectra_output_path "./temp/MONA-Export.mgf"
  """
}

process export_nist {
  conda "$params.conda_path"

  cache true

  input:
  each mgf_file

  output:
  path 'temp/*'

  """

  filename=\$(basename $mgf_file)
  filename_no_ext="\${filename%.*}"
  mkdir -p "./temp/"
  
  python3 $TOOL_FOLDER/MGF_ingest.py \
          --mgf_path "$mgf_file" \
          --summary_output_path "./temp/MGF-Export-\$filename_no_ext.csv" \
          --spectra_output_path "./temp/MGF-Export-\$filename_no_ext.mgf"
  """
  
}

// Merges all the exports together
process merge_export {
  //conda "$params.conda_path/"
  conda "$params.conda_path"

  input:
  path temp_files

  output:
  path './ALL_GNPS_merged.mgf'
  path './ALL_GNPS_merged.csv'

  """
  if [ "$params.include_riken" == "true" ]; then
    riken_flag="--include_riken"
  else
    riken_flag=""
  fi

  if [ "$params.include_massbank" == "true" ]; then
    massbank_flag="--include_massbank"
  else
    massbank_flag=""
  fi

  if [ -n "$params.include_mona" ]; then
    mona_flag="--include_mona"
  else
    mona_flag=""
  fi

  if [ -n "$params.path_to_nist" ]; then
    mgf_export_flag="--include_mgf_exports"
  else
    mgf_export_flag=""
  fi

  # Run the python script with the constructed flags
  python3 $TOOL_FOLDER/merge_files.py \$riken_flag \$massbank_flag \$mgf_export_flag \$mona_flag

  """
}

// Cleaning work - unifying the Controlled Vocabulary
process postprocess {
  //conda "$params.conda_path/"
  conda "$params.conda_path"

  publishDir "$params.output_dir", mode: 'copy'
  publishDir "$TOOL_FOLDER", mode: 'copy', pattern: "smiles_mapping_cache.json", saveAs: { filename -> "smiles_mapping_cache.json" } // The script will create this file and copy it back

  cache true

  input:
  path merged_csv
  path merged_parquet
  path adduct_mapping

  output:
  path "ALL_GNPS_cleaned.csv", emit: cleaned_csv
  path "ALL_GNPS_cleaned.parquet", emit: cleaned_parquet, optional: true
  path "ALL_GNPS_cleaned.mgf", emit: cleaned_mgf
  path "ALL_GNPS_cleaned.json", emit: cleaned_json
  path "smiles_mapping_cache.json", includeInputs: true

  script:
  if (params.include_massbank)
    """
    if $params.include_riken; then
      riken_flag="--includes_riken"
    else
      riken_flag=""
    fi

    if [ -n "$params.include_mona" ]; then
      mona_flag="--includes_mona"
    else
      mona_flag=""
    fi

    cp $TOOL_FOLDER/smiles_mapping_cache.json ./smiles_mapping_cache.json
    python3 $TOOL_FOLDER/GNPS2_Postprocessor.py --includes_massbank --smiles_mapping_cache "smiles_mapping_cache.json" \$riken_flag \$mona_flag
    """
  else
    """
    cp $TOOL_FOLDER/smiles_mapping_cache.json ./smiles_mapping_cache.json
    python3 $TOOL_FOLDER/GNPS2_Postprocessor.py --smiles_mapping_cache "smiles_mapping_cache.json"
    """
}

process integrate_api_info {
  conda "$params.conda_path"

  publishDir "$params.output_dir", mode: 'copy'

  cache false

  input:
  path cleaned_csv
  path chem_info_service_api_cache, stageAs: 'api_cache/ChemInfoService/*'  // TODO: Can we symlink the directory rather than the files to save time?
  path classyfire_api_cache, stageAs: 'api_cache/Classyfire/*'
  path npclassifier_api_cache, stageAs: 'api_cache/Npclassifier/*'

  output:
  path "ALL_GNPS_cleaned_enriched.csv", emit: cleaned_csv
  
  """
  python3 $TOOL_FOLDER/integrate_API_info.py   --input_csv_path ${cleaned_csv} \
                                               --output_csv_path "./ALL_GNPS_cleaned_enriched.csv" \
                                               --api_cache_path "./api_cache"
  """
}

// Incoperate MatchMS Filtering into the Pipeline
process matchms_filtering {
  conda "$params.matchms_conda_path"

  publishDir "$params.output_dir", mode: 'copy'
  publishDir "$TOOL_FOLDER/matchms", mode: 'copy', pattern: "compound_name_annotation.csv", saveAs: { filename -> "pubchem_names.csv" } // The script will create this file and copy it back

  cache true

  input:
  each cleaned_mgf_chunk
  // file cache_dummy // See TODO: in workflow
  
  output:
  path "matchms_output/*", optional: true
  path "compound_name_annotation.csv", optional: true

  """
  mkdir -p ./matchms_output/

  python3 $TOOL_FOLDER/matchms/matchms_cleaning.py  --input_mgf_path ${cleaned_mgf_chunk}\
                                                    --cached_compound_name_annotation_path "$TOOL_FOLDER/matchms/pubchem_names.csv" \
                                                    --output_path "./matchms_output/" \
  """
}

// Exports the output in JSON format
process export_full_json {
  conda "$params.conda_path"

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
  conda "$params.conda_path"
 
  publishDir "$params.output_dir", mode: 'copy'

  cache false

  input:
  path cleaned_csv
  // path cleaned_parquet
  path cleaned_mgf
  // val json_dummy          // Dummy input to force this process to run after the export_full_json process

  output:
  path "summary/*"
  path "spectra/*.parquet", emit: output_parquet, optional: true
  path "spectra/*.mgf", emit: output_mgf
  path "json_outputs/*.json", emit: output_json, optional: true
  path "util/*", optional: true

  """
  python3 $TOOL_FOLDER/GNPS2_Subset_Generator.py "$params.subset"
  """
}

process generate_mgf {
  conda "$params.conda_path"

  input:
  each parquet_file

  output:
  path "*.mgf", emit: output_mgf

  """
  python3 $TOOL_FOLDER/GNPS2_MGF_Generator.py -input_path "$parquet_file"
  """
}

process select_data_for_ml {
  conda "$params.conda_path"
  publishDir "${params.output_dir}/ML_ready_subset_${params.ion_mode}/", mode: 'copy'

  input:
  path csv_file
  path mgf_file

  cache false

  output:
  path 'selected_summary.csv',  emit: selected_summary
  path 'selected_spectra.mgf',  emit: selected_spectra
  path 'selected_spectra.json',  emit: selected_spectra_json

  """
  python3 $TOOL_FOLDER/select_data.py \
          --input_csv "$csv_file" \
          --input_mgf "$mgf_file" \
          --ion_mode $params.ion_mode
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
  //conda "$params.conda_path/"
  conda "$params.conda_path"

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

  if (params.include_riken) {
    riken_urls = Channel.of("http://prime.psc.riken.jp/compms/msdial/download/msp/MSMS-Neg-RikenOxPLs.msp",
                            "http://prime.psc.riken.jp/compms/msdial/download/msp/MSMS-Pos-PlaSMA.msp",
                            "http://prime.psc.riken.jp/compms/msdial/download/msp/MSMS-Neg-PlaSMA.msp",
                            "http://prime.psc.riken.jp/compms/msdial/download/msp/BioMSMS-Pos-PlaSMA.msp",
                            "http://prime.psc.riken.jp/compms/msdial/download/msp/BioMSMS-Neg-PlaSMA.msp",
                            "http://prime.psc.riken.jp/compms/msdial/download/msp/MSMS-Neg-PFAS_20200806.msp",
                            "http://prime.psc.riken.jp/compms/msdial/download/msp/MSMS-Pos-bmdms-np_20200811.msp"
                            )
    riken_files = export_riken(riken_urls)
    temp_files = temp_files.concat(riken_files)

  }

  if (params.include_mona) {
    mona_urls = Channel.of("https://mona.fiehnlab.ucdavis.edu/rest/downloads/retrieve/4ebf9d78-f0d4-470c-b8cb-afc5f6241584")
    mona_files = export_mona(mona_urls, environment_creation.out.dummy)
    temp_files = temp_files.concat(mona_files)
  }

  if (params.path_to_nist != "") {
    nist_ch = Channel.fromPath(params.path_to_nist)
    nist_files = export_nist(nist_ch)
    temp_files = temp_files.concat(nist_files)
  }
    
  (merged_mgf, merged_csv) = merge_export(temp_files.collect())

  

  // A python dictionary that maps GNPS adducts to a unified set of adducts used in the GNPS2 workflow
  adduct_mapping_ch = channel.fromPath("$TOOL_FOLDER/adduct_mapping.txt")

  postprocess(merged_csv, merged_mgf, adduct_mapping_ch)

  if (params.api_cache != "") {
    chem_info_service_api_cache = Channel.fromPath(params.api_cache + "/ChemInfoService/**/*.json")
    classyfire_api_cache        = Channel.fromPath(params.api_cache + "/Classyfire/**/*.json")
    npclassifier_api_cache      = Channel.fromPath(params.api_cache + "/Npclassifier/**/*.json")
    api_info_csv = integrate_api_info(  postprocess.out.cleaned_csv, 
                                        chem_info_service_api_cache.collect(),
                                        classyfire_api_cache.collect(),
                                        npclassifier_api_cache.collect()
                                      )
  } else{
    api_info_csv = postprocess.out.cleaned_csv
  }

  if (params.export_json) {
    export_full_json(api_info_csv, postprocess.out.cleaned_mgf)
  }

  if (params.matchms_pipeline) {
    matchms_filtering(postprocess.out.cleaned_mgf)
  }

  if ((params.select_clean_subset == false) && (params.subset != "None")) {
    // Warn user that subset is not being generated
    println "Subset generation is disabled, skipping subset generation"
  }

  if (params.select_clean_subset) {
    // Perform basic data selection and cleaning for ML
    select_data_for_ml(api_info_csv, postprocess.out.cleaned_mgf)

  }
}