#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Workflow Boiler Plate
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

TOOL_FOLDER_LS = "$moduleDir/bin"
FAST_SEARCH_LIBRARY_BIN = "$TOOL_FOLDER_LS/GNPS_FastSearch_Library/bin"

// Input data
params.input_csv = "https://fileserver.wanglab.science/p_ml_cleanup/to_upload/cleaning_outputs/ML_ready_subset_positive/selected_summary.csv"
params.input_mgf = "https://fileserver.wanglab.science/p_ml_cleanup/to_upload/cleaning_outputs/ML_ready_subset_positive/selected_spectra.mgf"

// Which task subset to use
params.subset = "Structural_Similarity_Prediction"
params.split_type = "sample_structures_smart_inchikey"//"structure_smart"  // 'basic_sampling_scheme', 'structure_smart', 'random' (random spectra), or 'structure' (random inchi14)

params.test_set_num = 500//2300 // random structure=16154  // An integer (or float) representing the number (or percentage of) data points to use as a test set

// Batch Global Generation Parameters
params.num_epochs = 150
params.split_size = 10  // How many epochs each worker has, num_epochs should be divisible by this number
params.batch_size = 32
params.num_turns = 2    // Number of times each inchikey shows up per epoch
params.save_dir = "./nf_output/${params.subset}/${params.split_type}/"

// Which mass analyzers to include in filtered data.
params.mass_analyzer_lst = "" // Default is "", which is all mass analyzers. Otherwise, a semicolon delimited list of mass analyzers to include

// Prebatching Parameters
params.prebatch = true
params.unfiltered = true
params.filtered = true
params.include_biased = false
params.include_strict_ce = false

// Generate test sets
params.generate_test_set = true

// Parallelism Settings
params.parallelism = 5

// Submodule for batch generation
include { sampleSeveralEpochs as unfiltered_train_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as unfiltered_train_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as unfiltered_val_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as unfiltered_val_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as filtered_train_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as filtered_train_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as filtered_val_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as filtered_val_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as unfiltered_biased_train_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as unfiltered_biased_train_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as unfiltered_biased_val_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as unfiltered_biased_val_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as filtered_biased_train_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as filtered_biased_train_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as filtered_biased_val_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as filtered_biased_val_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as filtered_strict_ce_train_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as filtered_strict_ce_train_assembler } from '../Prebatch/presample_and_assemble_data.nf'
include { sampleSeveralEpochs as filtered_strict_ce_val_sampler } from '../Prebatch/presample_and_assemble_data.nf'
include { assembleEpochs      as filtered_strict_ce_val_assembler } from '../Prebatch/presample_and_assemble_data.nf'


include { exhaustivePairEnumeration as unfiltered_test_enumerator } from '../Prebatch/presample_and_assemble_data.nf'
include { exhaustivePairEnumeration as filtered_test_enumerator } from '../Prebatch/presample_and_assemble_data.nf'
include { exhaustivePairEnumeration as filtered_strict_ce_test_enumerator } from '../Prebatch/presample_and_assemble_data.nf'


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

process calculate_pairwise_similarity {
  conda "$TOOL_FOLDER_LS/conda_env.yml"
  // publishDir "./nf_output", mode: 'copy'

  cache true

  input:
  path metadata_csv

  output:
  path 'pairwise_similarities.csv', emit: pairwise_similarities

  """
  python3 $TOOL_FOLDER_LS/calculate_pairwise_similarity.py --input_metadata $metadata_csv --output_filename pairwise_similarities.csv
  """
}

/*
 * Generates a random subset of params.test_set_num  data points to use as a test set
 * These data points should be consistent across splitting methodology
 */
process generate_test_set {
  conda "$TOOL_FOLDER_LS/conda_env.yml"
  publishDir "./nf_output/${params.subset}/${params.split_type}", mode: 'copy'

  cache true

  input:
  path csv_file
  path pairwise_similarities
  path mgf_file

  output:
  path 'test_rows.csv',  emit: test_rows_csv
  path 'test_rows.mgf',  emit: test_rows_mgf
  path 'test_rows.json', emit: test_rows_json
  path 'test_similarities.csv', emit: test_similarities_csv
  path 'train_rows.csv', emit: train_rows_csv
  path 'train_rows.mgf', emit: train_rows_mgf
  path 'train_rows.json', emit: train_rows_json
  path 'train_similarities.csv', emit: train_similarities_csv
  path 'val_rows.csv', emit: val_rows_csv
  path 'val_rows.mgf', emit: val_rows_mgf
  path 'val_rows.json', emit: val_rows_json
  path 'val_similarities.csv', emit: val_similarities_csv
  path 'train_test_similarities.csv', emit: train_test_similarities_csv

  """
  python3 $TOOL_FOLDER_LS/calc_test.py \
          --input_csv "$csv_file" \
          --input_similarities "$pairwise_similarities" \
          --input_mgf "$mgf_file" \
          --num_test_points $params.test_set_num \
          --sampling_strategy "$params.split_type" \
          --threshold 0.5 \
          --debug
  """
}

workflow {
  csv_file  = Channel.fromPath(params.input_csv)
  mgf_file  = Channel.fromPath(params.input_mgf)

  generate_subset(csv_file, mgf_file)
  calculate_pairwise_similarity(generate_subset.out.output_summary)
  generate_test_set(generate_subset.out.output_summary, calculate_pairwise_similarity.out.pairwise_similarities, generate_subset.out.output_mgf)  

  // To split the epochs
  epoch_ch = Channel.of(1..params.num_epochs/params.split_size)
  val_epoch_ch = Channel.of(1,)

  // Unfiltered
  if (params.prebatch && params.unfiltered) {
    unfiltered_train_sampler(
                        epoch_ch,
                        generate_test_set.out.train_rows_csv, 
                        generate_test_set.out.train_similarities_csv, 
                        tuple(  params.batch_size,
                                params.num_turns,
                                "standard",  // mode
                                "",  // mass_analyzer_lst
                                10,  // num_bins
                                "", // exponential_bins
                                false, // strict_collision_energy
                                false, // force_one_epoch
                                "train_unfiltered.hdf5", // output_name
                            ))
    unfiltered_train_assembler(unfiltered_train_sampler.out.prebatched_files.collect(), 
                              unfiltered_train_sampler.out.name_file.take(1))

    unfiltered_val_sampler(
                            val_epoch_ch,
                            generate_test_set.out.val_rows_csv, 
                            generate_test_set.out.val_similarities_csv, 
                            tuple(  params.batch_size,
                                    params.num_turns,
                                    "standard",  // mode
                                    "",   // mass_analyzer_lst
                                    10,  // num_bins
                                    "", // exponential_bins
                                    false, // strict_collision_energy
                                    true, // force_one_epoch
                                    "val_unfiltered.hdf5", // output_name
                                ))
    unfiltered_val_assembler(unfiltered_val_sampler.out.prebatched_files.collect(), 
                              unfiltered_val_sampler.out.name_file.take(1))

  // Unfiltered Biased
    if (params.include_biased) {
      unfiltered_biased_train_sampler(
                          epoch_ch,
                          generate_test_set.out.train_rows_csv, 
                          generate_test_set.out.train_similarities_csv, 
                          tuple(  params.batch_size,
                                  params.num_turns,
                                  "standard",  // filter
                                  "",  // mass_analyzer_lst
                                  20,  // num_bins
                                  0.3, // exponential_bins
                                  false, // strict_collision_energy
                                  false, // force_one_epoch
                                  "train_unfiltered_biased.hdf5", // output_name
                              ))
      unfiltered_biased_train_assembler(unfiltered_biased_train_sampler.out.prebatched_files.collect(), 
                                unfiltered_biased_train_sampler.out.name_file.take(1))

      unfiltered_biased_val_sampler(
                          val_epoch_ch,
                          generate_test_set.out.val_rows_csv,
                          generate_test_set.out.val_similarities_csv, 
                          tuple(  params.batch_size,
                                  params.num_turns,
                                  "standard",  // mode
                                  "",  // mass_analyzer_lst
                                  20,  // num_bins
                                  0.3, // exponential_bins
                                  false, // strict_collision_energy
                                  true, // force_one_epoch
                                  "val_unfiltered_biased.hdf5", // output_name
                              ))
      unfiltered_biased_val_assembler(unfiltered_biased_val_sampler.out.prebatched_files.collect(),
                              unfiltered_biased_val_sampler.out.name_file.take(1))
    }
  }

  // Filtered
  if (params.prebatch && params.filtered) {
    filtered_train_sampler(
                        epoch_ch,
                        generate_test_set.out.train_rows_csv, 
                        generate_test_set.out.train_similarities_csv, 
                        tuple(  params.batch_size,
                                params.num_turns,
                                "filter",  // mode
                                params.mass_analyzer_lst,
                                10,  // num_bins
                                "", // exponential_bins
                                false, // strict_collision_energy
                                false, // force_one_epoch
                                "train_filtered.hdf5", // output_name
                            ))
    filtered_train_assembler(filtered_train_sampler.out.prebatched_files.collect(), 
                              filtered_train_sampler.out.name_file.take(1))

    filtered_val_sampler(
                        val_epoch_ch,
                        generate_test_set.out.val_rows_csv,
                        generate_test_set.out.val_similarities_csv, 
                        tuple(  params.batch_size,
                                params.num_turns,
                                "filter",  // mode
                                params.mass_analyzer_lst,
                                10,  // num_bins
                                "", // exponential_bins
                                false, // strict_collision_energy
                                true, // force_one_epoch
                                "val_filtered.hdf5", // output_name
                            ))
    filtered_val_assembler(filtered_val_sampler.out.prebatched_files.collect(),
                            filtered_val_sampler.out.name_file.take(1))
    // Filtered Biased
    if (params.include_biased) {
      filtered_biased_train_sampler(
                      epoch_ch,
                      generate_test_set.out.train_rows_csv, 
                      generate_test_set.out.train_similarities_csv, 
                      tuple(  params.batch_size,
                              params.num_turns,
                              "filter",  // mode
                              params.mass_analyzer_lst,
                              20,  // num_bins
                              0.3, // exponential_bins
                              false, // strict_collision_energy
                              false, // force_one_epoch
                              "train_filtered_biased.hdf5", // output_name
                          ))
      filtered_biased_train_assembler(filtered_biased_train_sampler.out.prebatched_files.collect(), 
                                filtered_biased_train_sampler.out.name_file.take(1))  

      filtered_biased_val_sampler(
                          val_epoch_ch,
                          generate_test_set.out.val_rows_csv,
                          generate_test_set.out.val_similarities_csv, 
                          tuple(  params.batch_size,
                                  params.num_turns,
                                  "filter",  // mode
                                  params.mass_analyzer_lst,
                                  20,  // num_bins
                                  0.3, // exponential_bins
                                  false, // strict_collision_energy
                                  true, // force_one_epoch
                                  "val_filtered_biased.hdf5", // output_name
                              ))

      filtered_biased_val_assembler(filtered_biased_val_sampler.out.prebatched_files.collect(),
                            filtered_biased_val_sampler.out.name_file.take(1))

    }

    // Filtered Strict CE
    if (params.include_strict_ce) {
      filtered_strict_ce_train_sampler(
                      epoch_ch,
                      generate_test_set.out.train_rows_csv, 
                      generate_test_set.out.train_similarities_csv, 
                      tuple(  params.batch_size,
                              params.num_turns,
                              "filter",  // mode
                              params.mass_analyzer_lst,
                              10,  // num_bins
                              "", // exponential_bins
                              true, // strict_collision_energy
                              false, // force_one_epoch
                              "train_filtered_strict_collision_energy.hdf5", // output_name
                          ))
      filtered_strict_ce_train_assembler(filtered_strict_ce_train_sampler.out.prebatched_files.collect(), 
                            filtered_strict_ce_train_sampler.out.name_file.take(1))

      filtered_strict_ce_val_sampler(
                      val_epoch_ch,
                      generate_test_set.out.val_rows_csv,
                      generate_test_set.out.val_similarities_csv, 
                      tuple(  params.batch_size,
                              params.num_turns,
                              "filter",  // mode
                              params.mass_analyzer_lst,
                              10,  // num_bins
                              "", // exponential_bins
                              true, // strict_collision_energy
                              true, // force_one_epoch
                              "val_filtered_strict_collision_energy.hdf5", // output_name
                          ))

      filtered_strict_ce_val_assembler(filtered_strict_ce_val_sampler.out.prebatched_files.collect(),
                          filtered_strict_ce_val_sampler.out.name_file.take(1))
    }
  }

  if (params.generate_test_set) {
    // Unfiltered Test
    unfiltered_test_enumerator(generate_test_set.out.test_rows_csv, 
                              generate_test_set.out.test_similarities_csv,
                              generate_test_set.out.train_test_similarities_csv,
                              tuple( false, // mode
                                      "", // mass_analyzer_lst
                                      false, // strict_collision_energy
                                      "unfiltered_test.parquet"
                              )
                              )

    // Filtered Test
    filtered_test_enumerator(generate_test_set.out.test_rows_csv, 
                            generate_test_set.out.test_similarities_csv,
                            generate_test_set.out.train_test_similarities_csv,
                            tuple( true, // mode
                                    params.mass_analyzer_lst, // mass_analyzer_lst
                                    false, // strict_collision_energy
                                    "filtered_test.parquet"
                            )
                            )

    if (params.include_strict_ce) {
      // Filtered Test (Strict CE)
      filtered_strict_ce_test_enumerator(generate_test_set.out.test_rows_csv, 
                              generate_test_set.out.test_similarities_csv,
                              generate_test_set.out.train_test_similarities_csv,
                              tuple( true, // mode
                                      params.mass_analyzer_lst, // mass_analyzer_lst
                                      true, // strict_collision_energy
                                      "filtered_strict_collision_energy_test.parquet"
                              )
                              )
    }
  }
}