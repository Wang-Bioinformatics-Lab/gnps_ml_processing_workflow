#!/usr/bin/env nextflow
nextflow.enable.dsl=2

TOOL_FOLDER = "$moduleDir/bin"

params.parallelism = 5
params.metadata = ''
params.tanimoto_scores_path = ''
params.batch_size = 32
params.num_turns = ''
params.num_bins = ''
params.exponential_bins = ''
params.strict_collision_energy = false
params.mode = 'standard' // choices=['standard', 'filter', 'triplet']

params.save_dir = './nf_output'
params.temp_file_dir = '/home/user/LabData/michael_s/temp_work_dir/'
params.output_name = 'my_data.hdf5'
params.num_epochs = 1800

// Filtering Parameters
params.mass_analyzer_lst = ""   // A semicolon delimited list of mass analyzers to include, an empty string means all are included

params.force_one_epoch = false

split_size = 10     // How many epochs each worker has, num_epochs should be divisible by this number to avoid errors

// This process samples one epoch of data
process sampleSeveralEpochs {
    conda "$TOOL_FOLDER/conda_env.yml"

    // errorStrategy { sleep(Math.pow(2, task.attempt) * 60 as long); return 'retry' }
    // maxRetries 5

    maxForks params.parallelism
    memory '10GB'
    cpus 2

    input:
    each epoch_num 
    path metadata
    path tanimoto_scores_path
    tuple   val(batch_size),
            val(num_turns),
            val(mode),
            val(mass_analyzer_lst),
            val(num_bins),
            val(exponential_bins),
            val(strict_collision_energy),
            val(force_one_epoch),
            val(output_name)

    output:
    path '*.hdf5', emit: prebatched_files
    path 'name.txt', emit: name_file

    """
    # Echo conda path

    # For each parameter, if it is not '', then we assemble it and add it to the command.
    # This will allow nextflow to work as a flexible wrapper around the python script and
    # the python script will handle any parameter errors.
    echo "params.metadata: $metadata"
    if [ $metadata != '' ]; then
        metadata='--metadata $metadata'
    else
        metadata=''
    fi
    echo "tanimoto_scores_path: $tanimoto_scores_path"
    if [ $tanimoto_scores_path != '' ]; then
        tanimoto_scores_path='--tanimoto_scores_path $tanimoto_scores_path'
    else
        tanimoto_scores_path=''
    fi
    if [ $batch_size != '' ]; then
        batch_size='--batch_size $batch_size'
    else
        batch_size=''
    fi
    if [ $num_turns != '' ]; then
        num_turns='--num_turns $num_turns'
    else
        num_turns=''
    fi
    if [ $mass_analyzer_lst != '' ]; then
        mass_analyzer_lst='--mass_analyzer_lst $mass_analyzer_lst'
    else
        mass_analyzer_lst=''
    fi
    if [ $num_bins != '' ]; then
        num_bins='--num_bins $num_bins'
    else
        num_bins=''
    fi
    if [ $exponential_bins != '' ]; then
        exponential_bins='--exponential_bins $exponential_bins'
    else
        exponential_bins=''
    fi
    if [ $strict_collision_energy = true ]; then
        strict_collision_energy='--strict_collision_energy'
    else
        strict_collision_energy=''
    fi

    # For validation, we'll want exactly one epoch
    if [ $force_one_epoch = true ]; then
        num_epochs='--num_epochs 1'
    else
        num_epochs='--num_epochs $split_size'
    fi

    # Print all parameters
    echo "metadata: \$metadata"
    echo "tanimoto_scores_path: \$tanimoto_scores_path"
    echo "batch_size: \$batch_size"
    echo "num_turns: \$num_turns"

    mkdir -p logs

    # Make a file called "name.txt" that output_name, for use later renaming files
    echo $output_name > name.txt

    python $TOOL_FOLDER/presample_pairs_generic.py  \$metadata \
                                                    \$tanimoto_scores_path \
                                                    \$batch_size \
                                                    \$num_turns \
                                                    \$num_epochs \
                                                    --mode $mode \
                                                    \$mass_analyzer_lst \
                                                    \$exponential_bins \
                                                    \$strict_collision_energy \
                                                    \$num_bins \
                                                    --save_dir "./" \
                                                    --seed $epoch_num
    """
}

// This process assembles the sampled epochs into one hdf5 file
process assembleEpochs {
    publishDir "${params.save_dir}/prebatched_data", mode: "copy"

    memory '10GB'
    cpus 2

    conda "$TOOL_FOLDER/conda_env.yml"

    input: 
    path hdf5_file, stageAs: "hdf5_file_*.hdf5"
    path name_file

    output:
    path "*.hdf5"

    """
    # Contents of name_file contain final output name
    output_name=\$(cat $name_file)

    python3 $TOOL_FOLDER/assemble_epochs.py --output_name \$output_name
    """
}

process exhaustivePairEnumeration {
    publishDir "${params.save_dir}/test_pairs", mode: "copy"

    conda "$TOOL_FOLDER/conda_env.yml"

    memory '75GB'
    cpus params.parallelism

    input:
    path metadata
    path pairwise_sims
    path train_test_sims
    tuple   val(filter),
            val(mass_analyzer_lst),
            val(strict_collision_energy),
            val(output_name)

    output:
    path "*.parquet"

    """
    if [ "$filter" = true ]; then
        filter_flag='--filter'
    else
        filter_flag=''
    fi
    if [ "$mass_analyzer_lst" != '' ]; then
        mass_analyzer_lst='--mass_analyzer_lst $mass_analyzer_lst'
    else
        mass_analyzer_lst=''
    fi
    if [ "$strict_collision_energy" = true ]; then
        strict_collision_energy='--strict_collision_energy'
    else
        strict_collision_energy=''
    fi

    python3 $TOOL_FOLDER/exhaustive_generation.py   --n_jobs $params.parallelism \
                                                    --metadata_path $metadata \
                                                    --pairwise_similarities_path "$pairwise_sims" \
                                                    --train_test_similarities_path "$train_test_sims" \
                                                    --output_path "${output_name}" \
                                                    --temp_file_dir "${params.temp_file_dir}" \
                                                    \$filter_flag \
                                                    \$mass_analyzer_lst \
                                                    \$strict_collision_energy

    """
}

workflow {
    epoch_ch = Channel.of(1..params.num_epochs/split_size)

    if (params.metadata != '') {
        metadata_ch = Channel.fromPath(params.metadata)
    } else {
        metadata_ch = Channel.fromPath('NO_FILE')
    }
    if (params.tanimoto_scores_path != '') {
        tanimoto_scores_path_ch = Channel.fromPath(params.tanimoto_scores_path)
    } else {
        tanimoto_scores_path_ch = Channel.fromPath('NO_FILE2')
    }

    parameter_ch = Channel.of(
                            tuple(  params.batch_size,
                                    params.num_turns,
                                    params.mode,
                                    params.mass_analyzer_lst,
                                    params.num_bins,
                                    params.exponential_bins,
                                    params.strict_collision_energy,
                                    params.force_one_epoch,
                                    params.output_name
                            )
    )
    
    sampleSeveralEpochs(epoch_ch, metadata_ch, tanimoto_scores_path_ch, parameter_ch)
    assembleEpochs(sampleSeveralEpochs.out.prebatched_files.collect(), sampleSeveralEpochs.out.name_file.take(1))
}