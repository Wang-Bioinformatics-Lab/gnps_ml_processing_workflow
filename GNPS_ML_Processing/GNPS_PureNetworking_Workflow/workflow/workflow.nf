#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Input Data
params.inputspectra = 'data/specs_ms.mgf'

// Similarity Type
params.similarity = "gnps"

params.parallelism = 24
params.maxforks = 4

// GNPS parameters
params.min_matched_peaks = 6
params.ms2_tolerance = 0.5
params.pm_tolerance = 0.5
params.min_cosine = 0.7
params.max_shift = 1000

maxforks_int = params.maxforks.toInteger()

// Boiler plate
TOOL_FOLDER = "$baseDir/bin"
params.publishdir = "nf_output"

process prepGNPSParams {
    publishDir "$params.publishdir", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    file spectrum_file

    output:
    file "params/*"

    """
    mkdir params
    python $TOOL_FOLDER/prep_molecular_networking_parameters.py \
        "$spectrum_file" \
        "params" \
        --parallelism "$params.parallelism" \
        --min_matched_peaks "$params.min_matched_peaks" \
        --ms2_tolerance "$params.ms2_tolerance" \
        --pm_tolerance "$params.pm_tolerance" \
        --min_cosine "$params.min_cosine" \
        --max_shift "$params.max_shift"
    """
}

process calculateGNPSPairs {
    publishDir "$params.publishdir/temp_pairs", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    maxForks maxforks_int

    input:
    file spectrum_file
    each file(params_file)

    output:
    file "*_aligns.tsv" optional true

    """
    $TOOL_FOLDER/main_execmodule \
        ExecMolecularParallelPairs \
        "$params_file" \
        -ccms_INPUT_SPECTRA_MS2 $spectrum_file \
        -ccms_output_aligns ${params_file}_aligns.tsv
    """
}

process blinkQuantizeData {
    publishDir "$params.publishdir/temp_pairs", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    cache false
    echo true

    input:
    file spectrum_file

    output:
    file "*npz"

    """
    python $TOOL_FOLDER/blink/blink.py \
        $spectrum_file
    """
}

process calculatePairsBlink {
    publishDir "$params.publishdir/temp_pairs", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    cache false
    echo true

    input:
    file file_quantized

    output:
    file "*npz" 

    """
    python $TOOL_FOLDER/blink/blink.py \
        $file_quantized
    """
}

// process formatPairsBlink {
//     publishDir "$params.publishdir", mode: 'copy'

//     input:
//     file blink_pairs from _pairsblink_ch

//     output:
//     file "merged_pairs.tsv"

//     """
//     python $TOOL_FOLDER/blink_reformat.py \
//         $blink_pairs \
//         merged_pairs.tsv
//     """
// }

process calculatePairsSimile {
    publishDir "$params.publishdir/temp_pairs", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    cache false
    //echo true

    maxForks maxforks_int

    input:
    file spectrum_file
    each i

    output:
    file "*.tsv" optional true

    """
    python $TOOL_FOLDER/run_simile.py \
        $spectrum_file \
        ${i}_pairs.tsv \
        --nodenumber $i \
        --nodetotal $params.parallelism
    """
}

process calculatePairsEntropy {
    publishDir "$params.publishdir/temp_pairs", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env_entropy.yml"

    cache false
    echo true

    maxForks maxforks_int

    input:
    file spectrum_file
    each i

    output:
    file "*.tsv" optional true

    """
    python $TOOL_FOLDER/run_SpectralEntropy.py \
        $spectrum_file \
        ${i}_pairs.tsv \
        --nodenumber ${i} \
        --nodetotal ${params.parallelism} \
        --min_matched_peaks ${params.min_matched_peaks} \
        --ms2_tolerance ${params.ms2_tolerance} \
        --min_cosine ${params.min_cosine}

    """
}

// TODO: Adding MS2DeepScore

workflow {
    _spectra_ch = Channel.fromPath( params.inputspectra )
    params_top = params.parallelism.toInteger() - 1

    if(params.similarity == "gnps"){
        all_params_ch = prepGNPSParams(_spectra_ch)
        _pairs_ch = calculateGNPSPairs(_spectra_ch, all_params_ch.flatten())
        _pairs_ch.collectFile(name: "merged_pairs.tsv", storeDir: "$params.publishdir", keepHeader: true)
    }
    else if(params.similarity == "blink"){
        _blink_quantized_ch = blinkQuantizeData(_spectra_ch)
        _blink_pairs_ch = calculatePairsBlink(_blink_quantized_ch)
    }
    else if(params.similarity == "simile"){
        parallel_ch = Channel.from(0..params_top)
        _pairs_simile_ch = calculatePairsSimile(_spectra_ch, parallel_ch)
        _pairs_simile_ch.collectFile(name: "merged_pairs.tsv", storeDir: "$params.publishdir", keepHeader: true)
    }
    else if(params.similarity == "entropy"){
        entropy_parallel_ch = Channel.from(0..params_top)
        _pairs_entropy_ch = calculatePairsEntropy(_spectra_ch, entropy_parallel_ch)
        _pairs_entropy_ch.collectFile(name: "merged_pairs.tsv", storeDir: "$params.publishdir", keepHeader: true)
    }
}