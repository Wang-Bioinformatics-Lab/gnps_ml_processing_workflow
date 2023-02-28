#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// params.subset = "Bruker_Fragmentation_Prediction"
params.subset = "MH_MNA_Translation"

params.parallelism = 10
// Workflow Boiler Plate
params.OMETALINKING_YAML = "flow_filelinking.yaml"
params.OMETAPARAM_YAML = "job_parameters.yaml"

TOOL_FOLDER = "$baseDir/bin"

process processData {
    publishDir "./nf_output", mode: 'copy'
    publishDir "./GNPS_ml_exports", mode: 'copy'

    conda "$TOOL_FOLDER/conda_env.yml"

    input:
    path './GNPS_ml_exports'

    output:
    path './nf_output'
    path './GNPS_ml_exports'

    """
    python3 $TOOL_FOLDER/GNPS2_Processor.py -p "$params.parallelism"
    python3 $TOOL_FOLDER/GNPS2_Postprocessor.py 
    python3 $TOOL_FOLDER/GNPS2_Subset_Generator.py "$params.subset"
    """
}

workflow {
  data = channel.fromPath('./GNPS_ml_exports')
  processData(data)
}