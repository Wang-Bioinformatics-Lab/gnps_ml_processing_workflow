# GNPS ML Export Workflow

## Running The Workflow
Building the conda environment
```
mamba env create --file ./bin/conda_env.yml --prefix ./bin/gnps_ml_processing_env2/
```

To run the workflow to test simply do

```
nohup nextflow run nf_workflow.nf --path_to_provenance=<your_path_to_gnps_library_provenance> -bg &
```

(by default path_to_provenance is /home/user/LabData/GNPS_Library_Provenance/)

You can specify a specific subset using one of the following:
['Bruker_Fragmentation_Prediction', 'MH_MNA_Translation', 'Orbitrap_Fragmentation_Prediction','Thermo_Bruker_Translation', 'Structural_Modification', 'Structural_Similarity_Prediction', 'Spectral_Similarity_Prediction']

E.g: ```nohup nextflow run nf_workflow.nf --subset=MH_MNA_Translation -bg &```

## Memory Requirements
Right now, the workflow requires at least 100 GB of memory to run, 16 processors is recomended.

To learn NextFlow checkout this documentation:

https://www.nextflow.io/docs/latest/index.html

## Deployment to GNPS2

In order to deploy, we have a set of deployment tools that will enable deployment to the various gnps systems. To run the deployment, use the following commands from the deploy_gnps2 folder. 

```
make deploy-prod
```

You might need to checkout the module, do this by running

```
git submodule init
git submodule update
```
