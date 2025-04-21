## Repository Contents and how to run:
This repository contains four seperate workflows:

1) GNPS_ML_Processing: 
A metadata-harmonization pipe for the GNPS spectral libraries that optionally incorperates MassBank.

2) MassBank_Processing:
A workflow to download MassBank and process it into the format needed for harmonization.

3) Train_Test_Splits:
A workflow that generates train, validation, and test sets as specified in the Manuscript.

4) Prebatch:
A workflow to create training batches from datasets.

## Running the Workflows

**Dependencies**
* Python
* Conda (preferably Mamba)
* Nextflow

1) GNPS_ML_Processing
Part of the input for this workflow is the provenance spectra files that back the GNPS libraries.
For your convenience, this pipeline is run daily, and it's output is posted [here](https://external.gnps2.org/gnpslibrary) (see Preprocessed Data).
A dated Nextflow report for the most recent run of the pipeline is available [here](https://external.gnps2.org/admin/download_cleaning_report).

2) MassBank Processing
This workflow is run as a part of (1).

3) Train_Test_Splits
```
cd Train_Test_Splits
nextflow run Dataset_Splitting.nf
```
A full set of parameters and further details can be found in the readme inside the directory.

5) Prebatch:
This workflow is run as a part of (3).
