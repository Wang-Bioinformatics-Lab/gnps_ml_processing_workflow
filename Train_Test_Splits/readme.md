## Dataset Splitting
This workflow splits the source datasetss into training, validation, and test sets.

In addition, it partitions the training and validation sets into epochs of training data stored in an hdf format. The number of epochs can be specified with the `--num_epochs` argument. Batch size can be spcified with the `--batch_size` argument. 

### Requirements
* conda (Mamba Preferred)
* Nextflow

### Running the workflow
```
nextflow run Dataset_Splitting.nf
```

### Parameters
**Generating Training Data**
* `--prebatch <value>`: Enables training data generation, `true` by default
* `--unfiltered <value>`: Enables generation of the unfiltered set, `true` by default
* `--filtered <value>`: Enables generation of the filtered set, `true` by default
* `--include_biased <value>`: Enables generation of a training set with a skewed pairwise similarity distribution. `false` by default
* `--num_epochs <value>`: Number of pre-batched epochs to generate, `150` by default. Must be divisible by `split_size`
* `--batch_size <value>`: Number of values per batch, `32` by default
* `--split_size <value>`: Number of epochs for each worker to generate, `10` by default

**Generating Testing Data**
* `--generate_test_set <value>`: Enables generation of test sets, `true` by default

**Misc**
* `--include_strict_ce <value>`: Enables the generation of the strict collision energy sets. Applies to both training and test, `false` by default.
* `--parallelism <value>`: Max number of parallel workers.