#!/bin/bash

path="../Train_Test_Splits/nf_output/Structural_Similarity_Prediction/sample_structures_smart_inchikey/"

nextflow run presample_and_assemble_data.nf \
                --parallelism 2 \
                --metadata $path"train_rows.csv" \
                --tanimoto_scores_path $path"train_similarities.csv" \
                --batch_size 32 \
                --num_turns 2 \
                --num_bins 10 \
                --num_epochs 20
