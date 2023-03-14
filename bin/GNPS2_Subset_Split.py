import sys
import pandas as pd
from utils import build_tanimoto_similarity_list_precomputed, get_splits_index_only
import vaex

# argv should look like ['/home/user/SourceCode/GNPS_ML_Processing_Workflow/bin/GNPS2_Subset_Split.py', 'similarity_calculations/spectra_MH_MNA_Translation/merged_pairs.tsv', '/home/user/SourceCode/GNPS_ML_Processing_Workflow/nf_output/spectra_spectra_MH_MNA_Translation.parquet', '/home/user/SourceCode/GNPS_ML_Processing_Workflow/nf_output/summary_spectra_MH_MNA_Translation.csv']

def main():
    _, spectral_similarity, spectra, summary = sys.argv
    
    summary_df = pd.read_csv(summary)
    
    tanimoto_similarity = build_tanimoto_similarity_list_precomputed(summary_df, similarity_threshold = 0.5, output_dir=None)
    spectral_similarity = pd.read_table(spectral_similarity).drop(['DeltaMZ','MinRatio','AlignScore3','MatchedPeaks'], axis=1)

    # Similarity table will be 1 indexed and we will have to map it back to the original files
    ids = summary_df['spectrum_id']
    mapping = {i+1:val for i, val in enumerate(ids)}
    spectral_similarity['spectrumid1'] = spectral_similarity.CLUSTERID1.apply(lambda x: mapping[x])
    spectral_similarity['spectrumid2'] = spectral_similarity.CLUSTERID2.apply(lambda x: mapping[x])
    combined_similarity = pd.merge(spectral_similarity,tanimoto_similarity, on=['spectrumid1','spectrumid2'],how='outer')
    del spectral_similarity
    del tanimoto_similarity
    training_ids, testing_ids = get_splits_index_only(ids, combined_similarity,similarity_metrics=['Cosine','Tanimoto_Similarity'],similarity_threshold=[0.70,0.70],training_fraction=0.70, max_retries=3)
    print(training_ids[0:10],flush=True)
    
    # Build new summary files
    training_mask    = [(x in training_ids) for x in summary_df.spectrum_id]
    testing_mask     = [(x in testing_ids) for x in summary_df.spectrum_id]
    training_summary = summary_df.loc[training_mask]
    testing_summary  = summary_df.loc[testing_mask]
    training_summary.to_csv(summary[:-4]+'_train.csv')
    testing_summary.to_csv(summary[:-4]+'_test.csv')
    
    
    # Build new parquet files
    parquet_as_df = vaex.open(spectra)
    train_df = parquet_as_df[parquet_as_df.spectrum_id.isin(list(training_summary.spectrum_id))]
    train_df.export_parquet(spectra[:-8]+'_train.parquet')
    test_df = parquet_as_df[parquet_as_df.spectrum_id.isin(list(testing_summary.spectrum_id))]
    test_df.export_parquet(spectra[:-8]+'_test.parquet')
        
if __name__ == '__main__':
    main()