import os
import glob
import gzip
import pandas as pd
import numpy as np

def build_mutation_dataset():
    print("Starting 10_build_mutations.py...")
    maf_files = glob.glob('data/raw/mutations/*/*.maf.gz')
    print(f"Found {len(maf_files)} MAF files.")
    
    # We want non-silent / functionally relevant somatic mutations
    nonsilent_classes = {
        'Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 
        'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Splice_Site', 
        'Nonstop_Mutation', 'Translation_Start_Site'
    }
    
    parsed_mutations = []
    sample_summary = []
    
    cols_to_use = [
        'Hugo_Symbol', 'Variant_Classification', 'Variant_Type',
        'Chromosome', 'Start_Position', 'End_Position',
        'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode',
        'HGVSc', 'HGVSp_Short'
    ]
    
    for i, maf_path in enumerate(maf_files, 1):
        if i % 100 == 0 or i == len(maf_files):
            print(f"Processing MAF {i}/{len(maf_files)}...")
            
        try:
            df = pd.read_csv(
                maf_path,
                sep='\t',
                comment='#',
                usecols=lambda c: c in cols_to_use,
                low_memory=False
            )
        except Exception as e:
            print(f"Warning reading {maf_path}: {e}")
            continue
            
        if df.empty:
            continue
            
        sample_barcode = df['Tumor_Sample_Barcode'].iloc[0]
        patient_id = sample_barcode[:12]
        total_mutations = len(df)
        
        # Filter for non-silent alterations
        df_nonsilent = df[df['Variant_Classification'].isin(nonsilent_classes)].copy()
        nonsilent_count = len(df_nonsilent)
        
        sample_summary.append({
            'tumor_sample_barcode': sample_barcode,
            'patient_id': patient_id,
            'total_mutations': total_mutations,
            'nonsilent_mutations': nonsilent_count,
            'maf_file': os.path.basename(maf_path)
        })
        
        if not df_nonsilent.empty:
            df_nonsilent['patient_id'] = patient_id
            parsed_mutations.append(df_nonsilent[['patient_id', 'Tumor_Sample_Barcode', 'Hugo_Symbol', 'Variant_Classification', 'Variant_Type']])

    print("\nConcatenating mutation events...")
    all_mutations_df = pd.concat(parsed_mutations, ignore_index=True)
    summary_df = pd.DataFrame(sample_summary)
    
    print(f"Total mutation events extracted: {len(all_mutations_df):,}")
    print(f"Total unique samples: {summary_df['tumor_sample_barcode'].nunique()}")
    print(f"Total unique patients: {summary_df['patient_id'].nunique()}")
    
    # Save long format mutation table
    os.makedirs('data/processed', exist_ok=True)
    all_mutations_path = 'data/processed/somatic_mutations_long.parquet'
    all_mutations_df.to_parquet(all_mutations_path, index=False)
    print(f"Saved long mutation records to {all_mutations_path}")
    
    # Save sample summary
    sample_summary_path = 'data/processed/mutation_sample_summary.csv'
    summary_df.to_csv(sample_summary_path, index=False)
    print(f"Saved sample summary to {sample_summary_path}")
    
    # Find top 50 most frequently mutated genes in TCGA-BRCA
    patient_gene_pairs = all_mutations_df[['patient_id', 'Hugo_Symbol']].drop_duplicates()
    gene_freq = patient_gene_pairs['Hugo_Symbol'].value_counts()
    print("\nTop 15 Most Mutated Genes in TCGA-BRCA cohort:")
    for gene, count in gene_freq.head(15).items():
        pct = (count / summary_df['patient_id'].nunique()) * 100
        print(f"  {gene:10s}: {count:4d} patients ({pct:.1f}%)")
        
    # Build binary patient x gene mutation matrix for top 200 recurrently mutated genes
    top_genes = gene_freq.head(200).index.tolist()
    matrix_df = pd.crosstab(patient_gene_pairs['patient_id'], patient_gene_pairs['Hugo_Symbol'])
    matrix_df = (matrix_df > 0).astype(int)
    matrix_df = matrix_df.reindex(columns=top_genes, fill_value=0)
    
    # Reindex to all patients in summary_df to ensure zero-mutation patients are included
    all_patients = sorted(summary_df['patient_id'].unique())
    matrix_df = matrix_df.reindex(all_patients, fill_value=0)
    matrix_df.index.name = 'patient_id'
    
    mut_matrix_path = 'data/processed/mutation_matrix_top200.parquet'
    matrix_df.to_parquet(mut_matrix_path)
    print(f"\nSaved binary mutation matrix (patients x top 200 genes) to {mut_matrix_path} with shape {matrix_df.shape}")
    
    # Also save as CSV for easy inspection / Excel export
    matrix_df.to_csv('data/processed/mutation_matrix_top200.csv')
    print("Saved CSV version to data/processed/mutation_matrix_top200.csv")

if __name__ == '__main__':
    build_mutation_dataset()
