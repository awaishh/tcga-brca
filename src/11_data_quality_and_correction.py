import os
import glob
import pandas as pd
import numpy as np

def clean_and_integrate():
    print("Starting 11_data_quality_and_correction.py...")
    
    # -------------------------------------------------------------
    # 1. Gene Metadata & Annotation
    # -------------------------------------------------------------
    print("\n[Step 1] Extracting gene annotations (symbol, biotype)...")
    sample_tsv = glob.glob('data/raw/expression/*/*.augmented_star_gene_counts.tsv')[0]
    gene_meta = pd.read_csv(sample_tsv, sep='\t', skiprows=1, usecols=['gene_id', 'gene_name', 'gene_type'])
    # Filter out special summary rows like N_unmapped, N_multimapping, etc.
    gene_meta = gene_meta[~gene_meta['gene_id'].str.startswith('N_')].reset_index(drop=True)
    gene_meta.to_parquet('data/processed/gene_metadata.parquet', index=False)
    print(f"Extracted {len(gene_meta):,} genes. Protein coding genes: {(gene_meta['gene_type'] == 'protein_coding').sum():,}")

    # -------------------------------------------------------------
    # 2. Expression Matrix Normalization & Log2 Transformation
    # -------------------------------------------------------------
    print("\n[Step 2] Normalizing Gene Expression (CPM + log2(x+1))...")
    expr_df = pd.read_parquet('data/processed/expression_matrix.parquet')
    gene_ids = expr_df['gene_id']
    sample_cols = [c for c in expr_df.columns if c != 'gene_id']
    
    raw_counts = expr_df[sample_cols].astype(np.float32)
    # Library size (counts per sample)
    lib_sizes = raw_counts.sum(axis=0)
    
    # Counts Per Million (CPM)
    cpm = (raw_counts / lib_sizes) * 1e6
    # Log2(CPM + 1) transform for standard machine learning / statistical modeling
    log2_cpm = np.log2(cpm + 1.0)
    log2_cpm.insert(0, 'gene_id', gene_ids)
    
    log2_path = 'data/processed/expression_log2_cpm.parquet'
    log2_cpm.to_parquet(log2_path, index=False)
    print(f"Successfully saved Log2-CPM expression matrix to {log2_path} with shape {log2_cpm.shape}")
    
    # Also extract expression for top variable protein-coding genes (useful for ML baseline)
    pc_genes = set(gene_meta[gene_meta['gene_type'] == 'protein_coding']['gene_id'])
    is_pc = log2_cpm['gene_id'].isin(pc_genes)
    log2_pc = log2_cpm[is_pc].copy()
    
    # Calculate gene variances across samples
    variances = log2_pc[sample_cols].var(axis=1)
    top2000_idx = variances.nlargest(2000).index
    top2000_expr = log2_pc.loc[top2000_idx]
    top2000_path = 'data/processed/expression_top2000_var_genes.parquet'
    top2000_expr.to_parquet(top2000_path, index=False)
    print(f"Saved top 2,000 highly variable protein-coding genes to {top2000_path}")

    # -------------------------------------------------------------
    # 3. Clinical Data Cleansing & Harmonization
    # -------------------------------------------------------------
    print("\n[Step 3] Harmonizing and Cleansing Clinical & Survival Features...")
    patients_df = pd.read_csv('data/processed/patients_master.csv')
    survival_df = pd.read_csv('data/processed/survival.csv')
    
    # Merge survival outcomes
    clin_merged = pd.merge(patients_df, survival_df, on='patient_id', how='left')
    
    # Clean Stage: Simplify AJCC Stage (Stage I, II, III, IV, or Unknown)
    def simplify_stage(val):
        if pd.isna(val) or val == 'Stage X':
            return 'Unknown'
        val = str(val).strip()
        if val.startswith('Stage I') and not val.startswith('Stage II') and not val.startswith('Stage IV'):
            if val.startswith('Stage III'):
                return 'Stage III'
            return 'Stage I'
        elif val.startswith('Stage II') and not val.startswith('Stage III'):
            return 'Stage II'
        elif val.startswith('Stage III'):
            return 'Stage III'
        elif val.startswith('Stage IV'):
            return 'Stage IV'
        return 'Unknown'

    clin_merged['tumor_stage_simplified'] = clin_merged['ajcc_pathologic_stage'].apply(simplify_stage)
    
    # Impute missing age with median (for ML readiness) while preserving original
    median_age = clin_merged['age_at_diagnosis'].median()
    clin_merged['age_imputed'] = clin_merged['age_at_diagnosis'].fillna(median_age)
    
    # Standardize vital status
    clin_merged['is_deceased'] = (clin_merged['vital_status'].str.lower() == 'dead').astype(int)
    
    # Clean T, N, M categories
    clin_merged['pathologic_t_simple'] = clin_merged['ajcc_pathologic_t'].fillna('TX').str[:2]
    clin_merged['pathologic_n_simple'] = clin_merged['ajcc_pathologic_n'].fillna('NX').str[:2]
    clin_merged['pathologic_m_simple'] = clin_merged['ajcc_pathologic_m'].fillna('MX').str[:2]
    
    # -------------------------------------------------------------
    # 4. Integrate Primary Tumor Expression & Mutation Multi-Omics Table
    # -------------------------------------------------------------
    print("\n[Step 4] Assembling Patient-Level Multi-Omics Modeling Table...")
    # Load mutation matrix (top 200 genes)
    mut_df = pd.read_parquet('data/processed/mutation_matrix_top200.parquet')
    mut_df.columns = [f"mut_{c}" for c in mut_df.columns]
    
    # Load mutation summary for tumor mutation burden (TMB proxy)
    mut_summary = pd.read_csv('data/processed/mutation_sample_summary.csv')
    patient_tmb = mut_summary.groupby('patient_id')[['total_mutations', 'nonsilent_mutations']].max().reset_index()
    patient_tmb.rename(columns={
        'total_mutations': 'tmb_total_mutations',
        'nonsilent_mutations': 'tmb_nonsilent_mutations'
    }, inplace=True)
    
    # Merge clinical with TMB and key recurrent mutations (TP53, PIK3CA, CDH1, GATA3, etc.)
    multi_omics_df = pd.merge(clin_merged, patient_tmb, on='patient_id', how='left')
    multi_omics_df = pd.merge(multi_omics_df, mut_df, on='patient_id', how='left')
    
    # Fill mutation values with 0 for patients who had no mutations recorded or were not sequenced
    mut_cols = [c for c in multi_omics_df.columns if c.startswith('mut_')]
    multi_omics_df[mut_cols] = multi_omics_df[mut_cols].fillna(0).astype(int)
    multi_omics_df['tmb_total_mutations'] = multi_omics_df['tmb_total_mutations'].fillna(0).astype(int)
    multi_omics_df['tmb_nonsilent_mutations'] = multi_omics_df['tmb_nonsilent_mutations'].fillna(0).astype(int)
    
    # Add flag for expression data availability
    # Sample barcode prefix is patient ID
    expr_patients = set(c[:12] for c in sample_cols)
    multi_omics_df['has_rna_seq'] = multi_omics_df['patient_id'].isin(expr_patients).astype(int)
    
    # Save unified analysis-ready dataset
    model_ready_path = 'data/processed/tcga_brca_model_ready_cohort.csv'
    multi_omics_df.to_csv(model_ready_path, index=False)
    print(f"Unified model-ready cohort table saved to {model_ready_path} with shape {multi_omics_df.shape}")
    
    # Also save as Parquet for fast loading
    multi_omics_df.to_parquet('data/processed/tcga_brca_model_ready_cohort.parquet', index=False)
    print("Saved Parquet copy to data/processed/tcga_brca_model_ready_cohort.parquet")
    
    print("\nData cleansing and multi-omics alignment complete!")

if __name__ == '__main__':
    clean_and_integrate()
