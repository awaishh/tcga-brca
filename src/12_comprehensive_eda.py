import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA

# Styling configuration
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams['axes.edgecolor'] = '#333333'
plt.rcParams['axes.linewidth'] = 0.8
plt.rcParams['grid.alpha'] = 0.3

FIGURES_DIR = 'reports/figures'
os.makedirs(FIGURES_DIR, exist_ok=True)

def generate_eda_figures():
    print("Starting 12_comprehensive_eda.py...")
    
    # -------------------------------------------------------------
    # 1. Load Data
    # -------------------------------------------------------------
    print("Loading multi-omics cohort and processed tables...")
    cohort = pd.read_csv('data/processed/tcga_brca_model_ready_cohort.csv')
    expr_top2000 = pd.read_parquet('data/processed/expression_top2000_var_genes.parquet')
    gene_meta = pd.read_parquet('data/processed/gene_metadata.parquet')
    gene_map = dict(zip(gene_meta['gene_id'], gene_meta['gene_name']))
    
    # -------------------------------------------------------------
    # Figure 1: Clinical Demographics & Disease Staging Dashboard
    # -------------------------------------------------------------
    print("\n[Fig 1] Generating Clinical Demographics & Staging Dashboard...")
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    
    # Age distribution
    sns.histplot(cohort['age_imputed'], kde=True, ax=axes[0, 0], color='#2b5c8f', bins=25)
    axes[0, 0].set_title('A. Age at Diagnosis Distribution (Median: {:.1f} yrs)'.format(cohort['age_imputed'].median()), fontsize=12, fontweight='bold')
    axes[0, 0].set_xlabel('Age (Years)')
    axes[0, 0].set_ylabel('Patient Count')
    axes[0, 0].axvline(cohort['age_imputed'].median(), color='#d95f02', linestyle='--', label='Median')
    axes[0, 0].legend()
    
    # Tumor Stage
    stage_counts = cohort['tumor_stage_simplified'].value_counts()
    order = ['Stage I', 'Stage II', 'Stage III', 'Stage IV', 'Unknown']
    stage_counts = stage_counts.reindex(order).fillna(0)
    sns.barplot(x=stage_counts.index, y=stage_counts.values, hue=stage_counts.index, ax=axes[0, 1], palette='Blues_r', legend=False)
    axes[0, 1].set_title('B. Pathologic Tumor Stage Distribution', fontsize=12, fontweight='bold')
    axes[0, 1].set_xlabel('Tumor Stage')
    axes[0, 1].set_ylabel('Patient Count')
    for i, v in enumerate(stage_counts.values):
        axes[0, 1].text(i, v + 8, f"{int(v)} ({v/len(cohort)*100:.1f}%)", ha='center', fontsize=9)
        
    # Primary Diagnosis Subtypes
    top_histology = cohort['primary_diagnosis'].value_counts().head(5)
    sns.barplot(y=top_histology.index, x=top_histology.values, hue=top_histology.index, ax=axes[1, 0], palette='viridis', legend=False)
    axes[1, 0].set_title('C. Major Histological Subtypes', fontsize=12, fontweight='bold')
    axes[1, 0].set_xlabel('Patient Count')
    axes[1, 0].set_ylabel('')
    
    # Overall Survival Days by Stage (Boxplot)
    surv_stage = cohort.dropna(subset=['OS_time'])
    sns.boxplot(data=surv_stage, x='tumor_stage_simplified', y='OS_time', order=order, ax=axes[1, 1], palette='Set2')
    axes[1, 1].set_title('D. Overall Survival Time by Simplified Stage', fontsize=12, fontweight='bold')
    axes[1, 1].set_xlabel('Stage')
    axes[1, 1].set_ylabel('Overall Survival (Days)')
    
    plt.tight_layout()
    fig1_path = os.path.join(FIGURES_DIR, '01_clinical_demographics_and_staging.png')
    plt.savefig(fig1_path, dpi=300)
    plt.close()
    print(f"Saved: {fig1_path}")

    # -------------------------------------------------------------
    # Figure 2: Somatic Mutation Landscape & Co-occurrence
    # -------------------------------------------------------------
    print("\n[Fig 2] Generating Somatic Mutation Landscape & Oncoplot/Co-occurrence...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # Top 20 Mutated Genes Barplot
    mut_cols = [c for c in cohort.columns if c.startswith('mut_')]
    mut_freq = cohort[mut_cols].mean().sort_values(ascending=False).head(20)
    mut_freq.index = [c.replace('mut_', '') for c in mut_freq.index]
    
    sns.barplot(x=mut_freq.values * 100, y=mut_freq.index, hue=mut_freq.index, ax=axes[0], palette='magma', legend=False)
    axes[0].set_title('A. Top 20 Frequently Mutated Genes in TCGA-BRCA', fontsize=12, fontweight='bold')
    axes[0].set_xlabel('Patient Mutation Frequency (%)')
    axes[0].set_ylabel('Gene Symbol')
    for i, v in enumerate(mut_freq.values * 100):
        axes[0].text(v + 0.5, i, f"{v:.1f}%", va='center', fontsize=9)
        
    # Top 10 Drivers Co-occurrence / Mutual Exclusivity Heatmap
    top10_mut = [f"mut_{g}" for g in mut_freq.head(10).index]
    cooccur_corr = cohort[top10_mut].corr()
    cooccur_corr.columns = [c.replace('mut_', '') for c in cooccur_corr.columns]
    cooccur_corr.index = [c.replace('mut_', '') for c in cooccur_corr.index]
    
    mask = np.triu(np.ones_like(cooccur_corr, dtype=bool))
    sns.heatmap(cooccur_corr, mask=mask, cmap='vlag', center=0, annot=True, fmt='.2f', 
                square=True, ax=axes[1], cbar_kws={'label': 'Pearson Correlation'})
    axes[1].set_title('B. Somatic Mutation Co-occurrence / Exclusivity Correlation', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    fig2_path = os.path.join(FIGURES_DIR, '02_somatic_mutation_landscape.png')
    plt.savefig(fig2_path, dpi=300)
    plt.close()
    print(f"Saved: {fig2_path}")

    # -------------------------------------------------------------
    # Figure 3: Tumor Mutational Burden (TMB) vs Survival & Molecular Subtypes
    # -------------------------------------------------------------
    print("\n[Fig 3] Generating Tumor Mutation Burden (TMB) Distribution...")
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # TMB log distribution
    tmb_series = cohort[cohort['tmb_nonsilent_mutations'] > 0]['tmb_nonsilent_mutations']
    sns.histplot(np.log10(tmb_series), kde=True, ax=axes[0], color='#8856a7', bins=30)
    axes[0].set_title('A. Log10 Non-Silent Mutation Burden Distribution', fontsize=12, fontweight='bold')
    axes[0].set_xlabel('log10(Non-Silent Mutation Count)')
    axes[0].set_ylabel('Sample Count')
    
    # TMB by TP53 mutation status
    cohort['TP53_Status'] = cohort['mut_TP53'].map({1: 'TP53 Mutated', 0: 'TP53 Wild-Type'})
    sns.boxplot(data=cohort[cohort['tmb_nonsilent_mutations'] > 0], x='TP53_Status', y=np.log10(cohort['tmb_nonsilent_mutations'] + 1), ax=axes[1], palette=['#7570b3', '#1b9e77'])
    axes[1].set_title('B. Mutation Burden by TP53 Alteration Status', fontsize=12, fontweight='bold')
    axes[1].set_xlabel('')
    axes[1].set_ylabel('log10(Non-Silent Mutation Count + 1)')
    
    plt.tight_layout()
    fig3_path = os.path.join(FIGURES_DIR, '03_tumor_mutation_burden.png')
    plt.savefig(fig3_path, dpi=300)
    plt.close()
    print(f"Saved: {fig3_path}")

    # -------------------------------------------------------------
    # Figure 4: Gene Expression PCA & Molecular Diversity
    # -------------------------------------------------------------
    print("\n[Fig 4] Computing Principal Component Analysis (PCA) on Top Variable Genes...")
    # Get expression samples
    sample_cols = [c for c in expr_top2000.columns if c != 'gene_id']
    X = expr_top2000[sample_cols].values.T # shape: (samples, 2000 genes)
    
    pca = PCA(n_components=2, random_state=42)
    pca_coords = pca.fit_transform(X)
    var_exp = pca.explained_variance_ratio_ * 100
    
    pca_df = pd.DataFrame({
        'sample_barcode': sample_cols,
        'patient_id': [c[:12] for c in sample_cols],
        'PC1': pca_coords[:, 0],
        'PC2': pca_coords[:, 1],
    })
    
    # Classify sample tissue type from barcode (01=Primary Solid Tumor, 11=Solid Tissue Normal)
    pca_df['tissue_type'] = pca_df['sample_barcode'].apply(
        lambda s: 'Primary Tumor' if '-01' in s else ('Normal Tissue' if '-11' in s else 'Other')
    )
    
    # Merge TP53 mutation status
    pca_df = pd.merge(pca_df, cohort[['patient_id', 'TP53_Status', 'tumor_stage_simplified']], on='patient_id', how='left')
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    
    # PCA colored by Tissue Type (Tumor vs Normal Separation)
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue='tissue_type', palette={'Primary Tumor': '#e41a1c', 'Normal Tissue': '#377eb8', 'Other': '#999999'}, alpha=0.7, s=40, ax=axes[0])
    axes[0].set_title(f'A. Transcriptomic PCA: Primary Tumor vs Normal Separation\n(PC1: {var_exp[0]:.1f}%, PC2: {var_exp[1]:.1f}% var)', fontsize=12, fontweight='bold')
    axes[0].set_xlabel(f'PC1 ({var_exp[0]:.1f}%)')
    axes[0].set_ylabel(f'PC2 ({var_exp[1]:.1f}%)')
    axes[0].legend(title='Tissue Type')
    
    # PCA colored by TP53 status
    tumor_pca = pca_df[pca_df['tissue_type'] == 'Primary Tumor']
    sns.scatterplot(data=tumor_pca, x='PC1', y='PC2', hue='TP53_Status', palette={'TP53 Mutated': '#d95f02', 'TP53 Wild-Type': '#1b9e77'}, alpha=0.7, s=40, ax=axes[1])
    axes[1].set_title('B. Primary Tumor Transcriptome Stratified by TP53 Mutation', fontsize=12, fontweight='bold')
    axes[1].set_xlabel(f'PC1 ({var_exp[0]:.1f}%)')
    axes[1].set_ylabel(f'PC2 ({var_exp[1]:.1f}%)')
    axes[1].legend(title='TP53 Status')
    
    plt.tight_layout()
    fig4_path = os.path.join(FIGURES_DIR, '04_expression_pca_transcriptomics.png')
    plt.savefig(fig4_path, dpi=300)
    plt.close()
    print(f"Saved: {fig4_path}")

    # -------------------------------------------------------------
    # Figure 5: Key Breast Cancer Marker Gene Expression Heatmap
    # -------------------------------------------------------------
    print("\n[Fig 5] Generating Breast Cancer PAM50 & Canonical Markers Heatmap...")
    # Canonical breast cancer markers: ESR1 (ER), PGR (PR), ERBB2 (HER2), MKI67 (proliferation), CDH1, FOXA1, GATA3, TP53, EGFR, KRT5
    canonical_markers = ['ESR1', 'PGR', 'ERBB2', 'MKI67', 'CDH1', 'FOXA1', 'GATA3', 'TP53', 'EGFR', 'KRT5', 'BCL2', 'MYC']
    inv_gene_map = {v: k for k, v in gene_map.items()}
    marker_ids = [inv_gene_map[g] for g in canonical_markers if g in inv_gene_map]
    
    expr_log2_cpm = pd.read_parquet('data/processed/expression_log2_cpm.parquet')
    marker_expr = expr_log2_cpm[expr_log2_cpm['gene_id'].isin(marker_ids)].copy()
    marker_expr['gene_symbol'] = marker_expr['gene_id'].map(gene_map)
    marker_expr.set_index('gene_symbol', inplace=True)
    marker_expr.drop(columns=['gene_id'], inplace=True)
    
    # Select primary tumor samples (first 100 for clean visualization)
    tumor_samples = [c for c in marker_expr.columns if '-01' in c][:100]
    marker_matrix = marker_expr[tumor_samples]
    
    # Z-score normalize across samples
    z_marker = marker_matrix.apply(lambda x: (x - x.mean()) / (x.std() + 1e-6), axis=1)
    
    plt.figure(figsize=(14, 8))
    g = sns.clustermap(z_marker, cmap='coolwarm', center=0, col_cluster=True, row_cluster=True,
                       figsize=(15, 8), cbar_kws={'label': 'Expression Z-Score'},
                       yticklabels=True, xticklabels=False)
    g.fig.suptitle('Canonical Breast Cancer Biomarkers Expression Landscape (PAM50/Subtype Signatures)', fontsize=14, fontweight='bold', y=1.02)
    fig5_path = os.path.join(FIGURES_DIR, '05_canonical_biomarkers_clustermap.png')
    g.savefig(fig5_path, dpi=300)
    plt.close()
    print(f"Saved: {fig5_path}")

    # -------------------------------------------------------------
    # Figure 6: Overall Survival Kaplan-Meier Style Survival Curves
    # -------------------------------------------------------------
    print("\n[Fig 6] Generating Overall Survival Analysis...")
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    
    # 1. Survival stratified by Stage
    for stg in ['Stage I', 'Stage II', 'Stage III', 'Stage IV']:
        subset = cohort[cohort['tumor_stage_simplified'] == stg]
        times = np.sort(subset['OS_time'].dropna().values)
        if len(times) == 0:
            continue
        # Empirical survival function
        surv_prob = 1.0 - (np.arange(len(times)) / len(times))
        axes[0].step(times, surv_prob, where='post', label=f"{stg} (n={len(times)})", linewidth=2)
        
    axes[0].set_title('A. Overall Survival Stratified by Stage', fontsize=12, fontweight='bold')
    axes[0].set_xlabel('Time (Days)')
    axes[0].set_ylabel('Empirical Overall Survival Probability')
    axes[0].set_ylim(-0.05, 1.05)
    axes[0].legend()
    
    # 2. Survival stratified by TP53 mutation
    for status, color in [('TP53 Mutated', '#d95f02'), ('TP53 Wild-Type', '#1b9e77')]:
        subset = cohort[cohort['TP53_Status'] == status]
        times = np.sort(subset['OS_time'].dropna().values)
        if len(times) == 0:
            continue
        surv_prob = 1.0 - (np.arange(len(times)) / len(times))
        axes[1].step(times, surv_prob, where='post', label=f"{status} (n={len(times)})", color=color, linewidth=2)
        
    axes[1].set_title('B. Overall Survival Stratified by TP53 Mutation Status', fontsize=12, fontweight='bold')
    axes[1].set_xlabel('Time (Days)')
    axes[1].set_ylabel('Empirical Overall Survival Probability')
    axes[1].set_ylim(-0.05, 1.05)
    axes[1].legend()
    
    plt.tight_layout()
    fig6_path = os.path.join(FIGURES_DIR, '06_survival_analysis_kaplan_meier.png')
    plt.savefig(fig6_path, dpi=300)
    plt.close()
    print(f"Saved: {fig6_path}")

    print("\nAll EDA visualizations successfully produced in reports/figures/!")

if __name__ == '__main__':
    generate_eda_figures()
