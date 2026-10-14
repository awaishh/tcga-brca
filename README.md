# TCGA-BRCA Multi-Omics Data Processing, Quality Audit & Exploratory Analysis

An end-to-end computational pipeline for acquiring, harmonizing, quality-auditing, and exploring multi-omics data from The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA) cohort.

This repository compiles clinical demographics, biospecimen mappings, whole-transcriptome RNA-Seq gene expression, and somatic mutation profiles into unified, analysis-ready datasets tailored for machine learning and survival modeling.

---

## Architecture Overview

```
tcga-brca/
├── data/
│   ├── raw/
│   │   ├── clinical/          # Clinical & follow-up records (TSV)
│   │   ├── biospecimen/      # Sample, aliquot, portion tables (TSV)
│   │   ├── expression/        # 1,219 STAR-augmented gene counts folders
│   │   └── mutations/         # 992 Masked Somatic Mutation MAF files
│   └── processed/
│       ├── patients_master.csv               # Cleaned patient clinical data (1,098 cases)
│       ├── samples.csv                       # Biospecimen barcode mapping (3,397 samples)
│       ├── survival.csv                      # Curated overall survival endpoints (1,074 cases)
│       ├── expression_matrix.parquet         # Full RNA-seq count matrix (60,660 genes x 1,219 samples)
│       ├── expression_log2_cpm.parquet       # Log2(CPM + 1) normalized expression matrix
│       ├── expression_top2000_var_genes.parquet # Top 2,000 highly variable protein-coding genes
│       ├── somatic_mutations_long.parquet    # 67,546 non-silent somatic mutation records
│       ├── mutation_matrix_top200.parquet    # Binary mutation matrix (968 patients x 200 genes)
│       └── tcga_brca_model_ready_cohort.parquet # Multimodal cohort table (1,098 patients x 226 features)
├── reports/
│   ├── figures/                              # High-resolution publication-quality PNG figures
│   ├── TCGA_BRCA_Comprehensive_EDA.ipynb     # Fully executed interactive Jupyter notebook
│   ├── TCGA_BRCA_Comprehensive_EDA.html      # Self-contained standalone HTML EDA report
│   ├── TCGA_BRCA_EDA_and_Data_Correction_Report.md # Formal technical findings & audit report
│   └── data_dictionary.md                    # Detailed field-level data dictionary
├── src/                                      # Sequential ETL & analysis scripts (01 to 12)
├── requirements.txt                          # Python dependencies specification
└── .gitignore                                # Optimized ignore patterns for large genomic assets
```

---

## Dataset Acquisition & Processing Pipeline

The dataset was compiled from the Genomic Data Commons (GDC) repository through a deterministic 13-stage workflow:

### 1. Clinical & Biospecimen Data
* **Deduplication:** Aggregated primary clinical enrollment with longitudinal follow-up records into a single row per patient (`patients_master.csv`, 1,098 unique cases).
* **Sample Typology:** Resolved all TCGA sample codes (`01` = Primary Solid Tumor, `10` = Blood Derived Normal, `11` = Solid Tissue Normal), preserving technical replicates (`samples.csv`, 3,397 rows).
* **Survival Endpoints:** Derived `OS_time` and `OS_event` while systematically filtering non-positive follow-up durations (`survival.csv`, 1,074 patients).

### 2. Whole-Transcriptome RNA-Seq Expression
* **Acquisition:** Downloaded and structured all 1,219 STAR-augmented gene counts TSV files under `data/raw/expression/`.
* **Matrix Assembly:** Built a compressed Parquet matrix of 60,660 genes across 1,219 samples, handling aliquot collisions via unique suffixing (`_rep2`).
* **Normalization:** Computed library-size scaled Counts Per Million (CPM) and variance-stabilizing $\log_2(\text{CPM} + 1)$ transformations (`expression_log2_cpm.parquet`).
* **Feature Selection:** Filtered 19,962 protein-coding genes to the top 2,000 most variable genes across samples (`expression_top2000_var_genes.parquet`).

### 3. Somatic Mutation Landscape
* **Acquisition:** Downloaded and verified all 992 Masked Somatic Mutation MAF files against the GDC manifest using exact byte matching (`data/raw/mutations/`).
* **Variant Classification:** Extracted 67,546 non-silent functional mutations (Missense, Nonsense, Frame-shift, In-frame, Splice-site) across 968 patients.
* **Driver Alterations:** Quantified mutation frequencies for primary drivers (*TP53*: 34.2%, *PIK3CA*: 33.9%, *TTN*: 17.3%, *CDH1*: 13.2%, *GATA3*: 12.9%).
* **Binary Matrix:** Constructed an aligned patient-by-gene alteration matrix for the top 200 recurrently mutated genes (`mutation_matrix_top200.parquet`).

### 4. Multimodal Integration & Data Quality Corrections
* Assembled the unified `tcga_brca_model_ready_cohort.parquet` table (1,098 patients $\times$ 226 curated features).
* Imputed missing clinical age values with the cohort median (58.0 years) while preserving original attributes.
* Consolidated granular AJCC pathologic stages into standardized macro-categories (Stage I, II, III, IV, Unknown).

---

## Exploratory Data Analysis & Visualizations

The executed analysis produces six core visual diagnostics in `reports/figures/`:

1. **Clinical Demographics & Staging (`01_clinical_demographics_and_staging.png`):**
   * Distribution of age at diagnosis (median: 58.0 years, range: 26–90 years).
   * Breakdown of simplified tumor stages: Stage I (16.0%), Stage II (56.0%), Stage III (22.5%), Stage IV (1.8%).
   * Major histology representation: Infiltrating Ductal (~70%) vs. Infiltrating Lobular (~18%).
   * Overall survival distributions across pathologic stages.

2. **Somatic Mutation Landscape (`02_somatic_mutation_landscape.png`):**
   * Frequency barplot of the top 20 mutated genes in the TCGA-BRCA cohort.
   * Mutual exclusivity correlation matrix highlighting distinct separation between *TP53* (Basal/HER2) and *PIK3CA* / *GATA3* (Luminal) mutations.

3. **Tumor Mutation Burden (`03_tumor_mutation_burden.png`):**
   * Non-silent somatic mutation burden distribution on a $\log_{10}$ scale.
   * Stratification of genomic instability demonstrating significantly higher mutation burden in *TP53*-mutated tumors ($p < 10^{-15}$).

4. **Transcriptomic Architecture via PCA (`04_expression_pca_transcriptomics.png`):**
   * Unsupervised Principal Component Analysis of top 2,000 variable genes.
   * Distinct biological segregation of Primary Tumors from Solid Tissue Normal controls along PC1 (14.6% explained variance).
   * Subtype divergence along PC2 separating *TP53*-mutant from wild-type tumors.

5. **Canonical Biomarkers Clustermap (`05_canonical_biomarkers_clustermap.png`):**
   * Hierarchical biclustering of canonical PAM50 marker genes (*ESR1*, *PGR*, *ERBB2*, *MKI67*, *FOXA1*, *GATA3*, *EGFR*, *KRT5*, *CDH1*), delineating Luminal, HER2-enriched, and Basal-like signatures.

6. **Overall Survival Projections (`06_survival_analysis_kaplan_meier.png`):**
   * Empirical Kaplan-Meier survival curves stratified by AJCC Stage.
   * Survival divergence between *TP53*-mutant and *TP53*-wild-type cohorts over 10-year follow-up.

---

## Installation & Environment Setup

Clone the repository and install dependencies inside an isolated Python virtual environment:

```bash
# Clone repository
git clone <your-repository-url>
cd tcga-brca

# Create virtual environment
python -m venv venv

# Activate virtual environment
# On Windows:
.\venv\Scripts\activate
# On Linux/macOS:
source venv/bin/activate

# Install required packages
pip install -r requirements.txt
```

---

## Execution Guide

### Reproduce Figures & Pipeline
To re-run the complete EDA pipeline and regenerate all figures:
```bash
python src/12_comprehensive_eda.py
```

### Launch Interactive Analysis Notebook
Open the pre-computed notebook in Jupyter or your preferred IDE:
```bash
jupyter notebook reports/TCGA_BRCA_Comprehensive_EDA.ipynb
```

Alternatively, open `reports/TCGA_BRCA_Comprehensive_EDA.html` directly in any web browser to view the complete executed report without starting a kernel.

---

## Data Governance & Version Control Notes

Raw genomic data files (`*.tsv`, `*.maf.gz`, `*.parquet`) and client executables exceed GitHub's file size limits and are tracked via `.gitignore`. 

Processed tables and matrices can be regenerated deterministically from the scripts in `src/`. For downstream machine learning models, load the parquet files using standard pandas/pyarrow engines:

```python
import pandas as pd

# Load unified multimodal cohort
cohort = pd.read_parquet('data/processed/tcga_brca_model_ready_cohort.parquet')

# Load normalized gene expression (top 2,000 variable genes)
expr = pd.read_parquet('data/processed/expression_top2000_var_genes.parquet')
```
