# TCGA-BRCA Data Pipeline — Complete Lineage Report

> **Author**: Capstone
> **Pipeline**: 12 scripts (src/01_data_inventory.py to src/12_comprehensive_eda.py)
> **Final Output**: tcga_brca_model_ready_cohort.parquet — **1,098 patients x 226 columns**

---

## 1. Column Counts Per Dataset

| Dataset | Rows | Columns | File |
|---|---|---|---|
| patients_master.csv | 1,098 | **15** | data/processed/patients_master.csv |
| survival.csv | 1,074 | **3** | data/processed/survival.csv |
| samples.csv | 3,397 | **4** | data/processed/samples.csv |
| expression_matrix.parquet | 60,660 | **1,220** | data/processed/expression_matrix.parquet |
| expression_top2000_var_genes.parquet | 2,000 | **1,220** | data/processed/expression_top2000_var_genes.parquet |
| somatic_mutations_long.parquet | 67,546 | **5** | data/processed/somatic_mutations_long.parquet |
| mutation_matrix_top200.parquet | 968 | **200** | data/processed/mutation_matrix_top200.parquet |
| **tcga_brca_model_ready_cohort.parquet** | **1,098** | **226** | data/processed/tcga_brca_model_ready_cohort.parquet |

### Breakdown of the 226 Columns in the Final Cohort

| Column Group | Count | Source |
|---|---|---|
| Clinical demographics | 15 | patients_master.csv |
| Survival endpoints (OS_time, OS_event) | 2 | survival.csv |
| Derived clinical (tumor_stage_simplified, age_imputed, is_deceased, T/N/M simplified) | 6 | Step 11 |
| Tumor mutation burden (tmb_total, tmb_nonsilent) | 2 | mutation_sample_summary.csv |
| Binary mutation flags (mut_TP53, mut_PIK3CA, ... mut_HTT) | **200** | mutation_matrix_top200.parquet |
| RNA-Seq availability flag (has_rna_seq) | 1 | Step 11 |
| **Total** | **226** | |

---

## 2. Column Renames

All column renames happen in `src/02_patients_master.py` and `src/10_build_mutations.py`.

### Clinical Renames (Step 02)

| Original TCGA Column | Renamed To |
|---|---|
| `cases.submitter_id` | `patient_id` |
| `diagnoses.age_at_diagnosis` | `age_at_diagnosis` |
| `demographic.sex_at_birth` | `gender` |
| `demographic.race` | `race` |
| `demographic.ethnicity` | `ethnicity` |
| `diagnoses.ajcc_pathologic_stage` | `ajcc_pathologic_stage` |
| `diagnoses.ajcc_pathologic_t` | `ajcc_pathologic_t` |
| `diagnoses.ajcc_pathologic_n` | `ajcc_pathologic_n` |
| `diagnoses.ajcc_pathologic_m` | `ajcc_pathologic_m` |
| `diagnoses.primary_diagnosis` | `primary_diagnosis` |
| `diagnoses.morphology` | `morphology` |
| `demographic.vital_status` | `vital_status` |
| `demographic.days_to_death` | `days_to_death` |
| `diagnoses.days_to_last_follow_up` | `days_to_last_follow_up_clinical` |
| `follow_ups.days_to_follow_up` | `days_to_last_follow_up` |

### Mutation Column Renames (Step 10 / Step 11)

| Original | Renamed To |
|---|---|
| `total_mutations` | `tmb_total_mutations` |
| `nonsilent_mutations` | `tmb_nonsilent_mutations` |
| Gene names (e.g., `TP53`, `PIK3CA`) | Prefixed with `mut_` (e.g., `mut_TP53`, `mut_PIK3CA`) |

### Derived Columns Created in Step 11

| New Column | Logic |
|---|---|
| `tumor_stage_simplified` | Maps AJCC stages (IA, IB, IIA, etc.) to Stage I / II / III / IV / Unknown |
| `age_imputed` | `age_at_diagnosis` with NaN filled by median (21,560.5 days) |
| `is_deceased` | Binary: `vital_status == 'Dead'` mapped to 1, else 0 |
| `pathologic_t_simple` | First 2 characters of `ajcc_pathologic_t` (e.g., T1, T2); NaN becomes TX |
| `pathologic_n_simple` | First 2 characters of `ajcc_pathologic_n`; NaN becomes NX |
| `pathologic_m_simple` | First 2 characters of `ajcc_pathologic_m`; NaN becomes MX |
| `has_rna_seq` | Binary: 1 if patient has RNA-Seq expression data, else 0 |

---

## 3. Missing Value Assumptions

### What was replaced with NaN

In `src/02_patients_master.py`, the following raw string values were treated as missing:

```python
missing_vals = ["'--", "Not Reported", "not reported", "Unknown"]
# All replaced with np.nan
```

### What was filled / imputed

| Column | Action | Assumption |
|---|---|---|
| `age_imputed` | Filled 16 NaN values with **median** (21,560.5 days) | Median imputation is standard for ML; original `age_at_diagnosis` preserved alongside for reference |
| `pathologic_t_simple` | NaN filled with `"TX"` | TX = "stage cannot be evaluated" (standard AJCC convention) |
| `pathologic_n_simple` | NaN filled with `"NX"` | NX = "regional nodes cannot be assessed" |
| `pathologic_m_simple` | NaN filled with `"MX"` | MX = "distant metastasis cannot be assessed" |
| `tumor_stage_simplified` | NaN / Stage X mapped to `"Unknown"` | Preserves the information that stage was unavailable |
| All `mut_*` columns (200) | NaN filled with **0** | Assumption: if a patient has no MAF file or no mutation in that gene, they are treated as wild-type (unmutated) for modeling purposes |
| `tmb_total_mutations` | NaN filled with **0** | Same as above |
| `tmb_nonsilent_mutations` | NaN filled with **0** | Same as above |
| `sample_type` (Step 04) | Unmapped codes kept as raw code | No data fabricated; only known TCGA codes (01, 06, 10, 11) are labeled |

### What was left as NaN (not imputed)

| Column | NaN Count | Reason |
|---|---|---|
| `age_at_diagnosis` | 16 | Original raw value preserved; imputed copy in `age_imputed` |
| `gender` | 1 | Cannot assume sex; left for user to handle |
| `race` | 96 | Cannot assume race; left as-is |
| `ethnicity` | 175 | Cannot assume ethnicity; left as-is |
| `ajcc_pathologic_stage` | 12 | Simplified version in `tumor_stage_simplified` (= "Unknown"); raw value preserved |
| `days_to_death` | 947 | Expected: only 151 patients are deceased; the 947 alive patients have no death date |
| `days_to_last_follow_up_clinical` | 1,097 | This is the clinical.tsv follow-up; nearly all NaN because TCGA clinical rarely populates it |
| `OS_time` | 24 | Patients with no valid survival time (see Section 4) |
| `OS_event` | 24 | Same patients as above |

---

## 4. Rows Deleted (and Why)

**No rows were arbitrarily deleted. Every deletion has a documented reason.**

### Step 02 — Clinical Deduplication

| Action | Rows Affected | Reason |
|---|---|---|
| Filter to `diagnoses.diagnosis_is_primary_disease == True` | Removes secondary/metastatic diagnosis rows | Each patient should have one row for their primary breast cancer diagnosis |
| `drop_duplicates(subset='cases.submitter_id', keep='first')` | Removes residual duplicate treatment rows | Keeps the first clinical record per patient |
| Fallback merge for filtered-out patients | 0 lost | If a patient had no "primary" flag, their first row from the original was kept |

**Result**: Raw clinical.tsv (multiple rows per patient) collapsed to exactly **1,098 unique patients** (1 row each).

### Step 02 — Follow-up Deduplication

| Action | Rows Affected | Reason |
|---|---|---|
| Sort by `days_to_follow_up` descending, then `drop_duplicates(keep='first')` | Keeps only the **latest** follow-up per patient | Multiple follow-up visits exist; we want the most recent survival status |

### Step 03 — Survival Endpoint Cleaning

| Dropped Rows | Count | Reason |
|---|---|---|
| Null `OS_time` | ~22 | No death date and no follow-up date available; survival analysis impossible |
| Negative `OS_time` | ~0 | Data error if present |
| Zero `OS_time` | ~2 | Enrolled and immediately lost; uninformative for survival curves |
| **Total dropped** | **24** | From 1,098 to 1,074 patients with valid survival |

### Step 08 — Expression Matrix

| Action | Reason |
|---|---|
| Rows starting with `N_` (e.g., `N_unmapped`, `N_multimapping`) dropped | These are QC summary rows from STAR, not gene expression values |

### Step 10 — Somatic Mutations

| Action | Reason |
|---|---|
| Silent mutations filtered out | Only **non-silent** variant classes kept (Missense, Nonsense, Frame_Shift, Splice_Site, etc.) for functional relevance |
| Empty MAF files skipped | If a file had zero mutation calls, it contributes nothing |

---

## 5. Complete Data Flow

```
RAW INPUTS                         PROCESSING STEPS                    FINAL OUTPUTS
==========                         ================                    =============

clinical.tsv ──────┐
(multi-row/patient) │
                    ├─► Step 02 ──► patients_master.csv (1,098 x 15)
follow_up.tsv ──────┘    │              │
(multi-visit/patient)    │              │
                         │              ▼
                         │         Step 03 ──► survival.csv (1,074 x 3)
                         │              │
sample.tsv ──────────► Step 04 ──► samples.csv (3,397 x 4)
                         │
1,219 Expression   ──┐   │
TSVs (RNA-Seq)       ├─► Step 08 ──► expression_matrix.parquet (60,660 x 1,220)
gdc_sample_sheet.tsv ┘   │              │
                         │              ▼
992 MAF.gz files ────► Step 10 ──► somatic_mutations_long.parquet (67,546 x 5)
(somatic mutations)      │         mutation_matrix_top200.parquet (968 x 200)
                         │              │
                         │              ▼
                         └────────► Step 11 (Integration & Quality)
                                       │
                                       ├──► expression_log2_cpm.parquet (60,660 x 1,220)
                                       ├──► expression_top2000_var_genes.parquet (2,000 x 1,220)
                                       └──► tcga_brca_model_ready_cohort.parquet (1,098 x 226)
                                                        │
                                                        ▼
                                                   Step 12 (EDA)
                                                   6 publication-quality figures
                                                   Comprehensive EDA notebook
```

---

## 6. How Files Connect to Each Other

### Join Keys

Every merge uses **`patient_id`** (TCGA-XXXX-XXXX format, 12 characters) as the join key.

| Merge | Left | Right | Join Type | Key |
|---|---|---|---|---|
| Clinical + Follow-up | clinical.tsv | follow_up.tsv | LEFT | `cases.submitter_id` |
| Clinical + Survival | patients_master.csv | survival.csv | LEFT | `patient_id` |
| Clinical + TMB | clinical+survival | mutation_sample_summary.csv | LEFT | `patient_id` |
| Clinical + Mutations | clinical+TMB | mutation_matrix_top200.parquet | LEFT | `patient_id` |
| Expression UUID mapping | STAR count files | gdc_sample_sheet.tsv | INNER | `File ID` (UUID) |
| Expression to patient | sample barcode | patient_id | Derived | `barcode[:12]` extracts patient_id |

### Why LEFT joins everywhere?

**Assumption**: The clinical table (1,098 patients) is the master anchor. LEFT join ensures:
- All 1,098 patients remain in the final cohort
- Patients without mutation data (130 patients) get `mut_* = 0` and `tmb = 0` (filled)
- Patients without survival data (24 patients) get `OS_time = NaN` and `OS_event = NaN` (preserved as-is)
- No patient is silently dropped due to missing omics data

### The Expression Data Linkage

Expression files don't directly join to the cohort table. The connection works as:

1. Each expression TSV lives in a UUID-named folder
2. `gdc_sample_sheet.tsv` maps `File ID (UUID)` to `Sample ID` (TCGA barcode)
3. The first 12 characters of the Sample ID = `patient_id`
4. The `has_rna_seq` flag in the cohort table indicates whether that patient has expression data
5. For actual expression analysis, the expression matrices are loaded separately and joined by barcode prefix

---

## 7. Summary of Key Design Decisions

| Decision | Rationale |
|---|---|
| Keep `age_at_diagnosis` in **days** (not convert to years) | GDC stores age in days; converting loses precision. The EDA notebook computes years where needed for display. |
| Median imputation for age only | Age is the only numeric clinical feature where imputation makes statistical sense. Race, ethnicity, and gender cannot be imputed. |
| Non-silent mutations only | Silent mutations have no functional impact; they add noise for downstream ML tasks. |
| Top 200 genes (not all ~18,000 mutated genes) | Sparsity: most genes are mutated in <1% of patients. Top 200 captures all biologically meaningful driver and passenger genes. |
| Top 2,000 variable genes (not all 60,660) | Reduces dimensionality for ML. Variance-based selection retains genes with the most discriminative signal. Protein-coding only removes pseudogenes and ncRNAs. |
| CPM + log2(x+1) normalization | CPM corrects for library size differences between samples. Log2 transform compresses the dynamic range and approximates normality — both are standard for RNA-Seq. |
| LEFT join from clinical as anchor | Ensures no patient is silently lost. Missing omics data is flagged (has_rna_seq) or zero-filled (mutations), never row-dropped. |
