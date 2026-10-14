import pandas as pd
import numpy as np
import os

def clean_missing(df):
    missing_vals = ["'--", "Not Reported", "not reported", "Unknown"]
    return df.replace(missing_vals, np.nan)

def main():
    print("Running 02_patients_master.py")
    
    # Load data
    clinical = pd.read_csv("data/raw/clinical/clinical.tsv", sep='\t', dtype=str)
    follow_up = pd.read_csv("data/raw/clinical/follow_up.tsv", sep='\t', dtype=str)
    
    # Clean missing values
    clinical = clean_missing(clinical)
    follow_up = clean_missing(follow_up)
    
    # --- Deduplicate clinical ---
    print(f"Total clinical rows before dedup: {len(clinical)}")
    patients_before = clinical['cases.submitter_id'].nunique()
    print(f"Unique patients before dedup: {patients_before}")
    print(f"Average rows per patient: {len(clinical) / patients_before:.2f}")
    
    # Check 'diagnoses.diagnosis_is_primary_disease'
    print("\nUnique values in diagnoses.diagnosis_is_primary_disease:")
    print(clinical['diagnoses.diagnosis_is_primary_disease'].value_counts(dropna=False))
    
    # Filter for primary disease if the column exists and has values
    # Some datasets might use 'True', 'true', 'Yes', 'yes'
    if 'diagnoses.diagnosis_is_primary_disease' in clinical.columns:
        is_primary = clinical['diagnoses.diagnosis_is_primary_disease'].astype(str).str.lower().isin(['true', 'yes'])
        if is_primary.sum() > 0:
            clinical_filtered = clinical[is_primary].copy()
        else:
            clinical_filtered = clinical.copy()
    else:
        clinical_filtered = clinical.copy()
        
    # Still might have multiple rows (e.g. multiple treatments for the same primary diagnosis)
    clinical_dedup = clinical_filtered.drop_duplicates(subset=['cases.submitter_id'], keep='first')
    
    # If some patients were lost during filtering, fallback to keeping their first row from the original
    missing_patients = set(clinical['cases.submitter_id']) - set(clinical_dedup['cases.submitter_id'])
    if missing_patients:
        fallback = clinical[clinical['cases.submitter_id'].isin(missing_patients)].drop_duplicates(subset=['cases.submitter_id'], keep='first')
        clinical_dedup = pd.concat([clinical_dedup, fallback])
        
    print(f"\nTotal clinical rows after dedup: {len(clinical_dedup)}")
    print(f"Unique patients after dedup: {clinical_dedup['cases.submitter_id'].nunique()}")
    
    # --- Process follow_up ---
    print(f"\nTotal follow_up rows: {len(follow_up)}")
    print(f"Unique patients in follow_up: {follow_up['cases.submitter_id'].nunique()}")
    
    # Convert to numeric to find max
    follow_up['follow_ups.days_to_follow_up_num'] = pd.to_numeric(follow_up['follow_ups.days_to_follow_up'], errors='coerce')
    
    # Sort by days_to_follow_up and drop duplicates to keep the latest
    follow_up_dedup = follow_up.sort_values(by='follow_ups.days_to_follow_up_num', ascending=False, na_position='last')
    follow_up_dedup = follow_up_dedup.drop_duplicates(subset=['cases.submitter_id'], keep='first')
    
    # --- Merge ---
    merged = pd.merge(clinical_dedup, follow_up_dedup, on='cases.submitter_id', how='left')
    
    print(f"\nPatients in clinical but missing from follow_up: {merged['follow_ups.follow_up_id'].isna().sum()}")
    
    # --- Extract & Map Columns ---
    columns_map = {
        'cases.submitter_id': 'patient_id',
        'diagnoses.age_at_diagnosis': 'age_at_diagnosis',
        'demographic.sex_at_birth': 'gender',
        'demographic.race': 'race',
        'demographic.ethnicity': 'ethnicity',
        'diagnoses.ajcc_pathologic_stage': 'ajcc_pathologic_stage',
        'diagnoses.ajcc_pathologic_t': 'ajcc_pathologic_t',
        'diagnoses.ajcc_pathologic_n': 'ajcc_pathologic_n',
        'diagnoses.ajcc_pathologic_m': 'ajcc_pathologic_m',
        'diagnoses.primary_diagnosis': 'primary_diagnosis',
        'diagnoses.morphology': 'morphology',
        'demographic.vital_status': 'vital_status',
        'demographic.days_to_death': 'days_to_death',
        'diagnoses.days_to_last_follow_up': 'days_to_last_follow_up_clinical',
        'follow_ups.days_to_follow_up': 'days_to_last_follow_up'
    }
    
    # Ensure all columns exist before extracting
    for col in columns_map.keys():
        if col not in merged.columns:
            print(f"WARNING: Column {col} not found in merged data! Filling with NaN.")
            merged[col] = np.nan
            
    final_df = merged[list(columns_map.keys())].rename(columns=columns_map)
    
    os.makedirs("data/processed", exist_ok=True)
    out_path = "data/processed/patients_master.csv"
    final_df.to_csv(out_path, index=False)
    
    print(f"\nSaved {out_path}")
    print(f"Shape: {final_df.shape}")
    print("Head:")
    print(final_df.head())

if __name__ == "__main__":
    main()
