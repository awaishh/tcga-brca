import pandas as pd
import os

descriptions = {
    'patient_id': 'Unique identifier for the patient (TCGA case submitter ID).',
    'age_at_diagnosis': 'Age of the patient in days at the time of primary diagnosis.',
    'gender': 'Sex of the patient at birth.',
    'race': 'Race of the patient.',
    'ethnicity': 'Ethnicity of the patient.',
    'ajcc_pathologic_stage': 'AJCC pathologic stage of the tumor.',
    'ajcc_pathologic_t': 'AJCC pathologic T (tumor) stage.',
    'ajcc_pathologic_n': 'AJCC pathologic N (node) stage.',
    'ajcc_pathologic_m': 'AJCC pathologic M (metastasis) stage.',
    'primary_diagnosis': 'Text description of the primary diagnosis.',
    'morphology': 'ICD-O-3 morphology code.',
    'vital_status': 'Survival state of the patient (Alive/Dead).',
    'days_to_death': 'Number of days from diagnosis to death.',
    'days_to_last_follow_up_clinical': 'Number of days to last follow-up recorded in clinical data.',
    'days_to_last_follow_up': 'Number of days to latest follow-up recorded in follow_up data.',
    'OS_time': 'Overall survival time in days (days to death or last follow-up).',
    'OS_event': 'Overall survival event indicator (1 = Dead, 0 = Alive).',
    'sample_barcode': 'Unique identifier for the biological sample.',
    'sample_type_code': '2-digit sample type code extracted from the barcode.',
    'sample_type': 'Text description of the sample type (e.g., Primary Tumor, Blood Derived Normal).'
}

def main():
    print("Running 06_data_dictionary.py")
    
    files = {
        'patients_master.csv': 'data/processed/patients_master.csv',
        'survival.csv': 'data/processed/survival.csv',
        'samples.csv': 'data/processed/samples.csv'
    }
    
    dict_rows = []
    
    for filename, filepath in files.items():
        if os.path.exists(filepath):
            df = pd.read_csv(filepath)
            for col in df.columns:
                dtype_str = str(df[col].dtype)
                pct_missing = (df[col].isna().sum() / len(df)) * 100
                desc = descriptions.get(col, 'Description not available.')
                
                dict_rows.append({
                    'Column Name': col,
                    'Source File': filename,
                    'Data Type': dtype_str,
                    '% Missing': f"{pct_missing:.2f}%",
                    'Description': desc
                })
                
    dict_df = pd.DataFrame(dict_rows)
    
    report_path = 'reports/data_dictionary.md'
    with open(report_path, 'w') as f:
        f.write("# Data Dictionary\n\n")
        f.write("This dictionary describes the columns present in the processed dataset.\n\n")
        f.write(dict_df.to_markdown(index=False))
        
    print(f"Data dictionary generated at {report_path}")

if __name__ == "__main__":
    main()
