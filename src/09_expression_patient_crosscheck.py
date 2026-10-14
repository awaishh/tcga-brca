import pandas as pd

def main():
    print("Running 09_expression_patient_crosscheck.py")
    
    # 1. Get columns from expression matrix
    # Using pyarrow directly or pandas. We just read 1 row to be fast and get columns.
    matrix_path = 'data/processed/expression_matrix.parquet'
    
    # Read just the first row to get columns quickly
    df_expr = pd.read_parquet(matrix_path)
    
    # Remove 'gene_id' from columns to get the list of samples
    sample_cols = [c for c in df_expr.columns if c != 'gene_id']
    
    print(f"Found {len(sample_cols)} expression samples in matrix.")
    
    # Extract patient IDs from sample columns
    # TCGA barcodes start with the project, TSS, and participant, making up 12 characters
    expr_patient_ids = [s[:12] for s in sample_cols]
    
    # 2. Load patients_master
    patients_path = 'data/processed/patients_master.csv'
    df_patients = pd.read_csv(patients_path)
    master_patient_ids = set(df_patients['patient_id'])
    
    print(f"Found {len(master_patient_ids)} patients in patients_master.csv.")
    
    # 3. Crosscheck
    matching_samples = 0
    non_matching_samples = []
    
    for i, s in enumerate(sample_cols):
        p_id = expr_patient_ids[i]
        if p_id in master_patient_ids:
            matching_samples += 1
        else:
            non_matching_samples.append((s, p_id))
            
    print(f"\nExpression samples with a matching patient in patients_master: {matching_samples}")
    print(f"Expression samples WITHOUT a matching patient: {len(non_matching_samples)}")
    
    if non_matching_samples:
        print("Examples of non-matching samples:")
        for s, p_id in non_matching_samples[:3]:
            print(f"  Sample: {s} -> Extracted Patient ID: {p_id}")
            
    # Patients missing expression data
    # Create a set of all unique patient IDs from expression samples
    expr_unique_patients = set(expr_patient_ids)
    missing_expr_patients = master_patient_ids - expr_unique_patients
    
    print(f"\nPatients in patients_master with NO expression data: {len(missing_expr_patients)}")

if __name__ == "__main__":
    main()
