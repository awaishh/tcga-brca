import pandas as pd

def main():
    print("Running 04_samples.py")
    
    # Load data
    df = pd.read_csv("data/raw/biospecimen/sample.tsv", sep='\t', dtype=str)
    
    # Required columns from raw
    # patient_id -> cases.submitter_id
    # sample_barcode -> samples.submitter_id
    
    # Derive sample_type_code
    # TCGA barcodes: e.g. TCGA-A2-A04Q-10A-01D-A128-09
    # Split by '-', the 4th element (index 3) is '10A'
    # The first 2 characters are the sample type code (e.g., '10')
    df['sample_barcode'] = df['samples.submitter_id']
    df['patient_id'] = df['cases.submitter_id']
    
    # Extract the 4th group (index 3)
    def extract_code(barcode):
        if pd.isna(barcode):
            return barcode
        parts = barcode.split('-')
        if len(parts) >= 4:
            return parts[3][:2]
        return barcode
        
    df['sample_type_code'] = df['sample_barcode'].apply(extract_code)
    
    # Map sample_type
    code_map = {
        '01': 'Primary Tumor',
        '06': 'Metastatic',
        '10': 'Blood Derived Normal',
        '11': 'Solid Tissue Normal'
    }
    
    # If not in map, keep raw code
    df['sample_type'] = df['sample_type_code'].map(code_map).fillna(df['sample_type_code'])
    
    # Subset
    samples_df = df[['patient_id', 'sample_barcode', 'sample_type_code', 'sample_type']].copy()
    
    # There could be duplicates if the raw sample.tsv contains duplicated sample rows?
    samples_df = samples_df.drop_duplicates()
    
    out_path = "data/processed/samples.csv"
    samples_df.to_csv(out_path, index=False)
    
    print(f"Saved {out_path}")
    print(f"Shape: {samples_df.shape}")
    print("Head:")
    print(samples_df.head())
    print("\nsample_type_code distribution:")
    print(samples_df['sample_type_code'].value_counts())
    print("\nsample_type distribution:")
    print(samples_df['sample_type'].value_counts())

if __name__ == "__main__":
    main()
