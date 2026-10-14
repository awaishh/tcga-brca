import pandas as pd
import os

def main():
    print("Running 05_excel_export.py")
    
    files = {
        'patients_master': 'data/processed/patients_master.csv',
        'survival': 'data/processed/survival.csv',
        'samples': 'data/processed/samples.csv'
    }
    
    out_path = 'data/processed/tcga_brca.xlsx'
    
    with pd.ExcelWriter(out_path, engine='openpyxl') as writer:
        for sheet_name, file_path in files.items():
            if os.path.exists(file_path):
                df = pd.read_csv(file_path)
                df.to_excel(writer, sheet_name=sheet_name, index=False)
                print(f"Added {sheet_name} with shape {df.shape}")
            else:
                print(f"WARNING: {file_path} not found. Skipping.")
                
    print(f"\nSuccessfully created {out_path}")

if __name__ == "__main__":
    main()
