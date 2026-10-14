import pandas as pd
import os
import glob

def generate_inventory():
    print("Running 01_data_inventory.py")
    base_dir = "data/raw"
    report_path = "reports/data_inventory.md"
    
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    
    files = glob.glob(f"{base_dir}/**/*.tsv", recursive=True)
    
    with open(report_path, "w") as f:
        f.write("# Data Inventory\n\n")
        
        for file_path in sorted(files):
            print(f"\nProcessing file: {file_path}")
            # Use dtype=str to avoid mixed type warnings on large unstructured datasets
            df = pd.read_csv(file_path, sep='\t', dtype=str)
            
            # Print to stdout per PRD
            print(f"Shape: {df.shape}")
            print("Head:")
            print(df.head())
            
            # Write to markdown
            f.write(f"## File: `{file_path}`\n")
            f.write(f"- **Row count:** {df.shape[0]}\n")
            f.write(f"- **Column count:** {df.shape[1]}\n\n")
            
            f.write("### Columns\n")
            f.write(", ".join([f"`{col}`" for col in df.columns]) + "\n\n")
            
            f.write("### First 3 Rows\n")
            f.write(df.head(3).to_markdown(index=False))
            f.write("\n\n---\n\n")

    print(f"\nData inventory generated at {report_path}")

if __name__ == "__main__":
    generate_inventory()
