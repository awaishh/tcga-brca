import os
import glob
import pandas as pd

def main():
    print("Running 08_build_expression_matrix.py")
    
    # 1. Read sample sheet to create mapping
    # Look for the sample sheet dynamically if the date changes
    sample_sheets = glob.glob('gdc_sample_sheet*.tsv')
    if not sample_sheets:
        raise FileNotFoundError("Could not find gdc_sample_sheet*.tsv in the project root.")
    
    sample_sheet_path = sample_sheets[0]
    sheet_df = pd.read_csv(sample_sheet_path, sep='\t')
    
    # Create mapping: File ID (UUID) -> Sample ID
    uuid_to_sample = dict(zip(sheet_df['File ID'], sheet_df['Sample ID']))
    
    # 2. Iterate through folders in data/raw/expression/
    expression_dir = 'data/raw/expression'
    
    if not os.path.exists(expression_dir):
        raise FileNotFoundError(f"{expression_dir} does not exist.")
        
    on_disk_uuids = [d for d in os.listdir(expression_dir) if os.path.isdir(os.path.join(expression_dir, d))]
    
    # Count how many in sample sheet are missing on disk
    sheet_uuids = set(uuid_to_sample.keys())
    missing_on_disk = sheet_uuids - set(on_disk_uuids)
    print(f"Total entries in sample sheet: {len(sheet_uuids)}")
    print(f"Skipping {len(missing_on_disk)} files listed in sample sheet but missing on disk.")
    
    # Process files
    series_list = []
    processed_count = 0
    sample_counts = {}
    
    # Sort for deterministic processing order
    for uuid in sorted(on_disk_uuids):
        if uuid not in uuid_to_sample:
            # If there's a folder not in sample sheet, warn but skip
            print(f"WARNING: UUID {uuid} found on disk but not in sample sheet. Skipping.")
            continue
            
        sample_id = uuid_to_sample[uuid]
        folder_path = os.path.join(expression_dir, uuid)
        
        # Find the TSV file inside
        tsv_files = glob.glob(os.path.join(folder_path, '*.rna_seq.augmented_star_gene_counts.tsv'))
        if not tsv_files:
            print(f"WARNING: No counts TSV found in {folder_path}. Skipping.")
            continue
            
        tsv_path = tsv_files[0]
        
        # Read TSV
        # The first line is a comment starting with '#'
        df = pd.read_csv(tsv_path, sep='\t', comment='#')
        
        # Drop summary rows starting with 'N_'
        df = df[~df['gene_id'].str.startswith('N_')]
        
        # Deduplicate sample_id if multiple files exist for the same sample
        sample_counts[sample_id] = sample_counts.get(sample_id, 0) + 1
        if sample_counts[sample_id] > 1:
            unique_sample_id = f"{sample_id}_rep{sample_counts[sample_id]}"
        else:
            unique_sample_id = sample_id
            
        # Extract series and set name to unique_sample_id
        # We assume gene_id is unique per file
        s = df.set_index('gene_id')['unstranded']
        s.name = unique_sample_id
        
        series_list.append(s)
        processed_count += 1
        
        if processed_count % 25 == 0:
            print(f"Processed {processed_count} files...")
            
    print(f"Finished reading {processed_count} files. Concatenating into matrix...")
    
    # 3. Concatenate all series into a single dataframe
    # This aligns them by the index (gene_id)
    matrix_df = pd.concat(series_list, axis=1)
    
    # 4. Save to parquet
    os.makedirs('data/processed', exist_ok=True)
    out_path = 'data/processed/expression_matrix.parquet'
    
    # We can reset index so gene_id becomes a regular column, or just leave it as index.
    # Leaving it as index is fine for parquet, but standardizing to a column is safer for interoperability.
    matrix_df = matrix_df.reset_index()
    matrix_df.to_parquet(out_path, index=False)
    
    print(f"\nSaved expression matrix to {out_path}")
    print(f"Final shape: {matrix_df.shape}")
    print(f"Expected ~60000 genes x {processed_count + 1} columns (gene_id + samples)")
    
if __name__ == "__main__":
    main()
