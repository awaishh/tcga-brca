import pandas as pd
import numpy as np

def main():
    print("Running 03_survival.py")
    
    # Load data
    df = pd.read_csv("data/processed/patients_master.csv")
    
    initial_count = len(df)
    
    # 1. Create OS_event
    # Map Dead -> 1, Alive -> 0. Anything else (NaN) -> NaN
    event_map = {'Dead': 1, 'Alive': 0}
    df['OS_event'] = df['vital_status'].map(event_map)
    
    # 2. Create OS_time
    df['days_to_death'] = pd.to_numeric(df['days_to_death'], errors='coerce')
    df['days_to_last_follow_up'] = pd.to_numeric(df['days_to_last_follow_up'], errors='coerce')
    
    df['OS_time'] = np.where(df['vital_status'] == 'Dead', 
                             df['days_to_death'], 
                             df['days_to_last_follow_up'])
    
    # Subset to required columns
    surv_df = df[['patient_id', 'OS_time', 'OS_event']].copy()
    
    # 3. Identify rows to drop
    is_null = surv_df['OS_time'].isna()
    is_negative = surv_df['OS_time'] < 0
    is_zero = surv_df['OS_time'] == 0
    
    null_count = is_null.sum()
    negative_count = is_negative.sum()
    zero_count = is_zero.sum()
    
    print("\nDropped rows summary:")
    print(f"- Null OS_time: {null_count}")
    print(f"- Negative OS_time: {negative_count}")
    print(f"- Zero OS_time: {zero_count}")
    
    # Keep only valid rows
    valid_mask = ~(is_null | is_negative | is_zero)
    surv_df_clean = surv_df[valid_mask].copy()
    
    total_dropped = initial_count - len(surv_df_clean)
    print(f"- Total dropped: {total_dropped}")
    
    # 4. Save
    out_path = "data/processed/survival.csv"
    surv_df_clean.to_csv(out_path, index=False)
    
    print(f"\nSaved {out_path}")
    print(f"Shape: {surv_df_clean.shape}")
    print("Head:")
    print(surv_df_clean.head())

if __name__ == "__main__":
    main()
