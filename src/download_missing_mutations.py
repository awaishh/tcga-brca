import os
import glob
import time
import requests

def download_missing_mutations():
    manifest_path = 'data/raw/mutations/gdc_manifest.2026-09-21.101715.txt'
    mutations_dir = 'data/raw/mutations'
    
    with open(manifest_path, 'r') as f:
        lines = [line.strip().split('\t') for line in f.readlines() if line.strip()]
    
    header = lines[0]
    records = lines[1:]
    
    print(f"Total manifest records: {len(records)}")
    
    missing = []
    for r in records:
        uuid, filename, md5, size, state = r[0], r[1], r[2], r[3], r[4]
        target_dir = os.path.join(mutations_dir, uuid)
        target_file = os.path.join(target_dir, filename)
        
        # Check if file exists and has size > 0
        if not (os.path.exists(target_file) and os.path.getsize(target_file) > 0):
            missing.append((uuid, filename, int(size)))
            
    print(f"Total missing or incomplete files: {len(missing)}")
    if not missing:
        print("All files are already present and verified!")
        return

    session = requests.Session()
    adapter = requests.adapters.HTTPAdapter(max_retries=5)
    session.mount('https://', adapter)
    session.mount('http://', adapter)
    
    success = 0
    failed = []
    
    for idx, (uuid, filename, expected_size) in enumerate(missing, 1):
        url = f"https://api.gdc.cancer.gov/data/{uuid}"
        target_dir = os.path.join(mutations_dir, uuid)
        os.makedirs(target_dir, exist_ok=True)
        target_file = os.path.join(target_dir, filename)
        
        print(f"[{idx}/{len(missing)}] Downloading {uuid}/{filename} ({expected_size/1024:.1f} KB)...", end=" ", flush=True)
        
        downloaded = False
        for attempt in range(5):
            try:
                response = session.get(url, stream=True, timeout=60)
                if response.status_code == 200:
                    with open(target_file, 'wb') as out_f:
                        for chunk in response.iter_content(chunk_size=1024*64):
                            if chunk:
                                out_f.write(chunk)
                    
                    actual_size = os.path.getsize(target_file)
                    if actual_size == expected_size:
                        print(f"Done! ({actual_size} bytes)")
                        success += 1
                        downloaded = True
                        break
                    else:
                        print(f"Size mismatch: got {actual_size}, expected {expected_size}. Retrying...")
                else:
                    print(f"HTTP {response.status_code}. Retrying...")
            except Exception as e:
                print(f"Error ({e}). Retrying in 2s...")
                time.sleep(2)
        
        if not downloaded:
            print("FAILED.")
            failed.append(uuid)
            
    print("\n==========================================")
    print(f"Download complete: {success} succeeded, {len(failed)} failed.")
    print("==========================================")
    
    total_present = len(glob.glob(os.path.join(mutations_dir, '*', '*.maf.gz')))
    print(f"Total MAF files on disk: {total_present} / {len(records)}")

if __name__ == '__main__':
    download_missing_mutations()
