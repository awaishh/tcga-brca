import os
import re
import shutil

def main():
    print("Running 07_organize_expression.py")
    
    source_dir = '.'
    target_dir = 'data/raw/expression'
    
    os.makedirs(target_dir, exist_ok=True)
    
    uuid_pattern = re.compile(r'^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$')
    
    moved_count = 0
    
    for item in os.listdir(source_dir):
        if os.path.isdir(item) and uuid_pattern.match(item):
            src_path = os.path.join(source_dir, item)
            dst_path = os.path.join(target_dir, item)
            
            # Move the folder
            shutil.move(src_path, dst_path)
            moved_count += 1
            
    print(f"Found and moved {moved_count} UUID folders to {target_dir}")

if __name__ == "__main__":
    main()
