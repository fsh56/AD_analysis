import pandas as pd
import glob
import os

# Set directories
input_dir = "/gpfs/data/gao-lab/gtex_v10_eqtl"
output_dir = "/gpfs/data/gao-lab/people/Sihao/data"

# Process all parquet files
for filepath in glob.glob(os.path.join(input_dir, "*.parquet")):
    filename = os.path.basename(filepath)
    outname = filename.replace(".parquet", ".txt.gz")
    outpath = os.path.join(output_dir, outname)
    print(f"Processing {filename} to {outname}")
    df = pd.read_parquet(filepath)
    df.to_csv(outpath, sep="\t", index=False, compression="gzip")



import pandas as pd
import glob
import os
import time

# Set directories
input_dir = "/gpfs/data/gao-lab/gtex_v10_eqtl"
output_dir = "/gpfs/data/gao-lab/people/Sihao/data"

# Get all parquet files
files = glob.glob(os.path.join(input_dir, "*.parquet"))
print(f"Found {len(files)} parquet files to process")

# Process all parquet files
for i, filepath in enumerate(files, 1):
    filename = os.path.basename(filepath)
    outname = filename.replace(".parquet", ".txt.gz")
    outpath = os.path.join(output_dir, outname)
    
    print(f"[{i}/{len(files)}] Processing {filename}...")
    start_time = time.time()
    
    df = pd.read_parquet(filepath)
    df.to_csv(outpath, sep="\t", index=False, compression="gzip")
    
    elapsed = time.time() - start_time
    print(f"Completed in {elapsed:.1f} seconds")

print("All files processed successfully!")