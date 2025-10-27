import pandas as pd
import glob
import os
import time

# Set directories
input_dir = "/gpfs/data/gao-lab/gtex_v10_eqtl"
output_dir = "/gpfs/data/gao-lab/people/Sihao/data"

# Create output directory if needed
os.makedirs(output_dir, exist_ok=True)

# Get all brain-related parquet files
files = glob.glob(os.path.join(input_dir, "*Brain_*.parquet"))
print(f"Found {len(files)} brain parquet files to process")

# Process all files
for i, filepath in enumerate(files, 1):
    filename = os.path.basename(filepath)
    outname = filename.replace(".parquet", ".txt.gz")
    outpath = os.path.join(output_dir, outname)
    
    # Read parquet and save as txt.gz
    df = pd.read_parquet(filepath)
    df.to_csv(outpath, sep="\t", index=False, compression="gzip")
    
    # Print progress every 10 files
    if i % 10 == 0:
        print(f"Processed {i}/{len(files)} files")

print(f"All {len(files)} files completed!")