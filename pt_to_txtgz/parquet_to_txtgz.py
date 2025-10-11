import pandas as pd
import glob
import os

# Set the directory containing the Parquet files
input_dir = "/gpfs/data/gao-lab/gtex_v10_eqtl"
output_dir = "/gpfs/data/gao-lab/people/Sihao/ad_analysis"

# Loop over Parquet files (files 8-22)
for filepath in glob.glob(os.path.join(input_dir, "*Adipose_Subcutaneous*.parquet"))[7:22]:
    filename = os.path.basename(filepath)
    outname = filename.replace(".parquet", ".txt.gz")
    outpath = os.path.join(output_dir, outname)
    print(f"Processing {filename} → {outname}")
    df = pd.read_parquet(filepath)
    df.to_csv(outpath, sep="\t", index=False, compression="gzip")


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