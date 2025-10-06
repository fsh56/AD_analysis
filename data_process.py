import pandas as pd
import glob
import os
from pathlib import Path
import concurrent.futures
from datetime import datetime

input_dir = "/gpfs/data/gao-lab/gtex_v10_eqtl"
output_dir = "/gpfs/data/gao-lab/people/Sihao/data"

def process_single_parquet(filepath, output_dir):
    filename = os.path.basename(filepath)
    
    try:
        tissue_name = filename.replace('.parquet', '').replace('GTEx_v10_eQTL_', '')
        print(f"[{datetime.now().strftime('%H:%M:%S')}] Processing {tissue_name}...")
        df = pd.read_parquet(filepath)
        initial_rows = len(df)
        
        # ========================================
        # ADD YOUR PROCESSING LOGIC HERE
        # ========================================
        
        # Example processing - modify based on your needs:
        
        # 1. Filter by p-value if column exists
        if 'pval_nominal' in df.columns:
            df = df[df['pval_nominal'] < 0.05]  # Adjust threshold as needed
        
        # 2. Filter by MAF (Minor Allele Frequency) if exists
        if 'maf' in df.columns:
            df = df[df['maf'] >= 0.01]  # Keep variants with MAF >= 1%
        
        # 3. Add computed columns if needed
        if 'slope' in df.columns and 'slope_se' in df.columns:
            df['effect_size_zscore'] = df['slope'] / df['slope_se']
        
        # 4. Sort by significance if p-value exists
        if 'pval_nominal' in df.columns:
            df = df.sort_values('pval_nominal')
        
        # 5. Remove unnecessary columns if needed (uncomment and modify)
        # columns_to_drop = ['unnecessary_col1', 'unnecessary_col2']
        # df = df.drop(columns=[col for col in columns_to_drop if col in df.columns])
        
        # ========================================
        # END OF PROCESSING LOGIC
        # ========================================
        
        final_rows = len(df)
        
        # Create tissue-specific output directory
        tissue_output_dir = os.path.join(output_dir, tissue_name)
        Path(tissue_output_dir).mkdir(parents=True, exist_ok=True)
        
        # Save processed data as parquet
        output_path = os.path.join(tissue_output_dir, f"{tissue_name}_processed.parquet")
        df.to_parquet(output_path, compression='snappy', index=False)
        
        # Also save a summary statistics file
        stats = {
            'tissue': tissue_name,
            'original_file': filename,
            'initial_rows': initial_rows,
            'final_rows': final_rows,
            'rows_filtered': initial_rows - final_rows,
            'filter_percentage': ((initial_rows - final_rows) / initial_rows * 100) if initial_rows > 0 else 0,
            'columns': list(df.columns),
            'n_columns': len(df.columns),
            'file_size_mb': os.path.getsize(output_path) / (1024 * 1024)
        }
        
        stats_path = os.path.join(tissue_output_dir, f"{tissue_name}_processing_stats.txt")
        with open(stats_path, 'w') as f:
            for key, value in stats.items():
                if key != 'columns':
                    f.write(f"{key}: {value}\n")
                else:
                    f.write(f"{key}: {', '.join(value[:5])}...\n")  # Show first 5 columns
        
        return f"✓ {tissue_name}: {initial_rows:,} → {final_rows:,} rows ({stats['filter_percentage']:.1f}% filtered)"
        
    except Exception as e:
        return f"✗ {filename}: ERROR - {str(e)}"

def process_all_tissues(input_dir=INPUT_DIR, output_dir=OUTPUT_DIR, 
                       max_workers=4, tissue_filter=None):
    """
    Process all parquet files in the input directory
    
    Args:
        input_dir: Directory containing input parquet files
        output_dir: Directory to save processed files
        max_workers: Number of parallel workers
        tissue_filter: Optional list of tissue names to process (None = process all)
    """
    
    print(f"=" * 60)
    print(f"GTEx v10 Parquet Processing Script")
    print(f"=" * 60)
    print(f"Input directory:  {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Max workers: {max_workers}")
    print(f"=" * 60)
    
    # Create main output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Get all parquet files
    all_files = sorted(glob.glob(os.path.join(input_dir, "*.parquet")))
    
    if not all_files:
        print(f"No parquet files found in {input_dir}")
        return
    
    # Filter by specific tissues if requested
    if tissue_filter:
        files_to_process = []
        for tissue in tissue_filter:
            matching_files = [f for f in all_files if tissue in f]
            files_to_process.extend(matching_files)
        files_to_process = list(set(files_to_process))  # Remove duplicates
    else:
        files_to_process = all_files
    
    print(f"Found {len(files_to_process)} files to process")
    print(f"Starting processing at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"-" * 60)
    
    # Process files in parallel
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_file = {
            executor.submit(process_single_parquet, filepath, output_dir): filepath 
            for filepath in files_to_process
        }
        
        for future in concurrent.futures.as_completed(future_to_file):
            result = future.result()
            results.append(result)
            print(result)
    
    # Summary
    print(f"-" * 60)
    print(f"Processing completed at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    successful = sum(1 for r in results if r.startswith('✓'))
    failed = sum(1 for r in results if r.startswith('✗'))
    
    print(f"Successfully processed: {successful}/{len(files_to_process)} files")
    if failed > 0:
        print(f"Failed: {failed} files")
    
    # Create overall summary file
    summary_path = os.path.join(output_dir, "processing_summary.txt")
    with open(summary_path, 'w') as f:
        f.write(f"GTEx v10 Processing Summary\n")
        f.write(f"{'=' * 40}\n")
        f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Input directory: {input_dir}\n")
        f.write(f"Output directory: {output_dir}\n")
        f.write(f"Total files processed: {len(files_to_process)}\n")
        f.write(f"Successful: {successful}\n")
        f.write(f"Failed: {failed}\n")
        f.write(f"\nDetailed Results:\n")
        f.write(f"{'-' * 40}\n")
        for result in sorted(results):
            f.write(f"{result}\n")
    
    print(f"Summary saved to: {summary_path}")

def main():
    """Main execution function"""
    
    # Option 1: Process ALL tissues
    process_all_tissues(
        input_dir=INPUT_DIR,
        output_dir=OUTPUT_DIR,
        max_workers=4  # Adjust based on your system capacity
    )
    
    # Option 2: Process specific tissues only (uncomment to use)
    # tissues_to_process = [
    #     "Adipose_Subcutaneous",
    #     "Adipose_Visceral_Omentum",
    #     "Brain_Cortex",
    #     "Brain_Hippocampus",
    #     "Muscle_Skeletal",
    #     "Whole_Blood",
    #     "Liver",
    #     "Heart_Left_Ventricle"
    # ]
    # 
    # process_all_tissues(
    #     input_dir=INPUT_DIR,
    #     output_dir=OUTPUT_DIR,
    #     max_workers=4,
    #     tissue_filter=tissues_to_process
    # )

if __name__ == "__main__":
    main()xq