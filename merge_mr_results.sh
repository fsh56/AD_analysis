#!/bin/bash
# ==========================================================
# Merge MR result files (Egger / IVW / WeightedMedian)
# ==========================================================

base_dir="/gpfs/data/gao-lab/people/Sihao/amyloid_amydgala"
out_dir="${base_dir}/mr_results1027"
methods=("Egger" "IVW" "WeightedMedian")
mkdir -p "$out_dir"

for method in "${methods[@]}"; do
    out_file="${out_dir}/amyloid_amy_${method}_results_all.csv"
    rm -f "$out_file"

    echo "========== Merging ${method} results =========="

    first_file=1
    for chr in {1..22}; do
        infile="${base_dir}/ldByGene${chr}/BA9_chr${chr}_${method}_results.csv"

        if [[ -f "$infile" ]]; then
            echo "Adding chr=${chr} from $infile"

            if [[ $first_file -eq 1 ]]; then
                # 第一个文件：保留表头，并在首列添加“chr”
                awk -v c="$chr" 'NR==1{print "chr," $0} NR>1{print c "," $0}' "$infile" > "$out_file"
                first_file=0
            else
                # 后续文件：跳过表头，只加数据
                awk -v c="$chr" 'NR>1{print c "," $0}' "$infile" >> "$out_file"
            fi
        else
            echo "⚠️ Warning: $infile not found, skipping."
        fi
    done

    echo "Done for ${method}: ${out_file}"
    echo
done

echo "All MR result files merged successfully with 'chr' column added!"
