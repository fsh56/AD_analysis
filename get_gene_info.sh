#!/bin/bash

# 基础路径
base_dir="/gpfs/data/gao-lab/people/Sihao/amyloid_amydgala"
out_dir="${base_dir}/gene_info"
out_file="${out_dir}/amy_gene_summary_w50kb_p0.001_all.txt"

mkdir -p "$out_dir"

echo "=== Step 1: Adding chromosome column to each file ==="
# 遍历 1 到 22 号染色体，添加 chr 列
for chr in {1..22}; do
    infile="${base_dir}/ldByGene${chr}/BA9_gene_summary_w50kb_p0.001.txt"
    outfile="${base_dir}/ldByGene${chr}/BA9_gene_summary_w50kb_p0.001_chr.txt"
    
    if [[ -f "$infile" ]]; then
        echo "Processing chromosome ${chr}..."
        awk -v chr="$chr" 'NR==1 {print $0"\tchr"} NR>1 {print $0"\t"chr}' "$infile" > "$outfile"
    else
        echo "Warning: $infile not found, skipping."
    fi
done

echo ""
echo "=== Step 2: Merging all chromosome files ==="
# 合并所有染色体文件
for chr in {1..22}; do
    infile="${base_dir}/ldByGene${chr}/BA9_gene_summary_w50kb_p0.001_chr.txt"
    
    if [[ -f "$infile" ]]; then
        echo "Merging chromosome ${chr}..."
        if [[ $chr -eq 1 ]]; then
            # 第一个文件保留表头
            cat "$infile" > "$out_file"
        else
            # 之后的文件去掉表头再追加
            tail -n +2 "$infile" >> "$out_file"
        fi
    else
        echo "Warning: $infile not found, skipping."
    fi
done

echo ""
echo "✅ All done!"
echo "Files with chr column: ${base_dir}/ldByGene*/amy_gene_summary_w50kb_p0.001_chr.txt"
echo "Merged file saved to: $out_file"