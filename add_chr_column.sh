#!/bin/bash
# 基础路径
base_dir="/gpfs/data/gao-lab/people/Sihao/amyloid_ba9"

# 遍历 1 到 22 号染色体
for chr in {1..22}; do
    # 输入文件路径
    infile="${base_dir}/ldByGene${chr}/BA9_gene_summary_w50kb_p0.001.txt"
    # 输出文件路径（添加了新列）
    outfile="${base_dir}/ldByGene${chr}/BA9_gene_summary_w50kb_p0.001_chr.txt"

    # 确保输入文件存在
    if [[ -f "$infile" ]]; then
        echo "Processing chromosome ${chr}..."
        # 用 awk 添加新列 chr
        awk -v chr="$chr" 'NR==1 {print $0"\tchr"} NR>1 {print $0"\t"chr}' "$infile" > "$outfile"
    else
        echo "Warning: $infile not found, skipping."
    fi
done

echo "✅ Done! Files with '_chr.txt' suffix contain the new column."
