#!/bin/bash
# 基础路径
base_dir="/gpfs/data/gao-lab/people/Sihao/amyloid_ba9"
out_dir="${base_dir}/gene_info"
out_file="${out_dir}/BA9_gene_summary_w50kb_p0.001_all.txt"

# 逐个合并文件
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

echo "✅ Merge complete! Output saved to:"
echo "$out_file"
