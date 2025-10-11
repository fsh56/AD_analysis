# check all
for f in /gpfs/data/gao-lab/people/Sihao/data/processed_gtex/*.txt.gz; do 
    echo "File: $(basename $f)"; 
    total=$(zcat "$f" | tail -n +2 | wc -l); 
    non_snp=$(zcat "$f" | tail -n +2 | cut -f2 | awk -F'_' '{ref=$(NF-2); alt=$(NF-1); if(length(ref)>1 || length(alt)>1) print}' | wc -l); 
    pct=$(awk "BEGIN {printf \"%.2f\", ($non_snp/$total)*100}"); 
    printf "%-50s Total: %8d  Non-SNPs: %8d (%.2f%%)\n" "$f" $total $non_snp $pct; 
done


# for specific chr
for f in /gpfs/data/gao-lab/people/Sihao/data/processed_gtex/*chr22_with_rsid.txt.gz; do 
    echo "File: $(basename $f)"; 
    total=$(zcat "$f" | tail -n +2 | wc -l); 
    non_snp=$(zcat "$f" | tail -n +2 | cut -f2 | awk -F'_' '{ref=$(NF-2); alt=$(NF-1); if(length(ref)>1 || length(alt)>1) print}' | wc -l); 
    pct=$(awk "BEGIN {printf \"%.2f\", ($non_snp/$total)*100}"); 
    printf "%-50s Total: %8d  Non-SNPs: %8d (%.2f%%)\n" "$f" $total $non_snp $pct; 
done

# check if there is assigned rsid yesss ---------------------------
zcat /gpfs/data/gao-lab/people/Sihao/data/processed_gtex/Brain_Amygdala_chr1_with_rsid.txt.gz | tail -n +2 | awk -F'\t' '{
    split($2, v, "_"); 
    ref=v[length(v)-2]; 
    alt=v[length(v)-1]; 
    if(length(ref)>1 || length(alt)>1) {
        if($NF != "." && $NF != "") print
    }
}' | wc -l

# check is reference table has rsid for nonSNP genetic variant
zcat /gpfs/data/gao-lab/people/Sihao/data/GTEx_Analysis_v10_QTLs-GTEx_Analysis_v10_eQTL_all_associations-Brain_Amygdala.v10.allpairs.chr1.txt.gz | tail -n +2 | awk -F'\t' '{
    split($2, v, "_"); 
    ref=v[length(v)-2]; 
    alt=v[length(v)-1]; 
    if(length(ref)>1 || length(alt)>1) print $2
}' | sort -u > /tmp/nonsnp_variants.txt
awk 'NR==FNR{a[$1]; next} $1 in a' /tmp/nonsnp_variants.txt /gpfs/data/gao-lab/people/Sihao/draft/GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt | wc -l
#56298

# 1. check non-unique non-SNP variants
total=$(wc -l < /tmp/nonsnp_variants.txt)
echo "Total unique non-SNPs: $total"
echo "Non-SNPs in lookup table: 56298"
echo "Coverage: $(awk "BEGIN {printf \"%.2f%%\", (56298/$total)*100}")"

# 2. check some eg
awk 'NR==FNR{a[$1]; next} $1 in a' /tmp/nonsnp_variants.txt /gpfs/data/gao-lab/people/Sihao/draft/GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt | head -20
awk 'NR==FNR{a[$1]; next} $1 in a && $2 != "."' /tmp/nonsnp_variants.txt /gpfs/data/gao-lab/people/Sihao/draft/GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt | wc -l
