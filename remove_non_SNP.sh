cd /gpfs/data/gao-lab/people/Sihao/data/processed_gtex/

for file in *with_rsid.txt.gz; do
    echo "Processing: $file"
    newfile="${file/with_rsid.txt.gz/with_rsid_all_SNP.txt.gz}"
    zcat "$file" | awk -F'\t' 'NR==1 {print; next} {
        split($2, v, "_"); 
        ref=v[3]; 
        alt=v[4]; 
        if(length(ref)==1 && length(alt)==1) print
    }' | gzip > /gpfs/data/gao-lab/people/Sihao/data/gtex_with_rsid_all_SNP/"$newfile"
    
    echo "  -> Created: $newfile"
done
echo ""
echo "Done! All files saved to /gpfs/data/gao-lab/people/Sihao/data/gtex_with_rsid_all_SNP/"

