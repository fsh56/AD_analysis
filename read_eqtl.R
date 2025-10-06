library(data.table)
library(tidyverse)

brain_tissues = c("Brain_Amygdala")
tissue_short = c("Amygdala")

for (chr in 1:2) {
  cat(paste0("process chr", chr, "...\n"))
  
  for (j in 1:length(brain_tissues)) {
    cat(paste0("process tissue:", brain_tissues[j], "\n"))
    
    eqtl_path = paste0("/gpfs/data/gao-lab/people/Sihao/data/GTEx_Analysis_v10_QTLs-GTEx_Analysis_v10_eQTL_all_associations-", 
                       brain_tissues[j], ".v10.allpairs.chr", chr, ".txt.gz")
    eqtl_bk = fread(eqtl_path) %>% 
      separate(gene_id, into = c("gene_name", "transcript"), sep = "\\.") %>% 
      dplyr::rename(b=slope, se=slope_se, p=pval_nominal) %>% 
      dplyr::filter(se!=Inf, se!=0, p!=1)
    
    fwrite(eqtl_bk, paste0("eqtl_Brain_Amygdala_", chr))
  }
}
cat("Done!\n")

library(data.table)
brain_tissues = c("Brain_Amygdala")
tissue_short = c("Amygdala")

for (chr in 1:2) {
  cat(paste0("process chr", chr, "...\n"))
  
  for (j in 1:length(brain_tissues)) {
    cat(paste0("process tissue:", brain_tissues[j], "\n"))
    eqtl_path = paste0("/gpfs/data/gao-lab/people/Sihao/data/GTEx_Analysis_v10_QTLs-GTEx_Analysis_v10_eQTL_all_associations-", 
                       brain_tissues[j], ".v10.allpairs.chr", chr, ".txt.gz")
  
    eqtl_bk = fread(eqtl_path)[
      , c("gene_name", "transcript") := tstrsplit(gene_id, "\\.", fixed = FALSE)][
        , .(gene_name, transcript, variant_id, tss_distance, ma_samples, ma_count, 
            maf, pval_nominal, b = slope, se = slope_se, pval_nominal_threshold, 
            min_pval_nominal, pval_beta, ref, alt, chr, pos)][
              !is.infinite(se) & se != 0 & pval_nominal != 1]
    cat("Processed rows:", nrow(eqtl_bk), "\n")
    fwrite(eqtl_bk, paste0("eqtl_Brain_Amygdala_", chr, ".txt.gz"))
    rm(eqtl_bk)
    gc()
  }
}

cat("Done!\n")