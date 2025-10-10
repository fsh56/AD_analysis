library(data.table)
library(tidyverse)
library(ieugwasr)

clump <- function(dat_path, 
                  SNP_col = "variant_id", 
                  pval_col = "pval_nominal", 
                  clump_kb = 250, 
                  clump_r2 = 0.1, 
                  clump_p = 0.01, 
                  bfile = "/gpfs/data/linchen-lab/Yihao/Ke/education_AD_GWAS/EUR", 
                  plink_bin = genetics.binaRies::get_plink_binary(), 
                  pop="EUR") {
  dat <- fread(dat_path)
  dat = dat %>% dplyr::rename(slope=b, slope_se=se) %>% dplyr::filter(se != Inf, se > 0) %>% drop_na()
  df <- data.frame(rsid = dat[, ..SNP_col], pval = dat[,..pval_col])
  colnames(df) = c("rsid", "pval")
  out <- tryCatch({
    ieugwasr::ld_clump(df, clump_kb=clump_kb, clump_r2=clump_r2, clump_p=clump_p, bfile=bfile, plink_bin = plink_bin, pop = pop)
  }, silent = TRUE, error = function(x) return(NA)
  )
  if(length(out)==1) {
    return(NA)
  }
  MRdat <- dat[which(unlist(dat[,..SNP_col]) %in% out$rsid),]
  return(MRdat)
}

dat_path <- "/gpfs/data/gao-lab/people/Sihao/data/GTEx_Analysis_v10_QTLs-GTEx_Analysis_v10_eQTL_all_associations-Brain_Amygdala.v10.allpairs.chr1.txt.gz"

library(ieugwasr)

df <- data.frame(
  rsid = sub("_b38$", "", gsub("_", ":", dat$variant_id)),
  pval = dat$pval_nominal
)

res <- ieugwasr::ld_clump(
  df,
  pop = "EUR",          # 使用 1000G EUR 远程参考
  clump_kb = 250,
  clump_r2 = 0.1,
  clump_p = 0.01
)