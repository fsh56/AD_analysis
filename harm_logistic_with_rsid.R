#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)

# get args
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript harmonize_logistic_with_rsid.R <gwas_file_name> <tissue_name>\n")
  cat("Example: Rscript harmonize_logistic_with_rsid.R dcfdx_ad.assoc.logistic.with_rsid.txt.gz Brain_Amygdala\n")
  quit(status = 1)
}

gwas_file = args[1]
tissue_name = args[2]

# Define paths
gwas_dir = "/gpfs/data/gao-lab/people/Sihao/data/gwas_with_rsid"
eqtl_dir = "/gpfs/data/gao-lab/people/Sihao/data/eqtl_with_rsid"
output_dir = "/gpfs/data/gao-lab/people/Sihao/ad_analysis/harmonized_data"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# process GWAS data
cat("processing GWAS data:", gwas_file, "...\n")
gwas_path = file.path(gwas_dir, gwas_file)
if (!file.exists(gwas_path)) {
  stop("GWAS file not found: ", gwas_path)
}

# read GWAS data
gwas_data_tmp = fread(gwas_path)
gwas_data = gwas_data_tmp[
  , c("chr_tmp", "pos_tmp", "REF", "EFF") := tstrsplit(variant_id_original, ":", fixed = TRUE)][
    , b_gwas := log(OR)][
    , se_gwas := abs(b_gwas/STAT)][
      !is.infinite(se_gwas) & se_gwas != 0 & !is.na(se_gwas)][
        , .(rsid, REF, EFF, b_gwas, se_gwas, p_gwas=P)] %>% 
  distinct(rsid, .keep_all = TRUE)

# process eQTL data for each chr
for (chr in 1:22) {
  cat("processing chr", chr, "...\n")
  eqtl_file = paste0(tissue_name, "_chr", chr, "_with_rsid_SNP_only.txt.gz")
  eqtl_path = file.path(eqtl_dir, eqtl_file)
  
  if (!file.exists(eqtl_path)) {
    cat("Warning: eQTL file not found:", eqtl_path, "\n")
    cat("Skipping chr", chr, "\n\n")
    next
  }
  
  # read eQTL data
  cat("reading eQTL data...\n")
  eqtl_data = fread(eqtl_path)
  # process eQTL
  eqtl_bk = eqtl_data[
    , c("chr_eqtl", "pos_eqtl", "ref", "alt", "build") := tstrsplit(variant_id, "_", fixed = TRUE)][
      , .(gene_id, variant_id, rsid, ref, alt, tss_distance, af, 
          ma_samples, ma_count, p_eqtl = pval_nominal, b_eqtl = slope, se_eqtl = slope_se)]
  # remove invalid one
  eqtl_bk = eqtl_bk[
    !is.infinite(se_eqtl) & se_eqtl != 0 & !is.na(se_eqtl) & 
      !is.infinite(b_eqtl) & !is.na(b_eqtl) & 
      p_eqtl != 1 & !is.na(p_eqtl)]
  
  # harmonize
  cat("harmonizing...\n")
  tmp = merge(eqtl_bk, gwas_data, by = "rsid", all.x = FALSE, all.y = FALSE)
  tmp = tmp[
    , align := fcase(
      alt == EFF & ref == REF, 1,
      alt == REF & ref == EFF, -1,
      default = NA_real_
    )][
      !(is.na(align) & !is.na(b_gwas))][  # Remove SNPs that can't be aligned
        , b_gwas := b_gwas * align][
          !is.na(b_gwas) & !is.na(se_gwas)]  # Keep only SNPs with valid GWAS data
  
  # select cols
  tmp = tmp[, .(gene_id, variant_id, rsid, ref, alt, tss_distance, af,
                ma_samples, ma_count, p_eqtl, b_eqtl, se_eqtl, b_gwas, se_gwas, p_gwas)]
  
  cat("harmonized SNPs:", nrow(tmp), "\n")
  
  # Save results
  gwas_base = gsub("\\.with_rsid\\.txt\\.gz$", "", gwas_file)
  gwas_base = gsub("\\.assoc\\.logistic.*", "", gwas_base)
  
  output_file = file.path(output_dir, 
                          paste0("eqtl_", gwas_base, "_", tissue_name, "_", chr, ".txt.gz"))
  
  fwrite(tmp, output_file)
  cat("saved to:", output_file, "\n")
  # clean
  rm(eqtl_data, eqtl_bk, tmp)
  gc()
  cat("chr", chr, "done!\n\n")
}

cat("All done!\n")