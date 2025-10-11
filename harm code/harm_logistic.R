#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)

# process logistic gwas
# define args
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript harm_logistic.R <gwas_file_name> <tissue_name>\n")
  cat("Example: Rscript harm_logistic.R diagnosis.assoc.logistic.gz Brain_Amygdala\n")
  quit(status = 1)
}

gwas_file = args[1]
tissue_name = args[2]

# define paths
gwas_dir = "/gpfs/data/gao-lab/people/Bowei/rosmap_all_unzip_gao"
eqtl_dir = "/gpfs/data/gao-lab/people/Sihao/data"
output_dir = "/gpfs/data/gao-lab/people/Sihao/ad_analysis/harmonized_bk"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# process gwas data
cat("processing GWAS data (logistic):", gwas_file, "...\n")
gwas_path = file.path(gwas_dir, gwas_file)

if (!file.exists(gwas_path)) {
  stop("GWAS file not found: ", gwas_path)
}

out = fread(gwas_path, 
            header = FALSE, 
            col.names = c("chr", "rsid", "pos", "A1", "test", "nmiss", "or", "stat", "p"))[
              , beta := log(or)][                          # convert OR to beta
                , se_out := abs(beta/stat)][
                  !is.infinite(se_out) & se_out != 0 & !is.na(se_out)][
                    , c("chr_tmp", "pos_tmp", "REF", "ALT") := tstrsplit(rsid, ":", fixed = TRUE)][
                      , .(rsid, A1, REF, ALT, b_out = beta, se_out, p)] %>% 
  unique(by = "rsid")

cat("processed GWAS SNPs:", nrow(out), "\n\n")

# process eqtl
for (chr in 1:22) {
  cat("processing chr", chr, "\n")
  
  eqtl_file = paste0("GTEx_Analysis_v10_QTLs-GTEx_Analysis_v10_eQTL_all_associations-",
                     tissue_name, ".v10.allpairs.chr", chr, ".txt.gz")
  eqtl_path = file.path(eqtl_dir, eqtl_file)
  
  if (!file.exists(eqtl_path)) {
    cat("Warning: eQTL file not found:", eqtl_path, "\n")
    cat("Skipping chr", chr, "\n\n")
    next
  }
  
  cat("processing eqtl data...\n")
  eqtl_bk = fread(eqtl_path, header = TRUE)[
    , c("chr_eqtl", "pos_eqtl", "ref", "alt", "build") := tstrsplit(variant_id, "_", fixed = TRUE)][
      , rsid := paste0(gsub("chr", "", chr_eqtl), ":", pos_eqtl, ":", ref, ":", alt)][
        , .(gene_id, variant_id, rsid, ref, alt, tss_distance, af, 
            ma_samples, ma_count, pval_nominal, slope, slope_se)]
  setnames(eqtl_bk, c("pval_nominal", "slope", "slope_se"), c("p", "b", "se"))
  cat("eQTL SNPs:", nrow(eqtl_bk), "\n")
  
  # harmonize
  cat("harmonizing...\n")
  tmp = merge(eqtl_bk, out, by = "rsid", all.x = TRUE, all.y = FALSE, suffixes = c("", "_gwas"))[
    , align := fcase(
      alt == A1 & ref == REF, 1,
      alt == REF & ref == A1, -1,
      default = NA_real_
    )][
      !(is.na(align) & !is.na(b_out))][
        , b_out := b_out * align][
          !is.na(b_out) & !is.na(se_out)][
            , .(gene_id, variant_id, rsid, ref, alt, tss_distance, af,
                ma_samples, ma_count, p, b, se, b_out, se_out)]
  
  cat("harmonized SNPs:", nrow(tmp), "\n")
  
  # save results
  gwas_base = gsub("\\.gz$", "", gwas_file)
  gwas_base = gsub("\\..*$", "", gwas_base)
  output_file = file.path(output_dir, 
                          paste0("harmonized_", gwas_base, "_", tissue_name, "_chr", chr, ".txt.gz"))
  
  fwrite(tmp, output_file)
  cat("saved to:", output_file, "\n")
  
  # rm
  rm(eqtl_bk, tmp)
  gc()
  cat("chr", chr, "done!\n\n")
}

cat("All done\n")