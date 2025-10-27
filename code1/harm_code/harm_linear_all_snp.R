#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)

# Define args
args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript harm_linear_all_snp.R <gwas_file_name> <tissue_name>\n")
  cat("Example: Rscript harm_linear_all_snp.R amyloid.assoc.linear.gz Brain_Substantia_nigra\n")
  quit(status = 1)
}

gwas_file = args[1]
tissue_name = args[2]

# Define paths
gwas_dir = "/gpfs/data/gao-lab/people/Bowei/rosmap_all_unzip_gao"
eqtl_dir = "/gpfs/data/gao-lab/people/Sihao/data/processed_gtex"
output_dir = "/gpfs/data/gao-lab/people/Sihao/data/harmonized_bk"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# process GWAS data
cat("processing GWAS data (linear):", gwas_file, "...\n")
gwas_path = file.path(gwas_dir, gwas_file)
if (!file.exists(gwas_path)) {
  stop("GWAS file not found: ", gwas_path)
}

# Read GWAS with space-separated format
out = fread(gwas_path, 
            header = FALSE, 
            sep = " ",
            col.names = c("chr", "rsid", "pos", "A1", "test", "nmiss", "beta", "stat", "se_out"))

# Process GWAS
out = out[, .(
  rsid,
  A1,
  REF = sapply(strsplit(rsid, ":"), `[`, 3),
  ALT = sapply(strsplit(rsid, ":"), `[`, 4),
  b_out = beta,
  se_out = se_out,
  p = 2 * pnorm(-abs(stat))  # Calculate p-value from test statistic
)][
  !is.infinite(se_out) & se_out != 0 & !is.na(se_out)
] %>% 
  unique(by = "rsid")

cat("Processed GWAS SNPs:", nrow(out), "\n\n")

# Process eQTL by chromosome
for (chr in 1:22) {
  cat("Processing chr", chr, "\n")
  
  # New eQTL file naming pattern
  eqtl_file = paste0(tissue_name, "_chr", chr, "_with_rsid_SNP_only.txt.gz")
  eqtl_path = file.path(eqtl_dir, eqtl_file)
  
  if (!file.exists(eqtl_path)) {
    cat("Warning: eQTL file not found:", eqtl_path, "\n")
    cat("Skipping chr", chr, "\n\n")
    next
  }
  
  cat("Processing eQTL data...\n")
  
  # Read eQTL with comma-separated format
  eqtl_bk = fread(eqtl_path, 
                  header = TRUE,
                  sep = ",")[
                    , c("chr_eqtl", "pos_eqtl", "ref", "alt", "build") := tstrsplit(variant_id, "_", fixed = TRUE)
                  ]
  
  # Rename columns to match your analysis pipeline
  setnames(eqtl_bk, 
           old = c("pval_nominal", "slope", "slope_se"), 
           new = c("p", "b", "se"),
           skip_absent = TRUE)
  
  # Select relevant columns
  eqtl_bk = eqtl_bk[, .(gene_id, variant_id, rsid, ref, alt, tss_distance, 
                        af, ma_samples, ma_count, p, b, se)]
  
  # Filter out infinite/zero standard errors
  eqtl_bk = eqtl_bk[!is.infinite(se) & se != 0 & se > 0 & !is.na(se)]
  
  cat("eQTL SNPs:", nrow(eqtl_bk), "\n")
  
  # Harmonize
  cat("Harmonizing...\n")
  tmp = merge(eqtl_bk, out, by = "rsid", all.x = FALSE, all.y = FALSE)[
    , align := fcase(
      alt == A1 & ref == REF, 1,
      alt == REF & ref == A1, -1,
      default = NA_real_
    )][
      !(is.na(align) & !is.na(b_out))
    ][
      , b_out := b_out * align
    ][
      !is.na(b_out) & !is.na(se_out) & !is.infinite(se_out) & se_out > 0
    ][
      , .(gene_id, variant_id, rsid, ref, alt, tss_distance, af,
          ma_samples, ma_count, p, b, se, b_out, se_out)
    ]
  
  cat("Harmonized SNPs:", nrow(tmp), "\n")
  
  # Save results
  gwas_base = gsub("\\.gz$", "", gwas_file)
  gwas_base = gsub("\\.assoc\\.linear.*$", "", gwas_base)
  
  output_file = file.path(output_dir, 
                          paste0("harmonized_", gwas_base, "_", tissue_name, "_chr", chr, ".txt.gz"))
  
  fwrite(tmp, output_file)
  cat("Saved to:", output_file, "\n")
  
  # Clean up memory
  rm(eqtl_bk, tmp)
  gc()
  cat("Chr", chr, "done!\n\n")
}

cat("All done!\n")