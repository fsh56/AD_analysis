#!/usr/bin/env Rscript
# LD clumping for a single gene_id
# Usage: Rscript clumping_per_gene <input_file> <gene_id> <output_file> [clump_kb] [clump_r2] [clump_p]

library(data.table)
library(tidyverse)
library(ieugwasr)
library(genetics.binaRies)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript clumping_per_gene.R <input_file> <gene_id> <output_file> [clump_kb] [clump_r2] [clump_p]")
}

input_file <- args[1]
target_gene_id <- args[2]
output_file <- args[3]

# Set default parameters or use command line arguments
clump_kb_param <- ifelse(length(args) >= 4, as.numeric(args[4]), 20)
clump_r2_param <- ifelse(length(args) >= 5, as.numeric(args[5]), 0.1)
clump_p_param <- ifelse(length(args) >= 6, as.numeric(args[6]), 0.001)

# LD clumping function
clump <- function(dat, 
                  SNP_col = "rsid", 
                  pval_col = "p_eqtl", 
                  clump_kb = 250, 
                  clump_r2 = 0.1, 
                  clump_p = 0.001, 
                  bfile = "/gpfs/data/gao-lab/people/Sihao/data/1kg/EUR", 
                  plink_bin = genetics.binaRies::get_plink_binary(), 
                  pop = "EUR") {
  df <- data.frame(rsid = dat[[SNP_col]], pval = dat[[pval_col]])
  df <- df[complete.cases(df), ]
  # check if there's data to clump
  if (nrow(df) == 0) {
    cat("No valid SNPs for clumping\n")
    return(NA)
  }
  # LD clumping
  out <- tryCatch({
    ieugwasr::ld_clump(df, 
                       clump_kb = clump_kb, 
                       clump_r2 = clump_r2, 
                       clump_p = clump_p, 
                       bfile = bfile, 
                       plink_bin = plink_bin, 
                       pop = pop)
  }, error = function(e) {
    cat("Error in ld_clump:", e$message, "\n")
    return(NA)
  })
  if (length(out) == 1 && is.na(out)) {
    return(NA)
  }
  MRdat <- dat[dat[[SNP_col]] %in% out$rsid, ]
  return(MRdat)
}

# Main processing
cat("Reading input file:", input_file, "\n")

# Read harmonized data
dat <- fread(input_file)

cat("Total SNPs in file:", nrow(dat), "\n")

# Filter for target gene_id
cat("Filtering for gene_id:", target_gene_id, "\n")
gene_dat <- dat[dat$gene_id == target_gene_id, ]

if (nrow(gene_dat) == 0) {
  cat("Error: No SNPs found for gene_id:", target_gene_id, "\n")
  cat("Please check if the gene_id exists in the input file.\n")
  quit(status = 1)
}

cat("SNPs for", target_gene_id, ":", nrow(gene_dat), "\n\n")

# Perform LD clumping
cat("Starting LD clumping...\n")
clumped_df <- clump(gene_dat, 
                    SNP_col = "rsid", 
                    pval_col = "p_eqtl", 
                    clump_kb = clump_kb_param, 
                    clump_r2 = clump_r2_param, 
                    clump_p = clump_p_param, 
                    bfile = "/gpfs/data/gao-lab/people/Sihao/data/1kg/EUR", 
                    plink_bin = genetics.binaRies::get_plink_binary(),
                    pop = "EUR")
# Save results
if (length(clumped_df) == 1 && is.na(clumped_df)) {
  cat("LD clumping failed or returned no results\n")
  # Create empty output file to indicate completion
  writeLines("# LD clumping failed", output_file)
} else {
  fwrite(clumped_df, output_file, sep = ",", quote = FALSE)
  cat("Results saved to:", output_file, "\n")
}

cat("Done!\n")