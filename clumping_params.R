#!/usr/bin/env Rscript
# LD clumping for harmonized data with configurable parameters
# Usage: Rscript ld_clumping.R <input_file> <output_file> [clump_kb] [clump_r2] [clump_p]

library(data.table)
library(tidyverse)
library(ieugwasr)
library(genetics.binaRies)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript ld_clumping.R <input_file> <output_file> [clump_kb] [clump_r2] [clump_p]")
}

input_file <- args[1]
output_file <- args[2]

# Set default parameters or use command line arguments
clump_kb_param <- ifelse(length(args) >= 3, as.numeric(args[3]), 250)
clump_r2_param <- ifelse(length(args) >= 4, as.numeric(args[4]), 0.1)
clump_p_param <- ifelse(length(args) >= 5, as.numeric(args[5]), 0.01)

# Print parameters
cat("=== LD Clumping Parameters ===\n")
cat("Input file:", input_file, "\n")
cat("Output file:", output_file, "\n")
cat("clump_kb:", clump_kb_param, "\n")
cat("clump_r2:", clump_r2_param, "\n")
cat("clump_p:", clump_p_param, "\n")
cat("==============================\n\n")

# LD clumping function
clump <- function(dat, 
                  SNP_col = "rsid", 
                  pval_col = "p_eqtl", 
                  clump_kb = 250, 
                  clump_r2 = 0.1, 
                  clump_p = 0.01, 
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
  
  cat("Total SNPs for clumping:", nrow(df), "\n")
  
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
  
  # check if clumping was successful
  if (length(out) == 1 && is.na(out)) {
    return(NA)
  }
  
  # return clumped data
  MRdat <- dat[dat[[SNP_col]] %in% out$rsid, ]
  return(MRdat)
}

# main
cat("Processing file:", input_file, "\n")

# read harmonized data
dat <- fread(input_file)

cat("Total SNPs before clumping:", nrow(dat), "\n\n")

# Perform LD clumping with parameters from command line
clumped_df <- clump(dat, 
                    SNP_col = "rsid", 
                    pval_col = "p_eqtl", 
                    clump_kb = clump_kb_param, 
                    clump_r2 = clump_r2_param, 
                    clump_p = clump_p_param, 
                    bfile = "/gpfs/data/gao-lab/people/Sihao/data/1kg/EUR", 
                    plink_bin = genetics.binaRies::get_plink_binary(),
                    pop = "EUR")

# save results
if (length(clumped_df) == 1 && is.na(clumped_df)) {
  cat("LD clumping failed or returned no results\n")
  # Create empty output file to indicate completion
  writeLines("# LD clumping failed", output_file)
} else {
  cat("Number of SNPs after clumping:", nrow(clumped_df), "\n")
  fwrite(clumped_df, output_file, sep = "\t", quote = FALSE)
  cat("Results saved to:", output_file, "\n")
}

cat("Done!\n")