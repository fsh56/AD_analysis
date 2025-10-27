#!/usr/bin/env Rscript

# LD clumping for harmonized data
# Usage: Rscript ld_clumping.R <input_file> <output_file>

library(data.table)
library(tidyverse)
library(ieugwasr)
library(genetics.binaRies)


# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript ld_clumping.R <input_file> <output_file>")
}

input_file <- args[1]
output_file <- args[2]

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
cat("processing file:", input_file, "\n")
# read harmonized data
dat <- fread(input_file)
# Perform LD clumping
clumped_df <- clump(dat, 
                    SNP_col = "rsid", 
                    pval_col = "p_eqtl", 
                    clump_kb = 250, 
                    clump_r2 = 0.1, 
                    clump_p = 0.01, 
                    bfile = "/gpfs/data/gao-lab/people/Sihao/data/1kg/EUR", 
                    plink_bin = genetics.binaRies::get_plink_binary(),
                    pop = "EUR")

# save results
if (length(clumped_df) == 1 && is.na(clumped_df)) {
  cat("LD clumping failed or returned no results\n")
  # Create empty output file to indicate completion
  writeLines("# LD clumping failed", output_file)
} else {
  cat("number of SNPs after clumping:", nrow(clumped_df), "\n")
  fwrite(clumped_df, output_file, sep = "\t", quote = FALSE)
  cat("fesults saved to:", output_file, "\n")
}

cat("Done!\n")