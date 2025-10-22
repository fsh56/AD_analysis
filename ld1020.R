#!/usr/bin/env Rscript
# LD clumping for harmonized eQTL-GWAS data with configurable parameters

library(data.table)
library(tidyverse)
library(ieugwasr)
library(genetics.binaRies)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript ld_clumping_batch.R <input_file> <output_dir> <clump_kb> <clump_p>")
}

input_file <- args[1]
output_dir <- args[2]
clump_kb_param <- as.numeric(args[3])
clump_p_param <- as.numeric(args[4])

# Fixed r2 threshold
clump_r2_param <- 0.1

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Print parameters
cat("\n=== LD Clumping Parameters ===\n")
cat("Input file:", input_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("clump_kb:", clump_kb_param, "\n")
cat("clump_r2:", clump_r2_param, "\n")
cat("clump_p:", clump_p_param, "\n")
cat("==============================\n\n")

# LD clumping function (from original script)
clump <- function(dat, 
                  SNP_col = "rsid", 
                  pval_col = "p_eqtl", 
                  clump_kb = 250, 
                  clump_r2 = 0.1, 
                  clump_p = 0.01, 
                  bfile = "/scratch/sfeng56/data/1kg/EUR", 
                  plink_bin = genetics.binaRies::get_plink_binary(), 
                  pop = "EUR") {
  
  df <- data.frame(rsid = dat[[SNP_col]], pval = dat[[pval_col]])
  df <- df[complete.cases(df), ]
  
  # Check if there's data to clump
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
  
  # Check if clumping was successful
  if (length(out) == 1 && is.na(out)) {
    return(NA)
  }
  
  # Return clumped data
  MRdat <- dat[dat[[SNP_col]] %in% out$rsid, ]
  return(MRdat)
}

# Main processing
cat("Processing file:", input_file, "\n")

# Read harmonized data
dat <- fread(input_file)
cat("Total SNPs before clumping:", nrow(dat), "\n")
cat("Unique genes:", length(unique(dat$gene_id)), "\n\n")
output_file <- file.path(output_dir, 
                         sprintf("clumped_w%dkb_p%s.txt", 
                                 clump_kb_param, 
                                 format(clump_p_param, scientific = FALSE)))
cat("Output file:", output_file, "\n\n")

# Perform LD clumping with parameters from command line
clumped_df <- clump(dat, 
                    SNP_col = "rsid", 
                    pval_col = "p_eqtl", 
                    clump_kb = clump_kb_param, 
                    clump_r2 = clump_r2_param, 
                    clump_p = clump_p_param, 
                    bfile = "/scratch/sfeng56/data/1kg/EUR", 
                    plink_bin = genetics.binaRies::get_plink_binary(),
                    pop = "EUR")

# Save results
if (length(clumped_df) == 1 && is.na(clumped_df)) {
  cat("LD clumping failed or returned no results\n")
  # Create output file to indicate completion
  writeLines("# LD clumping failed", output_file)
} else {
  cat("Number of SNPs after clumping:", nrow(clumped_df), "\n")
  cat("Unique genes after clumping:", length(unique(clumped_df$gene_id)), "\n")
  fwrite(clumped_df, output_file, sep = ",", quote = FALSE)
  cat("Results saved to:", output_file, "\n")
}

cat("\nDone!\n")