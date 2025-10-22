#!/usr/bin/env Rscript
# Methods: IVW, MR-Egger

# Load packages
suppressPackageStartupMessages({
  library(MendelianRandomization)
  library(data.table)
  library(dplyr)
})

# Define args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript mr.R <input_file.txt.gz> <output_dir>")
}

input_file <- args[1]
output_dir <- args[2]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Reading data from:", input_file, "\n")
data <- fread(input_file)

# Extract p-value threshold from filename
p_threshold <- sub(".*_p(0\\.[0-9]+)\\.txt\\.gz", "\\1", basename(input_file))
cat("P-value threshold:", p_threshold, "\n")

unique_genes <- unique(data$gene_id)
cat("Number of unique genes:", length(unique_genes), "\n")

# Initialize result lists
all_ivw_results <- list()
all_egger_results <- list()

# Process each gene
for (i in seq_along(unique_genes)) {
  gene <- unique_genes[i]
  
  if (i %% 10 == 0) {
    cat(sprintf("Progress: %d/%d genes processed\n", i, length(unique_genes)))
  }
  
  gene_data <- data[gene_id == gene]
  n_snps <- nrow(gene_data)
  
  if (n_snps < 3) {
    next
  }
  
  # Create MR object
  b_exp <- gene_data$b_eqtl
  se_exp <- gene_data$se_eqtl
  b_out <- gene_data$b_gwas
  se_out <- gene_data$se_gwas
  snp_names <- gene_data$rsid
  
  mr_object <- try(MendelianRandomization::mr_input(
    bx = b_exp,
    bxse = se_exp,
    by = b_out,
    byse = se_out,
    snps = snp_names
  ), silent = TRUE)
  
  if (inherits(mr_object, 'try-error')) {
    next
  }
  
  # IVW
  ivw_result <- try(mr_ivw(mr_object), silent = TRUE)
  if (!inherits(ivw_result, 'try-error')) {
    het_Q <- ivw_result@Heter.Stat[1]
    het_df <- ivw_result@Heter.Stat[2]
    het_pval <- pchisq(het_Q, het_df, lower.tail = FALSE)
    het_I2 <- max(0, (het_Q - het_df) / het_Q)
    
    all_ivw_results[[length(all_ivw_results) + 1]] <- data.frame(
      gene_id = gene,
      p_threshold = p_threshold,
      method = "IVW",
      n_snps = n_snps,
      beta = ivw_result@Estimate,
      se = ivw_result@StdError,
      pvalue = ivw_result@Pvalue,
      heterogeneity_Q = het_Q,
      heterogeneity_Q_df = het_df,
      heterogeneity_Q_pval = het_pval,
      heterogeneity_I2 = het_I2
    )
  }
  
  # MR-Egger
  egger_result <- try(mr_egger(mr_object), silent = TRUE)
  if (!inherits(egger_result, 'try-error')) {
    all_egger_results[[length(all_egger_results) + 1]] <- data.frame(
      gene_id = gene,
      p_threshold = p_threshold,
      method = "MR-Egger",
      n_snps = n_snps,
      beta = egger_result@Estimate,
      se = egger_result@StdError.Est,
      pvalue = egger_result@Pvalue.Est,
      intercept = egger_result@Intercept,
      intercept_se = egger_result@StdError.Int,
      intercept_pvalue = egger_result@Pvalue.Int
    )
  }
  
  # Clean up memory
  rm(gene_data, b_exp, se_exp, b_out, se_out, snp_names, mr_object)
  if (exists("ivw_result")) rm(ivw_result)
  if (exists("egger_result")) rm(egger_result)
  
  if (i %% 50 == 0) {
    gc(verbose = FALSE)
  }
}

# Combine and save results
cat("\nSaving results...\n")

if (length(all_ivw_results) > 0) {
  ivw_combined <- do.call(rbind, all_ivw_results)
  output_file <- file.path(output_dir, paste0("IVW_p", p_threshold, ".csv"))
  write.csv(ivw_combined, output_file, row.names = FALSE)
  cat("IVW results saved:", output_file, "\n")
  cat("  Total genes with IVW results:", nrow(ivw_combined), "\n")
}

if (length(all_egger_results) > 0) {
  egger_combined <- do.call(rbind, all_egger_results)
  output_file <- file.path(output_dir, paste0("Egger_p", p_threshold, ".csv"))
  write.csv(egger_combined, output_file, row.names = FALSE)
  cat("Egger results saved:", output_file, "\n")
  cat("  Total genes with Egger results:", nrow(egger_combined), "\n")
}

# Clean up
rm(data, unique_genes, all_ivw_results, all_egger_results)
gc(verbose = FALSE)

cat("\nProcessing complete!\n")
cat("=" , rep("=", 70), "=\n", sep = "")