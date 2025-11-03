#!/usr/bin/env Rscript
# =============================================================================
# MR Analysis for BA9 Clumped Data
# Methods: IVW, MR-Egger, Weighted Median
# Input: Chromosome-specific clumped data
# Output: Method-specific results in the same directory as input
# =============================================================================

# Load required packages
suppressPackageStartupMessages({
  library(MendelianRandomization)
  library(data.table)
  library(dplyr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript mr_analysis_ba9.R <chromosome_number> <base_dir>
       Example: Rscript mr_analysis_ba9.R 1 /gpfs/data/gao-lab/people/Sihao/amyloid_ba9")
}

chr <- as.integer(args[1])
base_dir <- args[2]

# =============================================================================
# Setup paths
# =============================================================================
input_dir <- file.path(base_dir, paste0("ldByGene", chr))
input_file <- file.path(input_dir, "BA9_clumped_w50kb_p0.001.txt")

# Output will be saved in the same directory as input
output_dir <- input_dir

# Log header
cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("MR ANALYSIS FOR CHROMOSOME", chr, "\n")
cat(rep("=", 80), "\n", sep = "")
cat("Input file:", input_file, "\n")
cat("Output directory:", output_dir, "\n")
cat(rep("=", 80), "\n\n", sep = "")

# =============================================================================
# Read and validate data
# =============================================================================
if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

cat("Reading clumped data...\n")
data <- fread(input_file)
cat("Total rows loaded:", nrow(data), "\n")

# Check required columns
required_cols <- c("gene_id", "rsid", "b_eqtl", "se_eqtl", "b_gwas", "se_gwas")
missing_cols <- setdiff(required_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Get unique genes
unique_genes <- unique(data$gene_id)
cat("Number of unique genes:", length(unique_genes), "\n")
cat(rep("-", 80), "\n\n", sep = "")

# =============================================================================
# Initialize result storage
# =============================================================================
all_ivw_results <- list()
all_egger_results <- list()
all_wm_results <- list()

# =============================================================================
# Main MR analysis loop
# =============================================================================
cat("Starting MR analysis...\n\n")

for (i in seq_along(unique_genes)) {
  gene <- unique_genes[i]
  
  # Progress indicator
  if (i %% 50 == 0) {
    cat(sprintf("[%s] Progress: %d/%d genes (%.1f%%)\n", 
                Sys.time(), i, length(unique_genes), 100*i/length(unique_genes)))
  }
  
  # Extract gene-specific data
  gene_data <- data[gene_id == gene]
  n_snps <- nrow(gene_data)
  
  # Skip if too few SNPs (need at least 3 for MR-Egger)
  if (n_snps < 3) {
    next
  }
  
  # Prepare MR input data
  b_exp <- gene_data$b_eqtl
  se_exp <- gene_data$se_eqtl
  b_out <- gene_data$b_gwas
  se_out <- gene_data$se_gwas
  snp_names <- gene_data$rsid
  
  # Create MR object
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
  
  # ---------------------------------------------------------------------------
  # Method 1: Inverse Variance Weighted (IVW)
  # ---------------------------------------------------------------------------
  ivw_result <- try(mr_ivw(mr_object), silent = TRUE)
  if (!inherits(ivw_result, 'try-error')) {
    # Calculate heterogeneity statistics
    het_Q <- ivw_result@Heter.Stat[1]
    het_df <- ivw_result@Heter.Stat[2]
    het_pval <- pchisq(het_Q, het_df, lower.tail = FALSE)
    het_I2 <- max(0, (het_Q - het_df) / het_Q)
    
    all_ivw_results[[length(all_ivw_results) + 1]] <- data.frame(
      gene_id = gene,
      chr = chr,
      method = "IVW",
      n_snps = n_snps,
      beta = ivw_result@Estimate,
      se = ivw_result@StdError,
      pvalue = ivw_result@Pvalue,
      heterogeneity_Q = het_Q,
      heterogeneity_Q_df = het_df,
      heterogeneity_Q_pval = het_pval,
      heterogeneity_I2 = het_I2,
      stringsAsFactors = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Method 2: MR-Egger
  # ---------------------------------------------------------------------------
  egger_result <- try(mr_egger(mr_object), silent = TRUE)
  if (!inherits(egger_result, 'try-error')) {
    all_egger_results[[length(all_egger_results) + 1]] <- data.frame(
      gene_id = gene,
      chr = chr,
      method = "MR-Egger",
      n_snps = n_snps,
      beta = egger_result@Estimate,
      se = egger_result@StdError.Est,
      pvalue = egger_result@Pvalue.Est,
      intercept = egger_result@Intercept,
      intercept_se = egger_result@StdError.Int,
      intercept_pvalue = egger_result@Pvalue.Int,
      stringsAsFactors = FALSE
    )
  }
  
  # ---------------------------------------------------------------------------
  # Method 3: Weighted Median
  # ---------------------------------------------------------------------------
  wm_result <- try(mr_median(mr_object, weighting = "weighted", iterations = 10000), 
                   silent = TRUE)
  if (!inherits(wm_result, 'try-error')) {
    all_wm_results[[length(all_wm_results) + 1]] <- data.frame(
      gene_id = gene,
      chr = chr,
      method = "Weighted Median",
      n_snps = n_snps,
      beta = wm_result@Estimate,
      se = wm_result@StdError,
      pvalue = wm_result@Pvalue,
      stringsAsFactors = FALSE
    )
  }
  
  # Memory management
  rm(gene_data, b_exp, se_exp, b_out, se_out, snp_names, mr_object)
  if (exists("ivw_result")) rm(ivw_result)
  if (exists("egger_result")) rm(egger_result)
  if (exists("wm_result")) rm(wm_result)
  
  # Periodic garbage collection
  if (i %% 100 == 0) {
    gc(verbose = FALSE)
  }
}

# =============================================================================
# Save results
# =============================================================================
cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("SAVING RESULTS\n")
cat(rep("=", 80), "\n", sep = "")

# Save IVW results
if (length(all_ivw_results) > 0) {
  ivw_combined <- do.call(rbind, all_ivw_results)
  output_file <- file.path(output_dir, paste0("BA9_chr", chr, "_IVW_results.csv"))
  fwrite(ivw_combined, output_file)
  cat("✓ IVW results saved:", output_file, "\n")
  cat("  Total genes analyzed:", nrow(ivw_combined), "\n")
  cat("  Significant genes (p < 0.05):", sum(ivw_combined$pvalue < 0.05, na.rm = TRUE), "\n\n")
} else {
  cat("✗ No IVW results to save\n\n")
}

# Save MR-Egger results
if (length(all_egger_results) > 0) {
  egger_combined <- do.call(rbind, all_egger_results)
  output_file <- file.path(output_dir, paste0("BA9_chr", chr, "_Egger_results.csv"))
  fwrite(egger_combined, output_file)
  cat("✓ MR-Egger results saved:", output_file, "\n")
  cat("  Total genes analyzed:", nrow(egger_combined), "\n")
  cat("  Significant genes (p < 0.05):", sum(egger_combined$pvalue < 0.05, na.rm = TRUE), "\n")
  cat("  Significant intercept (p < 0.05):", 
      sum(egger_combined$intercept_pvalue < 0.05, na.rm = TRUE), "\n\n")
} else {
  cat("✗ No MR-Egger results to save\n\n")
}

# Save Weighted Median results
if (length(all_wm_results) > 0) {
  wm_combined <- do.call(rbind, all_wm_results)
  output_file <- file.path(output_dir, paste0("BA9_chr", chr, "_WeightedMedian_results.csv"))
  fwrite(wm_combined, output_file)
  cat("✓ Weighted Median results saved:", output_file, "\n")
  cat("  Total genes analyzed:", nrow(wm_combined), "\n")
  cat("  Significant genes (p < 0.05):", sum(wm_combined$pvalue < 0.05, na.rm = TRUE), "\n\n")
} else {
  cat("✗ No Weighted Median results to save\n\n")
}

# =============================================================================
# Summary statistics
# =============================================================================
cat(rep("=", 80), "\n", sep = "")
cat("ANALYSIS COMPLETE FOR CHROMOSOME", chr, "\n")
cat(rep("=", 80), "\n", sep = "")
cat("Total genes in input:", length(unique_genes), "\n")
cat("Genes with IVW results:", length(all_ivw_results), "\n")
cat("Genes with MR-Egger results:", length(all_egger_results), "\n")
cat("Genes with Weighted Median results:", length(all_wm_results), "\n")
cat(rep("=", 80), "\n\n", sep = "")

# Final cleanup
rm(data, unique_genes, all_ivw_results, all_egger_results, all_wm_results)
gc(verbose = FALSE)