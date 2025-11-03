#!/usr/bin/env Rscript
# =============================================================================
# MR Analysis for Brain_Amygdala - One GWAS per job
# Methods: IVW, MR-Egger, Weighted Median
# Input: Clumped data across all chromosomes (1-22)
# Output: Three method-specific result files (combining all chromosomes)
# =============================================================================

# Load required packages
suppressPackageStartupMessages({
  library(MendelianRandomization)
  library(data.table)
  library(dplyr)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("Usage: Rscript mr_analysis_amygdala.R <gwas_name> <tissue_name> <clump_dir> <output_dir> <clump_kb> <clump_p> <clump_r2>
       Example: Rscript mr_analysis_amygdala.R gpath Brain_Amygdala /path/to/clumped /path/to/output 50 0.001 0.1")
}

gwas_name <- args[1]
tissue_name <- args[2]
clump_dir <- args[3]
output_dir <- args[4]
clump_kb <- as.numeric(args[5])
clump_p <- as.numeric(args[6])
clump_r2 <- as.numeric(args[7])

# Format parameter strings for filenames
pval_str <- format(clump_p, scientific = FALSE)
r2_str <- format(clump_r2, scientific = FALSE)

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# =============================================================================
# Setup and logging
# =============================================================================
cat("\n")
cat(rep("=", 80), "\n", sep = "")
cat("MR ANALYSIS: ", gwas_name, " - ", tissue_name, "\n", sep = "")
cat(rep("=", 80), "\n", sep = "")
cat("GWAS:", gwas_name, "\n")
cat("Tissue:", tissue_name, "\n")
cat("Parameters: w=", clump_kb, "kb, p=", clump_p, ", r2=", clump_r2, "\n", sep = "")
cat("Clumped data directory:", clump_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat(rep("=", 80), "\n\n", sep = "")

# =============================================================================
# Read clumped data from all chromosomes
# =============================================================================
cat("Reading clumped data from all chromosomes...\n")

all_data <- list()
chromosomes_found <- c()

for (chr in 1:22) {
  input_file <- file.path(clump_dir, 
                         sprintf("clumped_%s_%s_chr%d_w%dkb_p%s_r%s.txt.gz", 
                                 gwas_name, tissue_name, chr, 
                                 clump_kb, pval_str, r2_str))
  
  if (file.exists(input_file)) {
    cat(sprintf("  Reading chr%d... ", chr))
    
    # Read data
    chr_data <- tryCatch({
      fread(input_file)
    }, error = function(e) {
      cat("ERROR\n")
      cat("    Error message:", e$message, "\n")
      return(NULL)
    })
    
    if (!is.null(chr_data) && nrow(chr_data) > 0) {
      # Add chromosome column if not present
      if (!"chr" %in% colnames(chr_data)) {
        chr_data$chr <- chr
      }
      all_data[[chr]] <- chr_data
      chromosomes_found <- c(chromosomes_found, chr)
      cat(sprintf("OK (%d SNPs, %d genes)\n", 
                  nrow(chr_data), 
                  length(unique(chr_data$gene_id))))
    } else {
      cat("EMPTY\n")
    }
  } else {
    cat(sprintf("  Chr%d: File not found\n", chr))
  }
}

# Check if we have any data
if (length(all_data) == 0) {
  stop("ERROR: No clumped data found for any chromosome")
}

# Combine all chromosome data
cat("\nCombining data from", length(chromosomes_found), "chromosomes...\n")
data <- rbindlist(all_data, use.names = TRUE, fill = TRUE)
cat("Total SNPs loaded:", nrow(data), "\n")

# Check required columns
required_cols <- c("gene_id", "rsid", "b_eqtl", "se_eqtl", "b_gwas", "se_gwas")
missing_cols <- setdiff(required_cols, colnames(data))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Get unique genes
unique_genes <- unique(data$gene_id)
cat("Number of unique genes:", length(unique_genes), "\n")
cat("Chromosomes included:", paste(chromosomes_found, collapse = ", "), "\n")
cat(rep("-", 80), "\n\n", sep = "")

# Clean up intermediate data
rm(all_data)
gc(verbose = FALSE)

# =============================================================================
# Initialize result storage
# =============================================================================
all_ivw_results <- list()
all_egger_results <- list()
all_wm_results <- list()

# =============================================================================
# Main MR analysis loop
# =============================================================================
cat("Starting MR analysis for", length(unique_genes), "genes...\n\n")

for (i in seq_along(unique_genes)) {
  gene <- unique_genes[i]
  
  # Progress indicator
  if (i %% 100 == 0) {
    cat(sprintf("[%s] Progress: %d/%d genes (%.1f%%)\n", 
                Sys.time(), i, length(unique_genes), 100*i/length(unique_genes)))
  }
  
  # Extract gene-specific data
  gene_data <- data[gene_id == gene]
  n_snps <- nrow(gene_data)
  
  # Get chromosome(s) for this gene
  gene_chrs <- unique(gene_data$chr)
  
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
      chr = paste(gene_chrs, collapse = ","),
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
      chr = paste(gene_chrs, collapse = ","),
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
      chr = paste(gene_chrs, collapse = ","),
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
  if (i %% 200 == 0) {
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

# Output file naming: {method}_{gwas}_{tissue}_w{kb}_p{pval}_r{r2}.csv
# Example: IVW_gpath_Brain_Amygdala_w50kb_p0.001_r0.1.csv

# Save IVW results
if (length(all_ivw_results) > 0) {
  ivw_combined <- do.call(rbind, all_ivw_results)
  output_file <- file.path(output_dir, 
                          sprintf("IVW_%s_%s_w%dkb_p%s_r%s.csv", 
                                  gwas_name, tissue_name,
                                  clump_kb, pval_str, r2_str))
  fwrite(ivw_combined, output_file)
  cat("✓ IVW results saved:", output_file, "\n")
  cat("  Total genes analyzed:", nrow(ivw_combined), "\n")
  cat("  Significant genes (p < 0.05):", sum(ivw_combined$pvalue < 0.05, na.rm = TRUE), "\n")
  cat("  Mean heterogeneity I²:", sprintf("%.2f%%", mean(ivw_combined$heterogeneity_I2 * 100, na.rm = TRUE)), "\n\n")
} else {
  cat("✗ No IVW results to save\n\n")
}

# Save MR-Egger results
if (length(all_egger_results) > 0) {
  egger_combined <- do.call(rbind, all_egger_results)
  output_file <- file.path(output_dir, 
                          sprintf("Egger_%s_%s_w%dkb_p%s_r%s.csv", 
                                  gwas_name, tissue_name,
                                  clump_kb, pval_str, r2_str))
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
  output_file <- file.path(output_dir, 
                          sprintf("WeightedMedian_%s_%s_w%dkb_p%s_r%s.csv", 
                                  gwas_name, tissue_name,
                                  clump_kb, pval_str, r2_str))
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
cat("ANALYSIS COMPLETE: ", gwas_name, " - ", tissue_name, "\n", sep = "")
cat(rep("=", 80), "\n", sep = "")
cat("Total genes in input:", length(unique_genes), "\n")
cat("Genes with IVW results:", length(all_ivw_results), "\n")
cat("Genes with MR-Egger results:", length(all_egger_results), "\n")
cat("Genes with Weighted Median results:", length(all_wm_results), "\n")
cat("Chromosomes processed:", paste(chromosomes_found, collapse = ", "), "\n")
cat(rep("=", 80), "\n\n", sep = "")

# Final cleanup
rm(data, unique_genes, all_ivw_results, all_egger_results, all_wm_results)
gc(verbose = FALSE)
