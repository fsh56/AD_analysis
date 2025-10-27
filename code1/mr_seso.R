#!/usr/bin/env Rscript
# MR analysis by gene from clumped data
# Methods: MR-SESO

# ===== 设置临时目录 (必须在最开始) =====
Sys.setenv(TMPDIR = "/scratch/sfeng56/tmp")
Sys.setenv(TEMP = "/scratch/sfeng56/tmp")
Sys.setenv(TMP = "/scratch/sfeng56/tmp")
dir.create("/scratch/sfeng56/tmp", showWarnings = FALSE, recursive = TRUE)
cat("Using temporary directory:", Sys.getenv("TMPDIR"), "\n\n")
# =========================================

# Load packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(mvtnorm)
  library(LaplacesDemon)
  library(invgamma)
})

# Source required files
cat("Compiling C++ code...\n")
Rcpp::sourceCpp("/scratch/sfeng56/draft/seso.cpp")
cat("C++ code compiled successfully\n")

source("/scratch/sfeng56/draft/helper_functions.R")

# Define args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript mr_seso.R <input_file.txt> <output_dir>")
}

input_file <- args[1]
output_dir <- args[2]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Extract p-value threshold from filename
p_threshold <- gsub(".*_p([0-9.]+)\\.txt$", "\\1", basename(input_file))
cat("Extracted p-value threshold:", p_threshold, "\n")

cat("Reading data from:", input_file, "\n")
data <- fread(input_file)

unique_genes <- unique(data$gene_id)
cat("Number of unique genes:", length(unique_genes), "\n\n")

data_clean <- data %>% 
  filter(se_eqtl != Inf, se_gwas != Inf, se_eqtl > 0, se_gwas > 0) %>% 
  drop_na()

all_seso_results <- list()

# Process each gene
for (gene in unique_genes) {
  cat("Processing gene:", gene, "\n")
  
  gene_data <- data[gene_id == gene]
  n_snps <- nrow(gene_data)
  cat("  Number of SNPs:", n_snps, "\n")
  
  if (n_snps < 2) {
    cat("  Skipping - insufficient variants (minimum 2 required)\n\n")
    next
  }
  
  if (n_snps >= 1000) {
    cat("  Skipping - too many variants (maximum 999)\n\n")
    next
  }
  
  # Extract data
  b_exp <- gene_data$b_eqtl
  se_exp <- gene_data$se_eqtl
  b_out <- gene_data$b_gwas
  se_out <- gene_data$se_gwas
  
  if (any(is.infinite(c(se_exp, se_out))) || any(c(se_exp, se_out) <= 0)) {
    cat("  Skipping - invalid standard errors\n\n")
    next
  }
  
  # Run SESO
  tryCatch({
    niter <- 20000
    beta_est <- gibbs_seso_uhp_only(niter, b_out, b_exp, se_out, se_exp)
    burnin <- floor(niter * 0.5)
    beta_est <- beta_est[(burnin + 1):niter]
    res_summary <- get_summary(beta_est)
    
    seso_result <- data.frame(
      beta = res_summary$beta_est,
      se = res_summary$beta_se,
      pvalue = res_summary$beta_pval
    )
    
    all_seso_results[[length(all_seso_results) + 1]] <- data.frame(
      gene_id = gene,
      p_threshold = p_threshold,
      method = "seso",
      n_snps = n_snps,
      beta = seso_result$beta,
      se = seso_result$se,
      pvalue = seso_result$pvalue
    )
    
    cat("  SESO completed successfully\n\n")
  }, error = function(e) {
    cat("  Error in SESO analysis:", conditionMessage(e), "\n\n")
  })
  
  # Clean up
  rm(gene_data, b_exp, se_exp, b_out, se_out)
  if (exists("seso_result")) rm(seso_result)
  if (exists("beta_est")) rm(beta_est)
  if (exists("res_summary")) rm(res_summary)
  
  # 每处理10个基因清理一次内存
  if (which(unique_genes == gene) %% 10 == 0) {
    gc(verbose = FALSE)
  }
}

# Combine and save results
cat("\nSaving results...\n")

if (length(all_seso_results) > 0) {
  seso_combined <- do.call(rbind, all_seso_results)
  output_file <- file.path(output_dir, paste0("amyloid_seso_p", p_threshold, ".csv"))
  write.csv(seso_combined, output_file, row.names = FALSE)
  cat("MR-SESO results saved:", output_file, "\n")
  cat("  Total genes with SESO results:", nrow(seso_combined), "\n")
} else {
  cat("Warning: No SESO results generated\n")
}

rm(data, unique_genes, data_clean)
gc(verbose = FALSE)

cat("\n====================================\n")
cat("All genes processed (SESO)!\n")
cat("====================================\n")