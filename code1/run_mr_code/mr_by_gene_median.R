#!/usr/bin/env Rscript
# MR analysis by gene from clumped data
# Methods: IVW, MR-Median, MR-Egger

# Load packages
suppressPackageStartupMessages({
  library(MendelianRandomization)
  library(data.table)
  library(dplyr)
})

# define args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript mr_analysis_by_gene.R <input_file.txt.gz> <output_dir>")
}

input_file <- args[1]
output_dir <- args[2]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("reading data from:", input_file, "\n")
data <- fread(input_file)
unique_genes <- unique(data$gene_id)
cat("Number of unique genes:", length(unique_genes), "\n\n")

# process each gene
for (gene in unique_genes) {
  cat("processing gene:", gene, "\n")
  gene_data <- data[gene_id == gene]
  n_snps <- nrow(gene_data)
  cat("  Number of SNPs:", n_snps, "\n")
  if (n_snps < 3) {
    cat("  Skipping - insufficient variants (minimum 3 required)\n\n")
    next
  }
  
  # create MR object
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
    cat("  Failed to create MR object\n\n")
    next
  }

  # Weighted Median
  cat("running Weighted Median...\n")
  median_result <- try(mr_median(mr_object, weighting = "weighted"), silent = TRUE)
  if (!inherits(median_result, 'try-error')) {
    median_summary <- data.frame(
      gene_id = gene,
      method = "Weighted_Median",
      n_snps = n_snps,
      beta = median_result@Estimate,
      se = median_result@StdError,
      pvalue = median_result@Pvalue
    )
    output_file <- file.path(output_dir, paste0(gene, "_Median.csv"))
    write.csv(median_summary, output_file, row.names = FALSE)
    cat("saved to:", output_file, "\n")
  } else {
    cat("MR Weighted Median failed\n")
  }
  
  rm(gene_data, b_exp, se_exp, b_out, se_out, snp_names, mr_object)
  if (exists("median_result")) rm(median_result, median_summary)
  gc(verbose = FALSE)
  cat("\n")
}
rm(data, unique_genes)
gc(verbose = FALSE)

cat("All genes processed!\n")
