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
  
  # MR-Egger
  cat("running MR-Egger...\n")
  egger_result <- try(mr_egger(mr_object), silent = TRUE)
  if (!inherits(egger_result, 'try-error')) {
    egger_summary <- data.frame(
      gene_id = gene,
      method = "MR-Egger",
      n_snps = n_snps,
      beta = egger_result@Estimate,
      se = egger_result@StdError.Est,
      pvalue = egger_result@Pvalue.Est,
      intercept = egger_result@Intercept,
      intercept_se = egger_result@StdError.Int,
      intercept_pvalue = egger_result@Pvalue.Int
    )
    output_file <- file.path(output_dir, paste0(gene, "_Egger.csv"))
    write.csv(egger_summary, output_file, row.names = FALSE)
    cat("saved to:", output_file, "\n")
  } else {
    cat("MR-Egger failed\n")
  }
  
  rm(gene_data, b_exp, se_exp, b_out, se_out, snp_names, mr_object)
  if (exists("egger_result")) rm(egger_result, egger_summary)
  gc(verbose = FALSE)
  cat("\n")
}
rm(data, unique_genes)
gc(verbose = FALSE)

cat("All genes processed!\n")
