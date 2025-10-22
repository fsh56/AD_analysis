#!/usr/bin/env Rscript
# Methods: IVW, MR-Median, MR-Egger, cML

# Load packages
suppressPackageStartupMessages({
  library(MendelianRandomization)
  library(data.table)
})

# Define args
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript mr_analysis.R <input_file.txt.gz> <output_prefix>")
}

input_file <- args[1]
output_prefix <- args[2]

# read data
cat("reading harmonized data after clumping from:\n")
cat("  ", input_file, "\n")
data <- fread(input_file)
cat("total selected variants:", nrow(data), "\n\n")

# Check for minimum number of variants
if (nrow(data) < 3) {
  stop("Insufficient number of variants for MR analysis (minimum 3 required)")
}

# Extract gene ID
gene_id <- sub(".*_(ENSG[0-9]+\\.[0-9]+)_.*", "\\1", basename(input_file))

b_exp <- data$b_eqtl
se_exp <- data$se_eqtl
b_out <- data$b_gwas
se_out <- data$se_gwas
snp_names <- data$rsid

# Create MR object
cat("running MR analysis...\n")
mr_object <- MendelianRandomization::mr_input(
  bx = b_exp,
  bxse = se_exp,
  by = b_out,
  byse = se_out,
  snps = snp_names
)

# IVW
cat("running IVW method...\n")
ivw_result <- try(mr_ivw(mr_object), silent = TRUE)
if (!inherits(ivw_result, 'try-error')) {
  # compute heterogeneity
  het_Q <- ivw_result@Heter.Stat[1]
  het_df <- ivw_result@Heter.Stat[2]
  het_pval <- pchisq(het_Q, het_df, lower.tail = FALSE)
  het_I2 <- max(0, (het_Q - het_df) / het_Q)
  
  ivw_summary <- data.frame(
    gene_id = gene_id,
    method = "IVW",
    n_snps = nrow(data),
    beta = ivw_result@Estimate,
    se = ivw_result@StdError,
    pvalue = ivw_result@Pvalue,
    heterogeneity_Q = het_Q,
    heterogeneity_Q_df = het_df,
    heterogeneity_Q_pval = het_pval,
    heterogeneity_I2 = het_I2
  )
  write.csv(ivw_summary, paste0(output_prefix, "_IVW_summary.csv"), row.names = FALSE)
  cat("  IVW result saved to:", paste0(output_prefix, "_IVW_summary.csv\n"))
} else {
  cat("  IVW analysis failed\n")
}

# MR-Egger
cat("running MR-Egger...\n")
egger_result <- try(mr_egger(mr_object), silent = TRUE)
if (!inherits(egger_result, 'try-error')) {
  egger_summary <- data.frame(
    gene_id = gene_id,
    method = "MR-Egger",
    n_snps = nrow(data),
    beta = egger_result@Estimate,
    se = egger_result@StdError.Est,
    pvalue = egger_result@Pvalue.Est,
    intercept = egger_result@Intercept,
    intercept_se = egger_result@StdError.Int,
    intercept_pvalue = egger_result@Pvalue.Int
  )
  write.csv(egger_summary, paste0(output_prefix, "_Egger_summary.csv"), row.names = FALSE)
  cat("  MR-Egger result saved to:", paste0(output_prefix, "_Egger_summary.csv\n"))
} else {
  cat("  MR-Egger analysis failed\n")
}

cat("\nMR analysis done!\n")