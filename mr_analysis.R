#!/usr/bin/env Rscript
# Methods: IVW, MR-Median, MR-Egger, cML

# Load packages
suppressPackageStartupMessages({
  library(TwoSampleMR)
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

b_exp <- data$b_eqtl
se_exp <- data$se_eqtl
b_out <- data$b_gwas
se_out <- data$se_gwas
snp_names <- data$rsid

results <- list()
cat("Running MR analysis...\n")

# prepare data
exposure_dat <- data.frame(
  SNP = snp_names,
  beta.exposure = b_exp,
  se.exposure = se_exp,
  effect_allele.exposure = data$alt,
  other_allele.exposure = data$ref,
  pval.exposure = data$p_eqtl
)

outcome_dat <- data.frame(
  SNP = snp_names,
  beta.outcome = b_out,
  se.outcome = se_out,
  effect_allele.outcome = data$alt,
  other_allele.outcome = data$ref,
  pval.outcome = data$p_gwas
)

harmonized <- suppressMessages(harmonise_data(exposure_dat, outcome_dat))
cat("Harmonized variants:", nrow(harmonized), "\n\n")
# Heterogeneity tests
cat("running heterogeneity test...\n")
heterogeneity <- mr_heterogeneity(harmonized)
# Pleiotropy test
cat("running pleiotropy test...\n")
pleiotropy <- mr_pleiotropy_test(harmonized)

# create MR object for MendelianRandomization package
mr_object <- MendelianRandomization::mr_input(
  bx = b_exp,
  bxse = se_exp,
  by = b_out,
  byse = se_out,
  snps = snp_names
)

# IVW
cat("Running IVW method...\n")
ivw_mr <- MendelianRandomization::mr_ivw(mr_object)
results$ivw_beta <- ivw_mr@Estimate
results$ivw_se <- ivw_mr@StdError
results$ivw_pvalue <- ivw_mr@Pvalue

# MR-Egger
cat("Running MR-Egger...\n")
egger <- try(mr_egger(mr_object), silent = TRUE)
if(!inherits(egger, 'try-error') && !is.null(egger)) {
  results$egger_beta <- egger@Estimate
  results$egger_se <- egger@StdError.Est
  results$egger_pvalue <- egger@Pvalue.Est
  results$egger_intercept <- egger@Intercept
} else {
  results$egger_beta <- NA
  results$egger_se <- NA
  results$egger_pvalue <- NA
  results$egger_intercept <- NA
  cat("  Warning: MR-Egger failed\n")
}

# Weighted Median
cat("Running MR-Median...\n")
median_mr <- try(mr_median(mr_object), silent = TRUE)
if(!inherits(median_mr, 'try-error') && !is.null(median_mr)) {
  results$median_beta <- median_mr@Estimate
  results$median_se <- median_mr@StdError.Est
  results$median_pvalue <- median_mr@Pvalue.Est
} else {
  results$median_beta <- NA
  results$median_se <- NA
  results$median_pvalue <- NA
  cat("  Warning: MR-Median failed\n")
}

# cML
cat("Running cML method...\n")
cml <- try(MendelianRandomization::mr_cML(
  mr_object, 
  MA = TRUE, 
  DP = FALSE, 
  num_pert = 100
), silent = TRUE)

if(!inherits(cml, 'try-error') && !is.null(cml)) {
  results$cml_beta <- cml@Estimate
  results$cml_se <- cml@StdError
  results$cml_pvalue <- cml@Pvalue
} else {
  results$cml_beta <- NA
  results$cml_se <- NA
  results$cml_pvalue <- NA
  cat("  Warning: cML failed\n")
}

results$n_snps <- nrow(data)
results$n_harmonized <- nrow(harmonized)
mr_results <- as.data.frame(results)

# save results
cat("saving results...\n")
write.csv(mr_results, 
          paste0(output_prefix, "_mr_results.csv"), 
          row.names = FALSE)
cat("  MR results:", paste0(output_prefix, "_mr_results.csv\n"))
write.csv(heterogeneity, 
          paste0(output_prefix, "_heterogeneity.csv"), 
          row.names = FALSE)
cat("  Heterogeneity:", paste0(output_prefix, "_heterogeneity.csv\n"))
write.csv(pleiotropy, 
          paste0(output_prefix, "_pleiotropy.csv"), 
          row.names = FALSE)
cat("  Pleiotropy:", paste0(output_prefix, "_pleiotropy.csv\n"))
cat("MR analysis done!\n")