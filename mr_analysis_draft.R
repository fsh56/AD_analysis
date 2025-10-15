#!/usr/bin/env Rscript

# Methods: IVW, MR-Median, MR-Egger, cML
# load packages
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(MendelianRandomization)
  library(data.table)
})

# define args
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript mr_analysis.R <input_file.txt.gz> <output_prefix>")
}

input_file <- args[1]
output_prefix <- args[2]

# read data
cat("reading harmonized data after clumping from:", input_file, "\n")
data <- fread(input_file)
cat("total selected variants:", nrow(data), "\n")

# check for minimum number of variants
if (nrow(data_clean) < 3) {
  stop("Insufficient number of variants for MR analysis after filtering")
}

b_exp <- data_clean$b_eqtl
se_exp <- data_clean$se_eqtl
b_out <- data_clean$b_gwas
se_out <- data_clean$se_gwas
snp_names <- data_clean$rsid

cat("running MR analysis...\n")
# format dataset
exposure_dat <- data.frame(
  SNP = snp_names,
  beta.exposure = b_exp,
  se.exposure = se_exp,
  effect_allele.exposure = data_clean$alt,
  other_allele.exposure = data_clean$ref,
  pval.exposure = data_clean$p_eqtl
)
outcome_dat <- data.frame(
  SNP = snp_names,
  beta.outcome = b_out,
  se.outcome = se_out,
  effect_allele.outcome = data_clean$alt,
  other_allele.outcome = data_clean$ref,
  pval.outcome = data_clean$p_gwas
)

harmonized <- suppressMessages(harmonise_data(exposure_dat, outcome_dat))
cat("Harmonized variants:", nrow(harmonized), "\n")
mr_results <- mr(harmonized, method_list = c(
  "ivw",
  "mr_egger",
  "mr_median"
))

# do heterogeneity tests
heterogeneity <- mr_heterogeneity(harmonized)
# do pleiotropy test
pleiotropy <- mr_pleiotropy_test(harmonized)

# create mr object 
mr_object <- mr_input(
  bx = b_exp,
  bxse = se_exp,
  by = b_out,
  byse = se_out,
  snps = snp_names
)

# IVW
cat("\running IVW method...\n")
ivw_mr <- mr_ivw(mr_object)
ivw_summary <- data.frame(
  method = "IVW",
  estimate = ivw_mr@Estimate,
  se = ivw_mr@StdError,
  pval = ivw_mr@Pvalue,
  ci_lower = ivw_mr@CILower,
  ci_upper = ivw_mr@CIUpper
)

# MR-Egger
cat("\running MR-Egger...\n")
egger_mr <- mr_egger(mr_object)
egger_summary <- data.frame(
  method = "MR-Egger",
  estimate = egger_mr@Estimate,
  se = egger_mr@StdError.Est,
  pval = egger_mr@Pvalue.Est,
  ci_lower = egger_mr@CILower.Est,
  ci_upper = egger_mr@CIUpper.Est,
  intercept = egger_mr@Intercept,
  intercept_pval = egger_mr@Pvalue.Int
)

# Weighted Median
cat("\running MR-Median...\n")
median_mr <- mr_median(mr_object)
median_summary <- data.frame(
  method = "MR-Median",
  estimate = median_mr@Estimate,
  se = median_mr@StdError,
  pval = median_mr@Pvalue,
  ci_lower = median_mr@CILower,
  ci_upper = median_mr@CIUpper
)

# cML
cat("\running cML-MA method...\n")
tryCatch({
  cml_result <- mr_conmix(mr_object, CIMin = -1, CIMax = 1, CIStep = 0.001)
  cml_summary <- data.frame(
    method = "cML-MA",
    estimate = cml_result@Estimate,
    se = cml_result@StdError,
    pval = cml_result@Pvalue,
    ci_lower = cml_result@CILower,
    ci_upper = cml_result@CIUpper
  )
}, error = function(e) {
  cat("Warning: cML-MA failed:", e$message, "\n")
  cml_summary <<- data.frame(
    method = "cML-MA",
    estimate = NA,
    se = NA,
    pval = NA,
    ci_lower = NA,
    ci_upper = NA
  )
})

# save results
mr_results <- rbind(
  ivw_summary[, c("method", "estimate", "se", "pval", "ci_lower", "ci_upper")],
  egger_summary[, c("method", "estimate", "se", "pval", "ci_lower", "ci_upper")],
  median_summary[, c("method", "estimate", "se", "pval", "ci_lower", "ci_upper")],
  cml_summary[, c("method", "estimate", "se", "pval", "ci_lower", "ci_upper")]
)
write.csv(mr_results, 
          paste0(output_prefix, "_mr_results.csv"), 
          row.names = FALSE)
write.csv(heterogeneity, 
          paste0(output_prefix, "_heterogeneity.csv"), 
          row.names = FALSE)
write.csv(pleiotropy, 
          paste0(output_prefix, "_pleiotropy.csv"), 
          row.names = FALSE)

cat("MR analysis done!")