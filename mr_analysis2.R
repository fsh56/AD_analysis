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

b_exp <- data$b_eqtl
se_exp <- data$se_eqtl
b_out <- data$b_gwas
se_out <- data$se_gwas
snp_names <- data$rsid

# run MR
cat("running MR analysis...\n")
mr_object <- MendelianRandomization::mr_input(
  bx = b_exp,
  bxse = se_exp,
  by = b_out,
  byse = se_out,
  snps = snp_names
)

cat("running IVW method...\n")
ivw_result <- mr_ivw(mr_object)

cat("running MR-Egger...\n")
egger_result <- try(mr_egger(mr_object), silent = TRUE)

cat("running MR-Median...\n")
median_result <- try(mr_median(mr_object, weighting = "weighted"), silent = TRUE)

cat("running cML method...\n")
cml_result <- try(mr_cML(mr_object, MA = TRUE, DP = FALSE, num_pert = 100), silent = TRUE)

# compute heter
het_Q <- ivw_result@Heter.Stat[1]
het_df <- ivw_result@Heter.Stat[2]
het_pval <- pchisq(het_Q, het_df, lower.tail = FALSE)
het_I2 <- max(0, (het_Q - het_df) / het_Q)

# get summary
summary_results <- data.frame(
  gene_id = sub(".*_(ENSG[0-9]+\\.[0-9]+)_.*", "\\1", basename(input_file)),
  n_snps = nrow(data),
  
  # IVW
  ivw_beta = ivw_result@Estimate,
  ivw_se = ivw_result@StdError,
  ivw_pvalue = ivw_result@Pvalue,
  
  # MR-Egger
  egger_beta = if(!inherits(egger_result, 'try-error')) egger_result@Estimate else NA,
  egger_se = if(!inherits(egger_result, 'try-error')) egger_result@StdError.Est else NA,
  egger_pvalue = if(!inherits(egger_result, 'try-error')) egger_result@Pvalue.Est else NA,
  egger_intercept = if(!inherits(egger_result, 'try-error')) egger_result@Intercept else NA,
  egger_intercept_se = if(!inherits(egger_result, 'try-error')) egger_result@StdError.Int else NA,
  egger_intercept_pval = if(!inherits(egger_result, 'try-error')) egger_result@Pvalue.Int else NA,
  
  # Median
  median_beta = if(!inherits(median_result, 'try-error')) median_result@Estimate else NA,
  median_se = if(!inherits(median_result, 'try-error')) median_result@StdError else NA,
  median_pvalue = if(!inherits(median_result, 'try-error')) median_result@Pvalue else NA,
  
  # cML
  cml_beta = if(!inherits(cml_result, 'try-error')) cml_result@Estimate else NA,
  cml_se = if(!inherits(cml_result, 'try-error')) cml_result@StdError else NA,
  cml_pvalue = if(!inherits(cml_result, 'try-error')) cml_result@Pvalue else NA,
  
  # heter stat
  heterogeneity_Q = het_Q,
  heterogeneity_Q_df = het_df,
  heterogeneity_Q_pval = het_pval,
  heterogeneity_I2 = het_I2
)

# save results
cat("saving results...\n")
write.csv(summary_results, paste0(output_prefix, "_mr_summary.csv"), row.names = FALSE)
cat("  mr result is saved to:", paste0(output_prefix, "_mr_summary.csv\n"))
cat("MR analysis done!\n")