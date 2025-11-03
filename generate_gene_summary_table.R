#!/usr/bin/env Rscript
# Generate gene summary table by chromosome

library(data.table)
library(dplyr)

# Set file paths
setwd("/Users/fengsihao/AD_analysis")
input_file  <- "amy_gene_summary_w50kb_p0.001_all.txt"
output_file <- "amy_gene_summary_table.csv"
#input_file <- "/gpfs/data/gao-lab/people/Sihao/amyloid_ba9/gene_info/BA9_gene_summary_w50kb_p0.001_all.txt"
#output_file <- "/gpfs/data/gao-lab/people/Sihao/amyloid_ba9/gene_info/gene_summary_table.csv"

# Read the data
cat("Reading data from:", input_file, "\n")
dat <- fread(input_file)

# Check if required columns exist
if (!"chr" %in% colnames(dat)) {
  stop("ERROR: Column 'chr' not found in the input file")
}
if (!"n_snps" %in% colnames(dat)) {
  stop("ERROR: Column 'n_snps' not found in the input file")
}

# Generate summary table by chromosome
cat("Generating summary statistics by chromosome...\n")

summary_table <- dat %>%
  group_by(chr) %>%
  summarise(
    min_n_snps = min(n_snps, na.rm = TRUE),
    max_n_snps = max(n_snps, na.rm = TRUE),
    median_n_snps = median(n_snps, na.rm = TRUE),
    mean_n_snps = mean(n_snps, na.rm = TRUE),
    total_gene_num = n(),
    pct_lt_3 = sum(n_snps < 3) / n() * 100,
    .groups = 'drop'
  ) %>%
  arrange(chr)

# Save to CSV
cat("\nSaving results to:", output_file, "\n")
fwrite(summary_table, output_file, sep = ",")
cat("\nDone!\n")