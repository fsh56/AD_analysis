#!/usr/bin/env Rscript
# LD per gene

library(data.table)
library(tidyverse)
library(ieugwasr)
library(genetics.binaRies)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript ld_per_gene.R <input_file> <output_dir> <clump_kb> <clump_p>")
}

input_file <- args[1]
output_dir <- args[2]
clump_kb_param <- as.numeric(args[3])
clump_p_param <- as.numeric(args[4])

# Fixed r2 threshold
clump_r2_param <- 0.1

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Print parameters
cat("\n=== LD Clumping Parameters (Per Gene) ===\n")
cat("Input file:", input_file, "\n")
cat("Output directory:", output_dir, "\n")
cat("clump_kb:", clump_kb_param, "\n")
cat("clump_r2:", clump_r2_param, "\n")
cat("clump_p:", clump_p_param, "\n")
cat("==========================================\n\n")

# LD clumping function (per gene)
clump_per_gene <- function(dat, 
                           SNP_col = "rsid", 
                           pval_col = "p_eqtl", 
                           clump_kb = 250, 
                           clump_r2 = 0.1, 
                           clump_p = 0.01, 
                           bfile = "/scratch/sfeng56/data/1kg/EUR", 
                           plink_bin = genetics.binaRies::get_plink_binary(), 
                           pop = "EUR") {
  
  # Prepare data frame for clumping
  df <- data.frame(rsid = dat[[SNP_col]], pval = dat[[pval_col]])
  df <- df[complete.cases(df), ]
  
  # Check if there's data to clump
  if (nrow(df) == 0) {
    return(data.frame())  # Return empty data frame instead of NA
  }
  
  # Filter by p-value threshold first
  df <- df[df$pval <= clump_p, ]
  
  if (nrow(df) == 0) {
    return(data.frame())  # No significant SNPs
  }
  
  # LD clumping
  out <- tryCatch({
    ieugwasr::ld_clump(df, 
                       clump_kb = clump_kb, 
                       clump_r2 = clump_r2, 
                       clump_p = clump_p, 
                       bfile = bfile, 
                       plink_bin = plink_bin, 
                       pop = pop)
  }, error = function(e) {
    return(NULL)
  })
  
  # Check if clumping was successful
  if (is.null(out) || nrow(out) == 0) {
    return(data.frame())
  }
  
  # Return clumped data
  MRdat <- dat[dat[[SNP_col]] %in% out$rsid, ]
  return(MRdat)
}

# Main processing
cat("Processing file:", input_file, "\n")

# Read harmonized data
dat <- fread(input_file)
cat("Total SNPs before clumping:", nrow(dat), "\n")
cat("Unique genes:", length(unique(dat$gene_id)), "\n\n")

# Get list of unique genes
gene_list <- unique(dat$gene_id)
cat("Starting per-gene clumping for", length(gene_list), "genes...\n\n")

# Initialize results list
all_clumped <- list()

# Process each gene separately
for (i in seq_along(gene_list)) {
  gene <- gene_list[i]
  
  # Progress indicator
  if (i %% 100 == 0) {
    cat(sprintf("Processed %d / %d genes (%.1f%%)\n", 
                i, length(gene_list), 100 * i / length(gene_list)))
  }
  
  # Extract data for this gene
  gene_dat <- dat[dat$gene_id == gene, ]
  
  # Skip if no data
  if (nrow(gene_dat) == 0) {
    next
  }
  
  # Perform clumping for this gene
  gene_clumped <- clump_per_gene(gene_dat, 
                                 SNP_col = "rsid", 
                                 pval_col = "p_eqtl", 
                                 clump_kb = clump_kb_param, 
                                 clump_r2 = clump_r2_param, 
                                 clump_p = clump_p_param, 
                                 bfile = "/scratch/sfeng56/data/1kg/EUR", 
                                 plink_bin = genetics.binaRies::get_plink_binary(),
                                 pop = "EUR")
  
  # Store results if not empty
  if (nrow(gene_clumped) > 0) {
    all_clumped[[gene]] <- gene_clumped
  }
}

cat("\nClumping completed for all genes!\n\n")

# Combine all clumped results
if (length(all_clumped) > 0) {
  clumped_df <- do.call(rbind, all_clumped)
  
  # Add summary statistics per gene
  gene_summary <- clumped_df %>%
    group_by(gene_id) %>%
    summarise(
      n_snps = n(),
      min_pval = min(p_eqtl, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # Print summary
  cat("=== Clumping Summary ===\n")
  cat("Total SNPs after clumping:", nrow(clumped_df), "\n")
  cat("Genes with clumped SNPs:", length(unique(clumped_df$gene_id)), "\n")
  cat("Genes with no SNPs after clumping:", 
      length(gene_list) - length(unique(clumped_df$gene_id)), "\n")
  cat("\nSNPs per gene statistics:\n")
  print(summary(gene_summary$n_snps))
  
  # Save main results
  output_file <- file.path(output_dir, 
                           sprintf("BA9_clumped_w%dkb_p%s.txt", 
                                   clump_kb_param, 
                                   format(clump_p_param, scientific = FALSE)))
  
  fwrite(clumped_df, output_file, sep = ",", quote = FALSE)
  cat("\nMain results saved to:", output_file, "\n")
  
  # Save gene summary
  summary_file <- file.path(output_dir, 
                            sprintf("BA9_gene_summary_w%dkb_p%s.txt", 
                                    clump_kb_param, 
                                    format(clump_p_param, scientific = FALSE)))
  
  fwrite(gene_summary, summary_file, sep = "\t", quote = FALSE)
  cat("Gene summary saved to:", summary_file, "\n")
  
  # Save detailed log
  log_file <- file.path(output_dir, 
                        sprintf("BA9_clumping_log_w%dkb_p%s.txt", 
                                clump_kb_param, 
                                format(clump_p_param, scientific = FALSE)))
  
  log_content <- c(
    "=== LD Clumping Log (Per Gene) ===",
    paste("Date:", Sys.time()),
    paste("Input file:", input_file),
    paste("clump_kb:", clump_kb_param),
    paste("clump_r2:", clump_r2_param),
    paste("clump_p:", clump_p_param),
    "",
    "=== Summary Statistics ===",
    paste("Total genes in input:", length(gene_list)),
    paste("Genes with significant eQTLs (p <", clump_p_param, "):", 
          length(unique(clumped_df$gene_id))),
    paste("Total SNPs before clumping:", nrow(dat)),
    paste("Total SNPs after clumping:", nrow(clumped_df)),
    paste("Reduction rate:", 
          sprintf("%.1f%%", 100 * (1 - nrow(clumped_df) / nrow(dat)))),
    "",
    "=== SNPs per Gene Distribution ===",
    capture.output(print(summary(gene_summary$n_snps))),
    "",
    "=== Genes with Most SNPs (Top 10) ===",
    capture.output(print(head(gene_summary[order(-gene_summary$n_snps), ], 10)))
  )
  
  writeLines(log_content, log_file)
  cat("Detailed log saved to:", log_file, "\n")
  
} else {
  cat("No SNPs passed clumping for any gene\n")
  
  # Create empty output file to indicate completion
  output_file <- file.path(output_dir, 
                           sprintf("clumped_w%dkb_p%s.txt", 
                                   clump_kb_param, 
                                   format(clump_p_param, scientific = FALSE)))
  writeLines("# No SNPs passed clumping", output_file)
}

cat("\n=== Process Completed ===\n")
cat("Done!\n")