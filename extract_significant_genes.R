#!/usr/bin/env Rscript
# Extract Bonferroni significant genes (simple version)
# Only output gene lists, no CSV files

library(data.table)
library(dplyr)

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript extract_significant_genes_simple.R <mr_results_dir> <output_dir>")
}

mr_dir <- args[1]
output_dir <- args[2]

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set significance threshold
p_threshold <- 0.05

cat("========================================\n")
cat("Extracting Bonferroni Significant Genes\n")
cat("========================================\n")
cat("P-value threshold:", p_threshold, "\n")
cat("MR results directory:", mr_dir, "\n")
cat("Output directory:", output_dir, "\n\n")

# Find all MR result files (exclude merged files)
mr_files <- list.files(mr_dir, pattern = "^(IVW|Egger|WeightedMedian)_.*\\.csv$", 
                       full.names = TRUE)
mr_files <- mr_files[!grepl("all_gwas", mr_files)]

cat(sprintf("Found %d MR result files\n\n", length(mr_files)))

# Parse filename function
parse_mr_filename <- function(filename) {
  base <- gsub("\\.csv$", "", basename(filename))
  parts <- strsplit(base, "_")[[1]]
  
  method <- parts[1]
  param_idx <- which(grepl("^w\\d+kb$", parts))
  
  if (length(param_idx) == 0) {
    stop("Cannot find parameter pattern w*kb in filename: ", filename)
  }
  
  gwas_tissue_parts <- parts[2:(param_idx - 1)]
  tissue_start <- which(grepl("^Brain", gwas_tissue_parts, ignore.case = TRUE))
  
  if (length(tissue_start) > 0) {
    gwas_name <- paste(gwas_tissue_parts[1:(tissue_start[1] - 1)], collapse = "_")
    tissue_name <- paste(gwas_tissue_parts[tissue_start[1]:length(gwas_tissue_parts)], 
                         collapse = "_")
  } else {
    gwas_name <- paste(gwas_tissue_parts[1:(length(gwas_tissue_parts) - 1)], 
                       collapse = "_")
    tissue_name <- gwas_tissue_parts[length(gwas_tissue_parts)]
  }
  
  window_size <- gsub("w|kb", "", parts[param_idx])
  clump_pval <- gsub("p", "", parts[param_idx + 1])
  
  return(list(
    method = method,
    gwas_name = gwas_name,
    tissue_name = tissue_name,
    window_size = window_size,
    clump_pval = clump_pval
  ))
}

# Process each file
total_files <- 0
total_genes <- 0

for (file in mr_files) {
  cat("Processing:", basename(file), "\n")
  
  tryCatch({
    # Parse filename
    parsed <- parse_mr_filename(file)
    method <- parsed$method
    gwas_name <- parsed$gwas_name
    tissue_name <- parsed$tissue_name
    window_size <- parsed$window_size
    clump_pval <- parsed$clump_pval
    
    # Read data
    data <- fread(file)
    
    # Calculate Bonferroni threshold
    bonf_threshold <- p_threshold / nrow(data)
    
    # Filter Bonferroni significant genes
    bonf_genes <- data %>% 
      filter(pvalue < bonf_threshold) %>%
      arrange(pvalue) %>%
      pull(gene_id)
    
    if (length(bonf_genes) > 0) {
      # Output filename
      output_file <- file.path(output_dir, 
                               sprintf("significant_genes_%s_%s_%s_w%s_p%s_list.txt",
                                       method, gwas_name, tissue_name, 
                                       window_size, clump_pval))
      
      # Save gene list
      write.table(bonf_genes, output_file, 
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
      
      cat("  -> ", length(bonf_genes), "Bonferroni significant genes\n")
      cat("  -> Saved to:", basename(output_file), "\n")
      
      total_files <- total_files + 1
      total_genes <- total_genes + length(bonf_genes)
    } else {
      cat("  -> No Bonferroni significant genes\n")
    }
    
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
  })
  
  cat("\n")
}

cat("========================================\n")
cat("Extraction Complete!\n")
cat("========================================\n")
cat(sprintf("Generated %d gene lists with %d total genes\n", total_files, total_genes))
cat("\nOutput directory:", output_dir, "\n")