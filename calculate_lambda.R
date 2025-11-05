#!/usr/bin/env Rscript
# =============================================================================
# Calculate Genomic Inflation Factor (Lambda) for MR Results
# =============================================================================
# This script calculates lambda (genomic inflation factor) for each MR method
# and generates a summary table
# =============================================================================

library(data.table)
library(dplyr)

# =============================================================================
# Function to calculate genomic inflation factor (lambda)
# =============================================================================
calculate_lambda <- function(pvalues) {
  # Remove NA values
  pvalues <- pvalues[!is.na(pvalues)]
  
  # Check if there are valid p-values
  if (length(pvalues) == 0) {
    return(NA)
  }
  
  # Convert p-values to chi-square statistics (1 df)
  chi_squared <- qchisq(pvalues, df = 1, lower.tail = FALSE)
  
  # Calculate lambda (median method)
  lambda <- median(chi_squared, na.rm = TRUE) / qchisq(0.5, df = 1)
  
  return(lambda)
}

# =============================================================================
# Function to parse filename and extract metadata
# =============================================================================
parse_filename <- function(filename) {
  # Remove .csv extension
  base_name <- gsub("\\.csv$", "", basename(filename))
  
  # Split by underscore
  parts <- strsplit(base_name, "_")[[1]]
  
  # Extract components
  # Format: {method}_{gwas}_{tissue1}_{tissue2}_w{kb}kb_p{pval}_r{r2}
  method <- parts[1]
  gwas_name <- parts[2]
  
  # Find where tissue name ends (look for pattern starting with 'w')
  tissue_parts <- c()
  param_start <- NULL
  for (i in 3:length(parts)) {
    if (grepl("^w\\d+kb$", parts[i])) {
      param_start <- i
      break
    }
    tissue_parts <- c(tissue_parts, parts[i])
  }
  
  tissue_name <- paste(tissue_parts, collapse = "_")
  
  # Extract parameters
  window_size <- gsub("w", "", gsub("kb", "", parts[param_start]))
  pval <- gsub("p", "", parts[param_start + 1])
  r2 <- gsub("r", "", parts[param_start + 2])
  
  return(list(
    method = method,
    gwas_name = gwas_name,
    tissue_name = tissue_name,
    window_size = window_size,
    pval = pval,
    r2 = r2
  ))
}

# =============================================================================
# Function to process a single MR result file
# =============================================================================
process_mr_file <- function(file_path) {
  cat(sprintf("Processing: %s\n", basename(file_path)))
  
  # Parse filename to extract metadata
  metadata <- parse_filename(file_path)
  
  # Read the file
  tryCatch({
    data <- fread(file_path)
    
    # Check if pvalue column exists
    if (!"pvalue" %in% colnames(data)) {
      cat("  Warning: 'pvalue' column not found\n")
      return(NULL)
    }
    
    # remove 
    n_snps_to_remove <- 3
    top_n_to_remove <- 0 
    data_filtered <- data %>% 
      filter(n_snps > n_snps_to_remove) %>%
      arrange(pvalue) %>%
      slice(-(1:top_n_to_remove))
    
    # Calculate lambda
    lambda <- calculate_lambda(data_filtered$pvalue)
    
    # Get number of tests
    n_genes <- nrow(data_filtered)
    n_valid_pvals <- sum(!is.na(data_filtered$pvalue))
    
    cat(sprintf("  Lambda: %.4f (n_genes: %d, valid_pvals: %d)\n", 
                lambda, n_genes, n_valid_pvals))
    
    # Create result data frame
    result <- data.frame(
      method = metadata$method,
      gwas_name = metadata$gwas_name,
      tissue_name = metadata$tissue_name,
      window_size = metadata$window_size,
      pval_threshold = metadata$pval,
      r2_threshold = metadata$r2,
      lambda = lambda,
      n_genes = n_genes,
      n_valid_pvals = n_valid_pvals,
      stringsAsFactors = FALSE
    )
    
    return(result)
    
  }, error = function(e) {
    cat(sprintf("  Error: %s\n", e$message))
    return(NULL)
  })
}

# =============================================================================
# Main function
# =============================================================================
calculate_all_lambda <- function(mr_dir, output_file = NULL) {
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("GENOMIC INFLATION FACTOR (LAMBDA) CALCULATION\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("MR results directory:", mr_dir, "\n")
  cat(rep("=", 80), "\n\n", sep = "")
  
  # Find all MR result files
  mr_files <- list.files(
    path = mr_dir,
    pattern = "^(IVW|Egger|WeightedMedian)_.*\\.csv$",
    full.names = TRUE
  )
  
  # Exclude merged files (containing "all_gwas")
  mr_files <- mr_files[!grepl("all_gwas", mr_files)]
  
  cat(sprintf("Found %d MR result files\n\n", length(mr_files)))
  
  if (length(mr_files) == 0) {
    stop("No MR result files found in the directory")
  }
  
  # Process each file
  all_results <- list()
  for (file in mr_files) {
    result <- process_mr_file(file)
    if (!is.null(result)) {
      all_results[[length(all_results) + 1]] <- result
    }
    cat("\n")
  }
  
  # Combine all results
  if (length(all_results) == 0) {
    stop("No valid results obtained")
  }
  
  lambda_table <- do.call(rbind, all_results)
  
  # Sort by method, gwas, tissue
  lambda_table <- lambda_table %>%
    arrange(method, gwas_name, tissue_name)
  
  # Save results
  if (!is.null(output_file)) {
    fwrite(lambda_table, output_file)
    cat(sprintf("Lambda table saved to: %s\n", output_file))
  } else {
    # Auto-generate output filename based on directory
    output_file <- file.path(mr_dir, "lambda_summary.csv")
    fwrite(lambda_table, output_file)
    cat(sprintf("Lambda table saved to: %s\n", output_file))
  }
  return(lambda_table)
}

# =============================================================================
# Command line interface
# =============================================================================
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 1) {
    cat("Usage: Rscript calculate_lambda.R <mr_results_dir> [output_file]\n")
    cat("\n")
    cat("Arguments:\n")
    cat("  mr_results_dir : Directory containing MR result CSV files\n")
    cat("  output_file    : (Optional) Output file path for lambda summary\n")
    cat("                   If not specified, saves to <mr_results_dir>/lambda_summary.csv\n")
    cat("\n")
    cat("Example:\n")
    cat("  Rscript calculate_lambda.R /path/to/mr_results\n")
    cat("  Rscript calculate_lambda.R /path/to/mr_results /path/to/lambda_results.csv\n")
    quit(status = 1)
  }
  
  mr_dir <- args[1]
  output_file <- if (length(args) >= 2) args[2] else NULL
  
  # Check if directory exists
  if (!dir.exists(mr_dir)) {
    stop(sprintf("Directory not found: %s", mr_dir))
  }
  
  # Run analysis
  lambda_table <- calculate_all_lambda(mr_dir, output_file)
}

# =============================================================================
# Example usage in R session
# =============================================================================
# # Load the script
# source("calculate_lambda.R")
# 
# # Calculate lambda for all MR results
# mr_dir <- "/gpfs/data/gao-lab/people/Sihao/AD/amygdala/mr_results"
# lambda_table <- calculate_all_lambda(mr_dir)
# 
# # Or specify custom output file
# lambda_table <- calculate_all_lambda(mr_dir, "my_lambda_results.csv")