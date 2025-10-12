#!/usr/bin/env Rscript
library(data.table)
library(tidyverse)

# parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript map_rsid_improved.R <tissue_name>")
}
tissue <- args[1]


# Define paths (make these configurable)
lookup_file <- "/gpfs/data/gao-lab/people/Sihao/data/GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz"
gtex_prefix <- "/gpfs/data/gao-lab/people/Sihao/data/gtex_old_files/"
output_dir <- "/gpfs/data/gao-lab/people/Sihao/data/processed_gtex/"

# Verify output directory exists, create if not
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("created output directory:", output_dir, "\n\n")
}

# Read and prepare lookup table
cat("reading lookup table...\n")
if (!file.exists(lookup_file)) {
  stop("ERROR: lookup table file not found at: ", lookup_file)
}

tryCatch({
  lookup_table <- fread(lookup_file, showProgress = FALSE)
  cat("✓ Lookup table loaded successfully\n")
  cat("  Total variants:", nrow(lookup_table), "\n")
  # clean lookup table
  lookup_clean <- lookup_table %>% 
    select(variant_id, rs_id_dbSNP155_GRCh38p13) %>%
    rename(rsid = rs_id_dbSNP155_GRCh38p13) %>%
    filter(rsid != ".")
  rm(lookup_table)
  gc(verbose = FALSE)
}, error = function(e) {
  stop("ERROR reading lookup table: ", e$message)
})

# Process each chromosome
cat("Step 2: processing chr...\n")
cat(rep("-", 72), "\n", sep = "")

chromosomes <- 1:22
for (chr in chromosomes) {
  cat(sprintf("processing chr %2d: ", chr))
  input_file <- paste0(
    gtex_prefix, 
    "GTEx_Analysis_v10_QTLs-GTEx_Analysis_v10_eQTL_all_associations-",
    tissue, ".v10.allpairs.chr", chr, ".txt.gz"
  )
  
  # process eqtl file
  tryCatch({
    # read eQTL data
    eqtl_data <- fread(input_file)
    # join with rsid lookup
    eqtl_with_rsid <- eqtl_data %>%
      inner_join(lookup_clean, by = "variant_id")
    # save output
    output_file <- paste0(output_dir, tissue, "_chr", chr, "_with_rsid_all_snp.txt.gz")
    fwrite(eqtl_with_rsid, output_file)
    # clean up memory
    rm(eqtl_data, eqtl_with_rsid)
    gc(verbose = FALSE)
  }, error = function(e) {
    cat("ERROR:", e$message, "\n")
  })
}

cat("Done with", tissue, "!\n")
