library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide tissue name as argument")
}

tissue <- args[1]
cat("processing tissue:", tissue, "...\n")

# Read reference table
cat("reading lookup table...\n")
lookup_table <- fread("/gpfs/data/gao-lab/people/Sihao/data/GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz")

lookup_clean <- lookup_table %>% 
  select(variant_id, rs_id_dbSNP155_GRCh38p13) %>%
  rename(rsid = rs_id_dbSNP155_GRCh38p13) %>%
  filter(rsid != ".") %>%
  # Filter for SNPs only: both ref and alt must be single nucleotides
  separate(variant_id, into = c("chr", "pos", "ref", "alt", "build"), 
           sep = "_", remove = FALSE) %>%
  filter(
    nchar(ref) == 1 & nchar(alt) == 1,  # Single nucleotide only
    ref %in% c("A", "C", "G", "T"),      # Valid nucleotides
    alt %in% c("A", "C", "G", "T")       # Valid nucleotides
  ) %>%
  select(variant_id, rsid)
cat("Total SNPs in lookup table:", nrow(lookup_clean), "\n")

# Define parameters
chromosomes <- 1:22
gtex_prefix <- "/gpfs/data/gao-lab/people/Sihao/data/gtex_old_files/"
output_dir <- "/gpfs/data/gao-lab/people/Sihao/data/processed_gtex/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

for (chr in chromosomes) {
  
  cat("chr", chr, "...")
  
  input_file <- paste0(gtex_prefix, 
                       "GTEx_Analysis_v10_QTLs-GTEx_Analysis_v10_eQTL_all_associations-",
                       tissue, ".v10.allpairs.chr", chr, ".txt.gz")
  
  if (!file.exists(input_file)) {
    cat("file not found\n")
    next
  }
  
  # read eQTL data
  eqtl_data <- fread(input_file)
  # join with lookup table (only SNPs will be retained)
  eqtl_with_rsid <- eqtl_data %>%
    inner_join(lookup_clean, by = "variant_id")
  
  cat("SNPs retained:", nrow(eqtl_with_rsid), "/", nrow(eqtl_data), 
      sprintf("(%.1f%%)", 100 * nrow(eqtl_with_rsid) / nrow(eqtl_data)), "...")
  
  # Save
  output_file <- paste0(output_dir, tissue, "_chr", chr, "_with_rsid_SNP_only.txt.gz")
  fwrite(eqtl_with_rsid, output_file)
  cat("✓\n")
  
  rm(lookup_table)
  rm(eqtl_data, eqtl_with_rsid)
  gc(verbose = FALSE)
}

cat("\nDone with", tissue, "!\n")
cat("All SNP-only files saved to:", output_dir, "\n")