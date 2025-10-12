library(data.table)
library(tidyverse)

# Get filename
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("no filename provided")
}
filename <- args[1]

# define directories
input_dir <- "/gpfs/data/gao-lab/people/Bowei/rosmap_all_unzip_gao/"
output_dir <- "/gpfs/data/gao-lab/people/Bowei/rosmap_all_unzip_gao/with_rsid/"
ref_table_path <- "/gpfs/data/gao-lab/people/Sihao/data/GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz"

# create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# read reference table
cat("reading reference table...\n")
ref_table <- fread(ref_table_path, select = c("variant_id", "rs_id_dbSNP155_GRCh38p13"))
cat("reference table loaded:", nrow(ref_table), "variants\n")

# read GWAS data
cat("reading GWAS file:", filename, "\n")
gwas <- fread(paste0(input_dir, filename), , header = FALSE)
if (grepl("linear", filename)) {
  # Quantitative trait association
  colnames(gwas) <- c("CHR", "SNP", "BP", "A1", "TEST", "NMISS", "BETA", "STAT", "P")
  cat("Detected linear regression format\n")
} else if (grepl("logistic", filename)) {
  # Case-control association
  colnames(gwas) <- c("CHR", "SNP", "BP", "A1", "TEST", "NMISS", "OR", "STAT", "P")
  cat("Detected logistic regression format\n")
}

# convert format and join
cat("converting and mapping...\n")
gwas_mapped <- gwas %>%
  rename(variant_id_original = SNP) %>%
  mutate(variant_id_ref = paste0("chr", gsub(":", "_", variant_id_original), "_b38")) %>%
  left_join(ref_table, by = c("variant_id_ref" = "variant_id")) %>%
  rename(rsid = rs_id_dbSNP155_GRCh38p13)

# filter NAs
gwas_mapped_filtered <- gwas_mapped %>% 
  filter(!is.na(rsid))

# save
output_filename <- gsub("\\.gz$", ".with_rsid.txt.gz", filename)
output_path <- paste0(output_dir, output_filename)

cat("saving to:", output_path, "\n")
fwrite(gwas_mapped_filtered, output_path)
cat("All done")