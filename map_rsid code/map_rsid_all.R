library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide tissue name as argument")
}
tissue <- args[1]
cat("processing tissue:", tissue, "...\n")

# read reference table
cat("reading lookup table...\n")
lookup_table <- fread("/gpfs/data/gao-lab/people/Sihao/data/GTEx_Analysis_2021-02-11_v10_WholeGenomeSeq_953Indiv.lookup_table.txt.gz")
lookup_clean <- lookup_table %>% 
  select(variant_id, rs_id_dbSNP155_GRCh38p13) %>%
  rename(rsid = rs_id_dbSNP155_GRCh38p13) %>%
  filter(rsid != ".")

# define params
chromosomes <- 1:22
gtex_prefix <- "/gpfs/data/gao-lab/people/Sihao/data/gtex_old_files"
output_dir <- "/gpfs/data/gao-lab/people/Sihao/data/processed_gtex/"
for (chr in chromosomes) {
  
  cat("chr", chr, "...")
  
  input_file <- paste0(gtex_prefix, 
                       "GTEx_Analysis_v10_QTLs-GTEx_Analysis_v10_eQTL_all_associations-",
                       tissue, ".v10.allpairs.chr", chr, ".txt.gz")
  if (!file.exists(input_file)) {
    cat("file not found\n")
    next
  }
  
  # read eqtl data
  eqtl_data <- fread(input_file)
  eqtl_with_rsid <- eqtl_data %>%
    inner_join(lookup_clean, by = "variant_id")
  
  # save
  output_file <- paste0(output_dir, tissue, "_chr", chr, "_with_rsid_all_snp.txt.gz")
  fwrite(eqtl_with_rsid, output_file)
  cat("✓\n")
}
cat("Done with", tissue, "!\n")