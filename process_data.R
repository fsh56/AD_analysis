#!/usr/bin/env Rscript
# Harmonize eQTL and GWAS Data 
library(data.table)
library(tidyverse)

# test
brain_tissues = c("Brain_Amygdala")
tissue_short = c("Amygdala")

# Paths
gwas_path = "/scratch/sfeng56/data/gwas_with_rsid/amyloid.assoc.linear.with_rsid.txt.gz"
eqtl_dir = "/scratch/sfeng56/data/eqtl_with_rsid"
output_dir = "/scratch/sfeng56/data/harmo_bk"

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message(paste0("Created output directory: ", output_dir))
}

# ============================================================================
# Step 1: Read and Process GWAS Data
# ============================================================================
message("Reading GWAS data...")
gwas_data = fread(gwas_path)

gwas_with_alleles = gwas_data %>%
  separate(variant_id_ref, 
           into = c("chr_temp", "pos_temp", "REF", "ALT", "build"), 
           sep = "_", 
           remove = FALSE,
           extra = "drop",
           fill = "right") %>%
  dplyr::select(rsid, A1, BETA, STAT, P, REF, ALT, variant_id_ref) %>%
  distinct(rsid, .keep_all = TRUE)

gwas_processed = gwas_with_alleles %>%
  dplyr::rename(EFF = A1) %>%
  # Calculate standard error: SE = BETA / STAT
  mutate(
    se_gwas = ifelse(STAT != 0, abs(BETA / STAT), NA_real_),
    b_gwas = BETA,
    p_gwas = P
  ) %>%
  # Quality filters
  dplyr::filter(
    !is.na(rsid),
    !is.na(b_gwas),
    !is.na(se_gwas),
    !is.infinite(se_gwas),
    se_gwas > 0,
    se_gwas != Inf
  ) %>%
  dplyr::select(rsid, REF, EFF, ALT, b_gwas, se_gwas, p_gwas)
message(paste0("Processed GWAS data: ", nrow(gwas_processed), " variants"))

# ============================================================================
# Step 2: Loop Through Chromosomes and Brain Tissues
# ============================================================================
message("\nStarting harmonization loop...")
message(paste0("Processing ", length(brain_tissues), " brain tissues across 22 chromosomes"))

# Track progress
total_files = length(brain_tissues) * 2
processed_files = 0
failed_files = 0

for (chr in 1:2) {
  message(paste0("\n========== Chr ", chr, " =========="))
  
  for (j in 1:length(brain_tissues)) {
    tissue = brain_tissues[j]
    tissue_name = tissue_short[j]
    message(paste0("Processing: ", j, "/", length(brain_tissues), " - ", tissue))
    
    # Construct eQTL file path
    eqtl_file = file.path(
      eqtl_dir,
      paste0(tissue, "_chr", chr, "_with_rsid_SNP_only.txt.gz")
    )
    
    # Check if file exists
    if (!file.exists(eqtl_file)) {
      message(paste0("  WARNING: File not found - ", eqtl_file))
      failed_files = failed_files + 1
      next
    }
    
    tryCatch({
      # Read eQTL data
      eqtl_data = fread(eqtl_file)
      
      # Process eQTL: extract gene name from gene_id (format: ENSG00000227232.5)
      eqtl_processed = eqtl_data %>%
        mutate(gene_name = str_replace(gene_id, "\\..*$", "")) %>%
        dplyr::rename(
          b_eqtl = slope,
          se_eqtl = slope_se,
          p_eqtl = pval_nominal
        ) %>%
        dplyr::filter(
          !is.na(se_eqtl),
          !is.infinite(se_eqtl),
          se_eqtl != 0,
          se_eqtl != Inf,
          !is.na(p_eqtl),
          p_eqtl != 1,
          !is.na(b_eqtl)
        )
      
      # Extract allele info from variant_id (format: chr1_13550_G_A_b38)
      eqtl_processed = eqtl_processed %>%
        separate(variant_id, 
                 into = c("chr_id", "pos", "ref", "alt", "build"), 
                 sep = "_", 
                 remove = FALSE, 
                 extra = "drop",
                 fill = "right") %>%
        dplyr::rename(variant_id_eqtl = variant_id)
      
      message(paste0("  eQTL variants: ", nrow(eqtl_processed)))
      
      # ============================================================================
      # Step 3: Harmonization - Join eQTL with GWAS
      # ============================================================================
      
      harmonized = eqtl_processed %>%
        left_join(gwas_processed, by = "rsid") %>%
        dplyr::filter(!is.na(b_gwas)) %>%
        # Harmonize alleles: compare eQTL (ref/alt) with GWAS (REF/EFF)
        mutate(
          align = case_when(
            # Direct match: eQTL alt = GWAS effect allele
            alt == EFF & ref == REF ~ 1,
            # Flipped: eQTL ref = GWAS effect allele (need to flip beta)
            ref == EFF & alt == REF ~ -1,
            # No match: strand issue or multi-allelic
            TRUE ~ NA_real_
          )
        ) %>%
        # Filter out mismatches
        dplyr::filter(!is.na(align)) %>%
        # Apply alignment to flip GWAS beta if needed
        mutate(b_gwas = b_gwas * align) %>%
        # Final quality check
        dplyr::filter(
          !is.na(b_gwas),
          !is.na(se_gwas),
          !is.infinite(b_gwas),
          !is.infinite(se_gwas)
        ) %>%
        dplyr::select(
          gene_name,
          gene_id,
          rsid,
          variant_id_eqtl,
          tss_distance,
          af,
          ma_samples,
          ma_count,
          b_eqtl,
          se_eqtl,
          p_eqtl,
          b_gwas,
          se_gwas,
          p_gwas
        )
      
      message(paste0("  Harmonized variants: ", nrow(harmonized)))
      
      # Write output file (format: eqtl_ad_Cortex_1.txt.gz)
      output_file = file.path(
        output_dir,
        paste0("eqtl_ad_", tissue_name, "_", chr, ".txt.gz")
      )
      
      fwrite(harmonized, output_file)
      message(paste0("  Saved: ", output_file))
      
      processed_files = processed_files + 1
      
    }, error = function(e) {
      message(paste0("  ERROR processing ", tissue, " chr", chr, ": ", e$message))
      failed_files = failed_files + 1
    })
    
  } # End tissue loop
  
} # End chromosome loop

# ============================================================================
# Summary
# ============================================================================

message("\n========== Harmonization Complete ==========")
message(paste0("Total expected files: ", total_files))
message(paste0("Successfully processed: ", processed_files))
message(paste0("Failed/Missing: ", failed_files))
message(paste0("Output directory: ", output_dir))
message("\nDone!")