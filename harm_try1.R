library(data.table)

# read gwas data
cat("reading GWAS data...\n")
out = fread("amyloid.assoc.linear.processed2.txt.gz")
setDT(out)
setkey(out, rsid)

# read eqtl
cat("processing eqtl data...\n")
eqtl_bk = fread("eqtl_Brain_Amygdala_1")[
  , c("chr_eqtl", "pos_eqtl", "ref", "alt", "build") := tstrsplit(variant_id, "_", fixed = TRUE)][
    , rsid := paste0(gsub("chr", "", chr_eqtl), ":", pos_eqtl, ":", ref, ":", alt)][
      , .(gene_name, variant_id, rsid, ref, alt, tss_distance, af, ma_samples, ma_count, p, b, se)]
setDT(eqtl_bk)
setkey(eqtl_bk, rsid)

cat("eQTL SNPs:", nrow(eqtl_bk), "\n")
cat("GWAS SNPs:", nrow(out), "\n")

# harmonize
cat("harmonizing...\n")
tmp = merge(eqtl_bk, out, by = "rsid", all.x = FALSE, all.y = FALSE)[
  , align := fcase(
    alt == A1 & ref == REF, 1,
    alt == REF & ref == A1, -1,
    default = NA_real_
  )][
    !(is.na(align) & !is.na(b_out))][
      , b_out := b_out * align][
        !is.na(b_out) & !is.na(se_out)][
          , .(gene_name, variant_id, rsid, ref, alt, tss_distance, af, 
              ma_samples, ma_count, p, b, se, b_out, se_out)]

cat("harmonized SNPs:", nrow(tmp), "\n")

# save and rm
fwrite(tmp, "harmonized_eqtl_gwas_chr1.txt.gz")
rm(eqtl_bk, out)
gc()

cat("Done!\n")