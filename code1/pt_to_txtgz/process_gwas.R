# process gwas data
library(data.table)
library(tidyverse)

cat("process GWAS data...\n")
out = fread("/gpfs/data/gao-lab/people/Bowei/rosmap_all_unzip_gao/amyloid.assoc.linear.gz", header = FALSE, 
            col.names = c("chr", "rsid", "pos", "A1", "test", "nmiss", "beta", "stat", "p"))[
              , se_out := abs(beta/stat)][
                !is.infinite(se_out) & se_out != 0 & !is.na(se_out)][
                  , .(rsid, A1, b_out = beta, se_out, p)] %>% unique(by = "rsid") 
cat("processed num of SNPs:", nrow(out), "\n")
fwrite(out, "amyloid.assoc.linear.processed.txt.gz")
cat("Done!\n")


# v2--correct
library(data.table)
library(tidyverse)

cat("process GWAS data...\n")
out = fread("/gpfs/data/gao-lab/people/Bowei/rosmap_all_unzip_gao/amyloid.assoc.linear.gz", 
            header = FALSE, 
            col.names = c("chr", "rsid", "pos", "A1", "test", "nmiss", "beta", "stat", "p"))[
              , se_out := abs(beta/stat)][
                !is.infinite(se_out) & se_out != 0 & !is.na(se_out)][
                  , c("chr_tmp", "pos_tmp", "REF", "ALT") := tstrsplit(rsid, ":", fixed = TRUE)][
                    , .(rsid, A1, REF, ALT, b_out = beta, se_out, p)] %>% unique(by = "rsid") 
cat("processed num of SNPs:", nrow(out), "\n")
fwrite(out, "amyloid.assoc.linear.processed2.txt.gz")
cat("Done!\n")