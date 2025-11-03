#!/usr/bin/env Rscript
library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)

setwd("/Users/fengsihao/AD_analysis/mr_results1027_amy")
input_file  <- "amy_gene_summary_w50kb_p0.001_all.txt"
output_file <- "nsnps_distribution_by_chr.png"
# Set file paths
#input_file <- "/gpfs/data/gao-lab/people/Sihao/amyloid_ba9/gene_info/BA9_gene_summary_w50kb_p0.001_all.txt"
#output_file <- "/gpfs/data/gao-lab/people/Sihao/amyloid_ba9/gene_info/gene_summary_table.csv"

dat <- fread(input_file)
stopifnot(all(c("chr","n_snps") %in% names(dat)))

gene_counts <- dat %>%
  group_by(chr) %>% summarise(n_gene = n(), .groups = "drop") %>%
  right_join(tibble::tibble(chr = 1:22), by = "chr") %>%
  mutate(n_gene = tidyr::replace_na(n_gene, 0)) %>%
  arrange(chr)

dat <- dat %>%
  mutate(chr = as.integer(chr)) %>%
  left_join(gene_counts, by = "chr") %>%
  mutate(chr_label = paste0("chr", chr, "\n(n=", n_gene, ")"))

levels_chr <- paste0("chr", gene_counts$chr, "\n(n=", gene_counts$n_gene, ")")
dat$chr_label <- factor(dat$chr_label, levels = levels_chr)

chr_means <- dat %>%
  group_by(chr_label) %>%
  summarise(mean_n_snps = mean(n_snps, na.rm = TRUE), .groups = "drop")

colors_22 <- c("#E6B3CC", "#B3D9E6", "#CCE6B3", "#E6CCB3",
               "#D9B3E6", "#B3E6CC", "#E6D9B3", "#CCB3E6",
               "#B3CCE6", "#E6E6B3", "#E6B3D9", "#B3E6D9",
               "#D9E6B3", "#E6D9CC", "#CCB3D9", "#D9CCE6",
               "#B3D9CC", "#CCE6D9", "#E6CCE6", "#D9E6CC",
               "#CCE6E6", "#E6CCD9")

p <- ggplot(dat, aes(x = chr_label, y = n_snps, fill = chr_label, color = chr_label)) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.6, show.legend = FALSE) +
  geom_boxplot(width = 0.22, alpha = 0.85, outlier.shape = NA, color = "black", show.legend = FALSE) +
  geom_point(data = chr_means, aes(y = mean_n_snps), inherit.aes = FALSE,
             x = chr_means$chr_label, color = "red", size = 2, shape = 16) +
  labs(x = "", y = "n_snps",
       title = "Distribution of SNPs per Gene across Chromosomes") +
  theme_classic() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y  = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "none"
  ) +
  scale_fill_manual(values = colors_22, drop = FALSE) +
  scale_color_manual(values = colors_22, drop = FALSE)

ggsave(output_file, plot = p, width = 14, height = 6, dpi = 300)
cat("Saved:", output_file, "\n")
