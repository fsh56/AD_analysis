library(data.table)
library(ggplot2)
library(dplyr)

setwd("/Users/fengsihao/AD_analysis/amyloid_Brain_Frontal_Cortex_BA9_1.txt.gz/")
output_dir <- "/Users/fengsihao/AD_analysis/plots/egger/"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

csv_files <- list.files(pattern = "*_Egger.csv$", recursive = TRUE, full.names = TRUE)
all_results <- rbindlist(lapply(csv_files, function(file) {
  tryCatch({
    fread(file)
  }, error = function(e) {
    message(paste("Error reading:", file))
    return(NULL)
  })
}))

print(head(all_results))

# compute inflation factor
chisq <- qchisq(1 - all_results$pvalue, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
print(paste("Inflation factor =", round(lambda, 3)))

# dist of p-vals
p1 <- ggplot(all_results, aes(x = pvalue)) +
  geom_histogram(bins = 50, fill = "#3182bd", color = "white", alpha = 0.8) +
  geom_hline(yintercept = nrow(all_results)/50, linetype = "dashed", 
             color = "red", linewidth = 0.8) +
  labs(title = "Distribution of p-values using MR_Egger",
       x = "p-value",
       y = "Frequency") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

png(paste0(output_dir, "pvalue_dist.png"), width = 8, height = 6, units = "in", res = 300)
print(p1)
dev.off()

# qqplot
all_results <- all_results %>%
  arrange(pvalue) %>%
  mutate(
    observed = -log10(pvalue),
    expected = -log10(ppoints(n()))
  )

N <- nrow(all_results)
all_results <- all_results %>%
  mutate(
    clower = -log10(qbeta(0.025, 1:N, N:1)),
    cupper = -log10(qbeta(0.975, 1:N, N:1))
  )

p2 <- ggplot(all_results, aes(x = expected, y = observed)) +
  geom_ribbon(aes(ymin = clower, ymax = cupper), 
              fill = "grey70", alpha = 0.5) +
  geom_abline(intercept = 0, slope = 1, 
              color = "black", linewidth = 0.8) +
  geom_point(color = "#3182bd", alpha = 0.6, size = 2, shape = 1) +
  labs(
    title = "Q-Q Plot using MR_Egger",
    x = "expected -log_10(p)",
    y = "observed -log_10(p)"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
    axis.line = element_blank()
  ) +
  annotate("text", x = 0.2, y = 11.5,
           label = paste("Lambda =", round(lambda, 3)),
           hjust = 0, size = 4.5)

png(paste0(output_dir, "qq_plot.png"), width = 7, height = 7, units = "in", res = 300)
print(p2)
dev.off()

cat("plots saved to:", output_dir, "...\n")