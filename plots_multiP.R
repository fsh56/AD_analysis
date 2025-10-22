library(data.table)
library(ggplot2)
library(dplyr)

setwd("/Users/fengsihao/AD_analysis/multi-pvals/")
output_dir <- "/Users/fengsihao/AD_analysis/plots/multi-pvals/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

files <- c(
  "Egger_p0.00001.csv",
  "Egger_p0.0001.csv", 
  "Egger_p0.001.csv",
  "IVW_p0.00001.csv",
  "IVW_p0.0001.csv",
  "IVW_p0.001.csv"
)

# 为每个文件生成QQ plot
for (file in files) {
  cat("\n===== Processing:", file, "=====\n")
  data <- fread(file)
  data <- data %>%
    filter(!is.na(pvalue), pvalue > 0, pvalue <= 1)
  
  if (nrow(data) == 0) {
    cat("Warning: No valid data in", file, "- skipping\n")
    next
  }
  
  method <- ifelse(grepl("Egger", file), "MR_Egger", "IVW")
  p_threshold <- gsub(".*_(p[0-9.]+)\\.csv", "\\1", file)
  
  chisq <- qchisq(1 - data$pvalue, 1)
  lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
  cat("Lambda:", round(lambda, 3), "\n")
  qq_data <- data %>%
    arrange(pvalue) %>%
    mutate(
      observed = -log10(pvalue),
      expected = -log10(ppoints(n()))
    )
  
  N <- nrow(qq_data)
  qq_data <- qq_data %>%
    mutate(
      clower = -log10(qbeta(0.025, 1:N, N:1)),
      cupper = -log10(qbeta(0.975, 1:N, N:1))
    )
  
  qq_data <- qq_data %>%
    filter(is.finite(observed), is.finite(expected))
  if (nrow(qq_data) == 0) {
    cat("Warning: No valid QQ data - skipping\n")
    next
  }
  y_max <- max(qq_data$observed, na.rm = TRUE)
  x_max <- max(qq_data$expected, na.rm = TRUE)
  lambda_x <- x_max * 0.6
  lambda_y <- y_max * 0.90
  
  p <- ggplot(qq_data, aes(x = expected, y = observed)) +
    geom_ribbon(aes(ymin = clower, ymax = cupper), 
                fill = "grey70", alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, 
                color = "black", linewidth = 0.8) +
    geom_point(color = "#3182bd", alpha = 0.6, size = 2, shape = 1) +
    annotate("text", 
             x = lambda_x, 
             y = lambda_y,
             label = paste("Lambda =", round(lambda, 3)),
             hjust = 0, 
             vjust = 1,
             size = 5,
             alpha = 0.8) +
    labs(
      title = paste0("Q-Q Plot using ", method, " (", p_threshold, ")"),
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
    )
  output_file <- paste0(output_dir, "qq_plot_", gsub("\\.csv", "", file), ".png")
  
  tryCatch({
    png(output_file, width = 7, height = 7, units = "in", res = 300)
    print(p)
    dev.off()
    cat("Plot saved to:", output_file, "\n")
  }, error = function(e) {
    cat("Error saving plot:", e$message, "\n")
    dev.off()
  })
}

cat("\n===== All QQ plots generated! =====\n")