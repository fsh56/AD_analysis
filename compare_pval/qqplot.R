# 简洁版QQ Plot：比较不同p_threshold的inflation factor
setwd("/Users/fengsihao/AD_analysis/compare_pval")

library(ggplot2)
library(dplyr)

# ============================================================================
# 简化的QQ Plot函数
# ============================================================================

plot_qq_comparison <- function(method, p_thresholds = c("0.0005", "0.0001")) {
  
  # 读取数据并计算lambda
  qq_data <- data.frame()
  lambda_values <- c()
  
  for (pt in p_thresholds) {
    # 构建文件名（根据实际文件结构）
    if (method == "WeightedMedian" && pt == "0.001") {
      # 特殊情况：WeightedMedian的0.001阈值文件名不同
      file <- paste0("amyloid_BA9_", method, "_results_p_", pt, "_all.csv")
    } else {
      file <- paste0("amyloid_BA9_", method, "_results_p", pt, "_all.csv")
    }
    
    if (!file.exists(file)) {
      warning(paste("文件不存在:", file))
      next
    }
    
    data <- read.csv(file)
    pval_col <- intersect(names(data), c("pvalue", "pval", "p", "P", "p.value"))[1]
    
    if (is.na(pval_col)) {
      warning(paste("在文件", file, "中找不到p-value列"))
      next
    }
    
    pvalues <- data[[pval_col]][!is.na(data[[pval_col]])]
    
    if (length(pvalues) == 0) {
      warning(paste("文件", file, "中没有有效的p-value"))
      next
    }
    
    # 计算QQ数据
    n <- length(pvalues)
    observed <- -log10(sort(pvalues))
    expected <- -log10(ppoints(n))
    
    # 计算lambda (inflation factor)
    chisq <- qchisq(1 - pvalues, 1)
    lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
    
    # 存储数据
    qq_data <- rbind(qq_data, data.frame(
      expected = expected,
      observed = observed,
      threshold = paste0("p < ", pt),
      lambda = round(lambda, 3)
    ))
    
    lambda_values <- c(lambda_values, paste0("p < ", pt, ": λ = ", round(lambda, 3)))
  }
  
  if (nrow(qq_data) == 0) {
    stop(paste("没有找到", method, "方法的任何有效数据"))
  }
  
  # 计算95%置信区间
  n_max <- max(table(qq_data$threshold))
  ci_data <- data.frame(
    expected = -log10(ppoints(n_max)),
    upper = -log10(qbeta(0.025, 1:n_max, n_max:1)),
    lower = -log10(qbeta(0.975, 1:n_max, n_max:1))
  )
  
  # 绘图
  p <- ggplot() +
    geom_ribbon(data = ci_data, aes(x = expected, ymin = lower, ymax = upper),
                fill = "gray80", alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) +
    geom_point(data = qq_data, aes(x = expected, y = observed, color = threshold),
               size = 2, alpha = 0.7) +
    scale_color_manual(values = c("#4A90E2", "#9B59B6")) +
    labs(title = paste0("Q-Q Plot: ", method, " Method"),
         x = "Expected -log10(p)",
         y = "Observed -log10(p)",
         color = "P-value Threshold") +
    annotate("text", x = 0, y = Inf, 
             label = paste(lambda_values, collapse = "\n"),
             hjust = -0.1, vjust = 1.5, size = 3.5) +
    theme_classic() +
    theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          legend.position = "right",
          panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
  
  return(p)
}

# ============================================================================
# 生成所有方法的QQ plots
# ============================================================================

methods <- c("IVW", "Egger", "WeightedMedian")

for (method in methods) {
  cat("生成", method, "的QQ plot...\n")
  tryCatch({
    p <- plot_qq_comparison(method)
    ggsave(paste0("qq_inflation_", method, ".png"), p, 
           width = 8, height = 6, dpi = 300)
    cat("  ✓ 成功生成 qq_inflation_", method, ".png\n", sep = "")
  }, error = function(e) {
    cat("  ✗ 生成", method, "失败:", e$message, "\n")
  })
}

cat("\n完成！\n")