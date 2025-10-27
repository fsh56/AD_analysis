# 对比QQ Plots：不同方法和不同p_threshold
# Comparative QQ Plots for MR Results
setwd("/Users/fengsihao/AD_analysis/mr_results1022/amyloid_ba9/mr_results")

library(ggplot2)
library(dplyr)
library(gridExtra)

# ============================================================================
# 1. QQ Plot 函数（支持多个数据集叠加）
# ============================================================================

create_comparative_qq_plot <- function(data_list, labels, colors, title, lambda_pos = "topleft") {
  # data_list: 包含多个pvalue向量的列表
  # labels: 每个数据集的标签
  # colors: 每个数据集的颜色
  
  # 准备所有数据
  all_qq_data <- data.frame()
  lambda_text <- c()
  
  for (i in 1:length(data_list)) {
    pvalues <- data_list[[i]][!is.na(data_list[[i]])]
    
    if (length(pvalues) == 0) {
      warning(paste("No valid p-values for", labels[i]))
      next
    }
    
    # 计算QQ plot数据
    n <- length(pvalues)
    observed <- -log10(sort(pvalues))
    expected <- -log10(ppoints(n))
    
    # 计算lambda
    chisq <- qchisq(1 - pvalues, 1)
    lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
    
    # 存储数据
    qq_data <- data.frame(
      expected = expected,
      observed = observed,
      method = labels[i],
      color = colors[i]
    )
    
    all_qq_data <- rbind(all_qq_data, qq_data)
    lambda_text <- c(lambda_text, paste0(labels[i], ": λ = ", round(lambda, 3)))
  }
  
  # 计算95%置信区间（使用第一个数据集的n）
  n_max <- max(table(all_qq_data$method))
  expected_ci <- -log10(ppoints(n_max))
  c95 <- data.frame(
    expected = expected_ci,
    upper = -log10(qbeta(0.025, 1:n_max, n_max:1)),
    lower = -log10(qbeta(0.975, 1:n_max, n_max:1))
  )
  
  # 创建图形
  p <- ggplot() +
    # 95% 置信区间
    geom_ribbon(data = c95, aes(x = expected, ymin = lower, ymax = upper),
                fill = "gray80", alpha = 0.5) +
    # 对角参考线
    geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) +
    # 数据点
    geom_point(data = all_qq_data, 
               aes(x = expected, y = observed, color = method),
               size = 2, alpha = 0.6) +
    # 手动设置颜色
    scale_color_manual(values = setNames(colors, labels)) +
    labs(
      title = title,
      x = "expected -log_10(p)",
      y = "observed -log_10(p)",
      color = "Method"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    coord_cartesian(clip = "off")
  
  # 添加Lambda文本
  lambda_label <- paste(lambda_text, collapse = "\n")
  
  if (lambda_pos == "topleft") {
    p <- p + annotate("text", x = 0, y = Inf, 
                      label = lambda_label,
                      hjust = -0.1, vjust = 1.5, size = 4)
  } else {
    p <- p + annotate("text", x = Inf, y = -Inf, 
                      label = lambda_label,
                      hjust = 1.1, vjust = -0.5, size = 4)
  }
  
  return(p)
}

# ============================================================================
# 2. 读取数据的辅助函数
# ============================================================================

read_mr_results <- function(method, p_threshold) {
  # 构建文件名
  filename <- paste0("amyloid_", method, "_p", p_threshold, ".csv")
  
  if (!file.exists(filename)) {
    warning(paste("File not found:", filename))
    return(NULL)
  }
  
  data <- read.csv(filename)
  
  # 假设pvalue列名为 "pvalue" 或 "pval" 或 "p"
  pval_col <- intersect(names(data), c("pvalue", "pval", "p", "P", "p.value"))
  
  if (length(pval_col) == 0) {
    warning(paste("No p-value column found in", filename))
    return(NULL)
  }
  
  return(data[[pval_col[1]]])
}

# ============================================================================
# 3. 生成对比图 - 同一p_threshold，不同方法
# ============================================================================

compare_methods_by_pthreshold <- function(p_threshold) {
  # 读取三种方法的数据
  ivw_pval <- read_mr_results("IVW", p_threshold)
  egger_pval <- read_mr_results("Egger", p_threshold)
  median_pval <- read_mr_results("WeightedMedian", p_threshold)
  
  # 准备数据列表
  data_list <- list(ivw_pval, egger_pval, median_pval)
  labels <- c("IVW", "MR-Egger", "Weighted Median")
  colors <- c("#4A90E2", "#E27D60", "#85C88A")  # 蓝色、橙红色、绿色
  
  # 创建标题
  title <- paste0("Q-Q Plot Comparing Methods (p_threshold = ", p_threshold, ")")
  
  # 创建图形
  p <- create_comparative_qq_plot(data_list, labels, colors, title)
  
  return(p)
}

# ============================================================================
# 4. 生成对比图 - 同一方法，不同p_threshold
# ============================================================================

compare_pthresholds_by_method <- function(method) {
  # 读取不同p_threshold的数据
  p0001_pval <- read_mr_results(method, "0.001")
  p00001_pval <- read_mr_results(method, "0.0001")
  
  # 准备数据列表
  data_list <- list(p0001_pval, p00001_pval)
  labels <- c("p = 0.001", "p = 0.0001")
  colors <- c("#4A90E2", "#9B59B6")  # 蓝色、紫色
  
  # 创建标题
  title <- paste0("Q-Q Plot using ", method)
  
  # 创建图形
  p <- create_comparative_qq_plot(data_list, labels, colors, title)
  
  return(p)
}

# ============================================================================
# 5. 生成所有对比图
# ============================================================================

generate_all_comparative_plots <- function() {
  cat("开始生成对比QQ plots...\n\n")
  
  # ========================================
  # 部分1: 同一p_threshold，不同方法
  # ========================================
  cat("生成部分1: 同一p_threshold下的方法对比\n")
  
  p_thresholds <- c("0.001", "0.0001")
  method_comparison_plots <- list()
  
  for (pt in p_thresholds) {
    cat(paste0("  处理 p_threshold = ", pt, "...\n"))
    p <- compare_methods_by_pthreshold(pt)
    method_comparison_plots[[pt]] <- p
    
    # 单独保存
    ggsave(paste0("qq_methods_comparison_p", pt, ".png"), p,
           width = 8, height = 6, dpi = 300)
  }
  
  # 组合保存
  combined_methods <- grid.arrange(grobs = method_comparison_plots, ncol = 2)
  ggsave("qq_all_methods_comparison.png", combined_methods,
         width = 16, height = 6, dpi = 300)
  
  # ========================================
  # 部分2: 同一方法，不同p_threshold
  # ========================================
  cat("\n生成部分2: 同一方法下的p_threshold对比\n")
  
  methods <- c("IVW", "Egger", "WeightedMedian")
  pthreshold_comparison_plots <- list()
  
  for (method in methods) {
    cat(paste0("  处理方法: ", method, "...\n"))
    p <- compare_pthresholds_by_method(method)
    pthreshold_comparison_plots[[method]] <- p
    
    # 单独保存
    ggsave(paste0("qq_pthreshold_comparison_", method, ".png"), p,
           width = 8, height = 6, dpi = 300)
  }
  
  # 组合保存
  combined_pthresholds <- grid.arrange(grobs = pthreshold_comparison_plots, ncol = 3)
  ggsave("qq_all_pthreshold_comparison.png", combined_pthresholds,
         width = 18, height = 6, dpi = 300)
  
  cat("\n所有QQ plots已生成！\n")
  cat("\n生成的文件:\n")
  cat("========================================\n")
  cat("方法对比 (同一p_threshold):\n")
  cat("  - qq_methods_comparison_p0.001.png\n")
  cat("  - qq_methods_comparison_p0.0001.png\n")
  cat("  - qq_all_methods_comparison.png (组合)\n")
  cat("\np_threshold对比 (同一方法):\n")
  cat("  - qq_pthreshold_comparison_IVW.png\n")
  cat("  - qq_pthreshold_comparison_Egger.png\n")
  cat("  - qq_pthreshold_comparison_WeightedMedian.png\n")
  cat("  - qq_all_pthreshold_comparison.png (组合)\n")
  cat("========================================\n")
  
  return(list(
    method_comparison = method_comparison_plots,
    pthreshold_comparison = pthreshold_comparison_plots
  ))
}

# ============================================================================
# 6. 使用示例
# ============================================================================

# 运行完整分析
results <- generate_all_comparative_plots()

# 查看单个图形
# results$method_comparison[["0.001"]]
# results$pthreshold_comparison[["IVW"]]

# ============================================================================
# 7. 自定义单个对比图
# ============================================================================

# 示例：自定义p_threshold = 0.001的方法对比
# custom_plot <- compare_methods_by_pthreshold("0.001")
# ggsave("custom_methods_comparison.png", custom_plot, width = 8, height = 6, dpi = 300)

# 示例：自定义IVW方法的p_threshold对比
# custom_plot2 <- compare_pthresholds_by_method("IVW")
# ggsave("custom_pthreshold_comparison.png", custom_plot2, width = 8, height = 6, dpi = 300)

# ============================================================================
# 8. 额外功能：三合一QQ plot（所有方法和p_threshold）
# ============================================================================

create_comprehensive_qq_plot <- function() {
  # 读取所有组合的数据
  all_data <- list()
  all_labels <- c()
  all_colors <- c()
  
  methods <- c("IVW", "Egger", "WeightedMedian")
  p_thresholds <- c("0.001", "0.0001")
  
  color_palette <- list(
    IVW = c("#4A90E2", "#2E5F8E"),
    Egger = c("#E27D60", "#C85A47"),
    WeightedMedian = c("#85C88A", "#5FA368")
  )
  
  for (method in methods) {
    for (i in 1:length(p_thresholds)) {
      pt <- p_thresholds[i]
      pval <- read_mr_results(method, pt)
      if (!is.null(pval)) {
        all_data <- c(all_data, list(pval))
        all_labels <- c(all_labels, paste0(method, " (p=", pt, ")"))
        all_colors <- c(all_colors, color_palette[[method]][i])
      }
    }
  }
  
  p <- create_comparative_qq_plot(
    all_data, 
    all_labels, 
    all_colors, 
    "Comprehensive Q-Q Plot: All Methods and P-thresholds",
    lambda_pos = "topleft"
  )
  
  ggsave("qq_comprehensive_all.png", p, width = 10, height = 8, dpi = 300)
  
  return(p)
}

# 生成综合QQ plot
# comprehensive_plot <- create_comprehensive_qq_plot()