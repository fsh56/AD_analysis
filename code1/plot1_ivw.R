# 分层分析脚本：QQ plots 和分布图
# Stratified analysis by p_threshold

library(ggplot2)
library(dplyr)
library(gridExtra)

# 读取数据
# data <- read.csv("your_data.csv")  # 替换为你的数据文件路径
# 或者
# data <- read.table("your_data.txt", header = TRUE, sep = "\t")

# ============================================================================
# 1. QQ Plot 函数（类似于MR_Egger风格）
# ============================================================================

create_qq_plot <- function(data, p_threshold_value, lambda = NULL) {
  # 筛选指定p_threshold的数据
  df <- data %>% filter(p_threshold == p_threshold_value)
  
  # 去除NA值
  pvalues <- df$pvalue[!is.na(df$pvalue)]
  
  if(length(pvalues) == 0) {
    warning(paste("No valid p-values for p_threshold =", p_threshold_value))
    return(NULL)
  }
  
  # 计算观察值和期望值
  n <- length(pvalues)
  observed <- -log10(sort(pvalues))
  expected <- -log10(ppoints(n))
  
  # 计算lambda（基因组膨胀因子）
  if(is.null(lambda)) {
    chisq <- qchisq(1 - pvalues, 1)
    lambda <- median(chisq) / qchisq(0.5, 1)
  }
  
  # 创建数据框
  qq_data <- data.frame(
    expected = expected,
    observed = observed
  )
  
  # 计算95%置信区间
  # 使用beta分布计算置信区间
  c95 <- data.frame(
    expected = expected,
    upper = -log10(qbeta(0.025, 1:n, n:1)),
    lower = -log10(qbeta(0.975, 1:n, n:1))
  )
  
  # 创建QQ plot
  p <- ggplot(qq_data, aes(x = expected, y = observed)) +
    geom_ribbon(data = c95, aes(x = expected, ymin = lower, ymax = upper),
                fill = "gray80", alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) +
    geom_point(color = "#4A90E2", size = 2, alpha = 0.6) +
    labs(
      title = paste0("Q-Q Plot using IVW"),
      subtitle = paste0("p_threshold = ", p_threshold_value),
      x = "expected -log_10(p)",
      y = "observed -log_10(p)"
    ) +
    annotate("text", x = Inf, y = -Inf, 
             label = paste("Lambda =", round(lambda, 3)),
             hjust = 1.1, vjust = -1, size = 5) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    coord_cartesian(clip = "off")
  
  return(p)
}

# ============================================================================
# 2. n_snps 分布直方图
# ============================================================================

create_nsnps_histogram <- function(data, p_threshold_value) {
  df <- data %>% filter(p_threshold == p_threshold_value)
  
  p <- ggplot(df, aes(x = n_snps)) +
    geom_histogram(binwidth = 1, fill = "#4A90E2", color = "white", alpha = 0.7) +
    labs(
      title = paste0("Distribution of n_snps"),
      subtitle = paste0("p_threshold = ", p_threshold_value),
      x = "Number of SNPs",
      y = "Frequency"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    geom_vline(aes(xintercept = median(n_snps, na.rm = TRUE)),
               color = "red", linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("Median = ", median(df$n_snps, na.rm = TRUE),
                            "\nMean = ", round(mean(df$n_snps, na.rm = TRUE), 1)),
             hjust = 1.1, vjust = 1.5, size = 4, color = "red")
  
  return(p)
}

# ============================================================================
# 3. heterogeneity_Q 分布图
# ============================================================================

create_hetQ_distribution <- function(data, p_threshold_value) {
  df <- data %>% filter(p_threshold == p_threshold_value)
  
  # 创建组合图：直方图 + 密度图
  p <- ggplot(df, aes(x = heterogeneity_Q)) +
    geom_histogram(aes(y = after_stat(density)), 
                   bins = 30, fill = "#E27D60", color = "white", alpha = 0.7) +
    geom_density(color = "#C38D9E", linewidth = 1.2) +
    labs(
      title = paste0("Distribution of Heterogeneity Q"),
      subtitle = paste0("p_threshold = ", p_threshold_value),
      x = "Heterogeneity Q",
      y = "Density"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    geom_vline(aes(xintercept = median(heterogeneity_Q, na.rm = TRUE)),
               color = "red", linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = Inf, y = Inf,
             label = paste0("Median = ", round(median(df$heterogeneity_Q, na.rm = TRUE), 2),
                            "\nMean = ", round(mean(df$heterogeneity_Q, na.rm = TRUE), 2)),
             hjust = 1.1, vjust = 1.5, size = 4, color = "red")
  
  return(p)
}

# ============================================================================
# 4. 主函数：对所有p_threshold进行分析
# ============================================================================

stratified_analysis <- function(data) {
  # 定义p_threshold值
  p_thresholds <- c(0.001, 0.0001)
  
  # 存储所有图形
  qq_plots <- list()
  nsnps_plots <- list()
  hetQ_plots <- list()
  
  # 对每个p_threshold生成图形
  for (pt in p_thresholds) {
    cat(paste0("\n处理 p_threshold = ", pt, "...\n"))
    
    # QQ plot
    qq_plots[[as.character(pt)]] <- create_qq_plot(data, pt)
    
    # n_snps histogram
    nsnps_plots[[as.character(pt)]] <- create_nsnps_histogram(data, pt)
    
    # heterogeneity_Q distribution
    hetQ_plots[[as.character(pt)]] <- create_hetQ_distribution(data, pt)
  }
  
  # 组合并保存QQ plots
  qq_combined <- grid.arrange(grobs = qq_plots, ncol = 3)
  ggsave("qq_plots_combined.png", qq_combined, 
         width = 18, height = 6, dpi = 300)
  
  # 组合并保存n_snps histograms
  nsnps_combined <- grid.arrange(grobs = nsnps_plots, ncol = 3)
  ggsave("nsnps_histograms_combined.png", nsnps_combined, 
         width = 18, height = 6, dpi = 300)
  
  # 组合并保存heterogeneity_Q distributions
  hetQ_combined <- grid.arrange(grobs = hetQ_plots, ncol = 3)
  ggsave("hetQ_distributions_combined.png", hetQ_combined, 
         width = 18, height = 6, dpi = 300)
  
  # 也可以单独保存每个图
  for (pt in p_thresholds) {
    pt_str <- as.character(pt)
    ggsave(paste0("qq_plot_", pt_str, ".png"), qq_plots[[pt_str]], 
           width = 6, height = 6, dpi = 300)
    ggsave(paste0("nsnps_hist_", pt_str, ".png"), nsnps_plots[[pt_str]], 
           width = 6, height = 6, dpi = 300)
    ggsave(paste0("hetQ_dist_", pt_str, ".png"), hetQ_plots[[pt_str]], 
           width = 6, height = 6, dpi = 300)
  }
  
  cat("\n所有图形已保存！\n")
  
  return(list(
    qq_plots = qq_plots,
    nsnps_plots = nsnps_plots,
    hetQ_plots = hetQ_plots
  ))
}

# ============================================================================
# 5. 使用示例
# ============================================================================

# 读取你的数据
data <- read.csv("/Users/fengsihao/Desktop/amyloid_Brain_Frontal_Cortex_BA9_1.txt.gz/all_IVW_combined.csv")

# 运行分层分析
results <- stratified_analysis(data)

# 查看单个图形
results$qq_plots[["0.0001"]]
results$nsnps_plots[["0.001"]]
results$hetQ_plots[["0.0001"]]

# ============================================================================
# 6. 额外：生成描述性统计表
# ============================================================================

generate_summary_table <- function(data) {
  summary_table <- data %>%
    group_by(p_threshold) %>%
    summarise(
      n_genes = n(),
      mean_nsnps = mean(n_snps, na.rm = TRUE),
      median_nsnps = median(n_snps, na.rm = TRUE),
      sd_nsnps = sd(n_snps, na.rm = TRUE),
      mean_hetQ = mean(heterogeneity_Q, na.rm = TRUE),
      median_hetQ = median(heterogeneity_Q, na.rm = TRUE),
      sd_hetQ = sd(heterogeneity_Q, na.rm = TRUE),
      mean_pvalue = mean(pvalue, na.rm = TRUE),
      median_pvalue = median(pvalue, na.rm = TRUE),
      n_significant = sum(pvalue < 0.05, na.rm = TRUE),
      .groups = "drop"
    )
  
  print(summary_table)
  write.csv(summary_table, "summary_statistics.csv", row.names = FALSE)
  
  return(summary_table)
}

# 使用示例：
# summary_stats <- generate_summary_table(data)