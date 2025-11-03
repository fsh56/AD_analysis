library(ggplot2)
library(dplyr)
library(gridExtra)

setwd("/Users/fengsihao/AD_analysis/mr_results1027")

# 读取数据
ivw <- read.csv("amyloid_BA9_IVW_results_all.csv")
egger <- read.csv("amyloid_BA9_Egger_results_all.csv")
wm <- read.csv("amyloid_BA9_WeightedMedian_results_all.csv")

n_snps_to_remove <- 2
top_n_to_remove <- 0 
ivw_filtered <- ivw %>% 
  filter(n_snps > n_snps_to_remove) %>%
  arrange(pvalue) %>%
  slice(-(1:top_n_to_remove))  # 去掉p值最小的前N个

egger_filtered <- egger %>% 
  filter(n_snps > n_snps_to_remove) %>%
  arrange(pvalue) %>%
  slice(-(1:top_n_to_remove))

wm_filtered <- wm %>% 
  filter(n_snps > n_snps_to_remove) %>%
  arrange(pvalue) %>%
  slice(-(1:top_n_to_remove))


# 定义QQ plot函数（带95% CI和lambda）
make_qqplot <- function(data, title, color = "#4A90E2") {
  data <- data %>% filter(!is.na(pvalue) & pvalue > 0)
  n <- nrow(data)
  
  # 计算observed和expected
  observed <- -log10(sort(data$pvalue))
  expected <- -log10(ppoints(n))
  
  # 计算lambda
  chisq <- qchisq(1 - data$pvalue, 1)
  lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
  lambda_text <- paste0("λ = ", round(lambda, 3))
  
  # 计算95%置信区间
  c95 <- data.frame(
    expected = expected,
    upper = -log10(qbeta(0.025, 1:n, n:1)),
    lower = -log10(qbeta(0.975, 1:n, n:1))
  )
  
  df <- data.frame(expected = expected, observed = observed)
  
  # 创建图形
  p <- ggplot() +
    # 95% 置信区间
    geom_ribbon(data = c95, aes(x = expected, ymin = lower, ymax = upper),
                fill = "gray80", alpha = 0.5) +
    # 对角参考线
    geom_abline(intercept = 0, slope = 1, color = "black", linewidth = 0.8) +
    # 数据点
    geom_point(data = df, aes(x = expected, y = observed), 
               alpha = 0.6, size = 2, color = color) +
    # Lambda标注
    annotate("text", x = 0, y = Inf, 
             label = lambda_text,
             hjust = -0.1, vjust = 1.5, size = 5) +
    labs(x = "Expected -log10(p)", y = "Observed -log10(p)", title = title) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
    ) +
    coord_cartesian(clip = "off")
  
  return(p)
}

# 画QQ plots
p_ivw <- make_qqplot(ivw_filtered, "QQ Plot - IVW", color = "#4A90E2")
p_egger <- make_qqplot(egger_filtered, "QQ Plot - MR Egger", color = "#E27D60")
p_wm <- make_qqplot(wm_filtered, "QQ Plot - Weighted Median", color = "#85C88A")
# 拼图展示
p_combined <- grid.arrange(p_ivw, p_egger, p_wm, ncol = 3)
ggsave("qqplot_combined.png", p_combined, width = 18, height = 6, dpi = 300)

