all_results <- read.csv("/Users/fengsihao/Desktop/amyloid_Brain_Frontal_Cortex_BA9_1.txt.gz/all_IVW_combined.csv")
all_results <- all_results %>% filter(p_threshold == 0.001)
chisq <- qchisq(1 - all_results$pvalue, 1)
lambda <- median(chisq) / qchisq(0.5, 1)
print(paste("Inflation factor =", round(lambda, 3)))

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
    title = "Q-Q Plot, IVW, p_threshold=0.001",
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
print(p2)