library(ggplot2)
# Load data (replace with your actual file path)
library(readr)
data <- read_csv("Fstats_Amy50_GeneSummary_p0.001.csv")


# (1) Histogram of n_snps
ggplot(data, aes(x = n_snps)) +
  geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of n_snps", x = "n.snps", y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

# (2) Histogram of mean_F
ggplot(data, aes(x = mean_F)) +
  geom_histogram(bins = 30, fill = "coral", color = "black", alpha = 0.7) +
  labs(title = "Distribution of mean_F", x = "mean_F", y = "Count") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

# (3) Summary statistics for strength_flag
table(data$strength_flag)
prop.table(table(data$strength_flag))

# (4) Summary statistics for pct_weak_lt10
summary(data$pct_weak_lt10)
sd(data$pct_weak_lt10, na.rm = TRUE)
