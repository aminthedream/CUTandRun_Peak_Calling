# Benchmarking Peak Calling Method for CUT&RUN
# Author: Elias Orouji
# Date: November 6, 2024
# Manuscript Title: "Benchmarking Peak Calling Methods for CUT&RUN"
# Figure Number: [Figure 5]

# Description: This script generates the plots for precision, recall and F1 score for all four peak calling methods based on the histone mark tested.

# Notes:
# - This script is intended for use with data from mouse brain samples with three different histone marks.
# - The peak calling methods compared include MACS2, SEACR, GoPeaks, and LanceOtron.
# - Ensure the appropriate input files are loaded for each figure, and customize the parameters as needed.

# Clean up environment
rm(list = ls())

library(tidyverse)
library(readxl)
library(ggpubr)

# Set working directory
setwd("~/Desktop/R_scripts/")

# Load the data
file_path <- "all_performance_metrics.xlsx"
data <- read_excel(file_path)

# Check data structure
head(data)

# Rename columns to match expected names
data <- data %>% rename(Score = F1)

# Split the data for each metric and histone mark
data_f1 <- data %>% select(Method, Score, HistoneMark) %>% mutate(Metric = "F1")
data_precision <- data %>% select(Method, Precision, HistoneMark) %>% rename(Score = Precision) %>% mutate(Metric = "Precision")
data_recall <- data %>% select(Method, Recall, HistoneMark) %>% rename(Score = Recall) %>% mutate(Metric = "Recall")

# Function to create plot for each metric with statistical tests, divided by histone mark
plot_metric_by_histone <- function(data, metric_name) {
  ggboxplot(data, x = "Method", y = "Score", fill = "Method", palette = "jco", facet.by = "HistoneMark") +
    stat_compare_means(method = "anova", label = "p.format", label.y = max(data$Score, na.rm = TRUE) * 1.1) +
    stat_compare_means(comparisons = list(c("GOPEAKS", "LANCEOTRON"), c("GOPEAKS", "MACS2"), c("GOPEAKS", "SEACR"),
                                          c("LANCEOTRON", "MACS2"), c("LANCEOTRON", "SEACR"), c("MACS2", "SEACR")),
                       method = "t.test", label = "p.signif", na.rm = TRUE, label.y.npc = "top", label.size = 3) +
    labs(title = paste("Average", metric_name, "Score by Method and Histone Mark"), y = paste(metric_name, "Score"), x = "Method") +
    theme_minimal()
}

# Plot each metric divided by histone mark
plot_f1_histone <- plot_metric_by_histone(data_f1, "F1")
plot_precision_histone <- plot_metric_by_histone(data_precision, "Precision")
plot_recall_histone <- plot_metric_by_histone(data_recall, "Recall")

# Arrange the plots
ggarrange(plot_f1_histone, plot_precision_histone, plot_recall_histone, ncol = 1, nrow = 3)

# Save the plots to a file
ggsave("performance_metrics_comparison_by_histone.png", width = 12, height = 18)

