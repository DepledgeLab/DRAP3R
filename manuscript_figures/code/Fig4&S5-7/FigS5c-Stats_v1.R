
library(dplyr)
library(ggplot2)
library(patchwork)  # for arranging plots in a grid

# --- 1. Read data ---
data <- read.csv("pseU-stats.txt", check.names = FALSE, sep = "\t")

data$`ARPE19 #2 poly(A)` <- NULL
data$`CRO-AP5` <- NULL

# Columns to compare vs IVT MIX
conditions <- setdiff(names(data), c("range_start", "IVT MIX"))

# Define z-score cutoff
z_cutoff <- 2

plot_list <- list()

for(cond in conditions) {
  
  df <- data
  
  # --- 1. Compute per-bin differences ---
  df$diff <- df[[cond]] - df$`IVT MIX`
  
  # --- 2. Compute z-scores ---
  df$z <- (df$diff - mean(df$diff)) / sd(df$diff)
  
  # --- 3. Flag significant bins ---
  df$significant <- abs(df$z) > z_cutoff
  
  # --- 4. Plot with z-score cutoff lines ---
  p <- ggplot(df, aes(x = range_start, y = diff, fill = significant)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("FALSE" = "grey", "TRUE" = "deeppink")) +
    geom_hline(yintercept = mean(df$diff) + z_cutoff*sd(df$diff), linetype = "dashed", color = "navy") +
    geom_vline(xintercept = 0.985, linetype = "dotted", color = "deeppink", size = 1) +
    theme_classic() +
    ylim(-0.005, 0.02) +
    labs(title = paste(cond, "vs IVT MIX"),
         x = "modification probability",
         y = "z-score of difference") +
    theme(axis.text = element_text(color = "black"),
        text = element_text(size = 24, color = "black"),
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.position = "none",
        plot.title = element_text(size = 12, hjust = 0.5)) 
  
  plot_list[[cond]] <- p
}

# --- 5. Arrange plots in a grid and save PDF ---
pdf("histogram_differences_zscore_with_thresholds.pdf", width = 5, height = 20)
wrap_plots(plot_list, ncol = 1)
dev.off()