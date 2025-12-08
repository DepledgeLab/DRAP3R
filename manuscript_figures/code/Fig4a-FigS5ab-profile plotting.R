library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)

setwd("")

# Import the tab-delimited file
pseU <- read.delim("pseU.txt", header = TRUE, sep = "\t")
U <- read.delim("U.txt", header = TRUE, sep = "\t")
m6A <- read.delim("m6A.txt", header = TRUE, sep = "\t")
A <- read.delim("A.txt", header = TRUE, sep = "\t")

# Extract only the required columns
plot_pseU <- pseU[, c(2, 4:9,12)]
plot_U <- U[, c(2, 4:9,12)]
plot_m6A <- m6A[, c(2, 4:9,12)]
plot_A <- A[, c(2, 4:9,12)]

# Reshape the data to long format for ggplot2
plot_pseU_long <- melt(plot_pseU, id.vars = colnames(pseU)[2], variable.name = "Series", value.name = "Value")
plot_U_long <- melt(plot_U, id.vars = colnames(U)[2], variable.name = "Series", value.name = "Value")
plot_m6A_long <- melt(plot_m6A, id.vars = colnames(m6A)[2], variable.name = "Series", value.name = "Value")
plot_A_long <- melt(plot_A, id.vars = colnames(A)[2], variable.name = "Series", value.name = "Value")

custom_colors <- c("forestgreen", "cadetblue", "cadetblue", "goldenrod1", "darkorange2", "cadetblue", "grey")


### PLOTTING FOR FIG 4A

# Create the plot
p<-ggplot(plot_pseU_long, aes_string(x = colnames(pseU)[2], y = "Value", color = "Series", group = "Series")) +
  geom_line(size = 0.75, linetype = "solid", show.legend = FALSE) +  
  scale_color_manual(values = custom_colors) +  
  labs(x = expression(psi ~ "probability"), y = expression(psi ~ "fraction")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.1)) +  
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 24, colour = "black")) +
  theme(axis.line = element_line(size = 1),  
        axis.ticks = element_line(size = 1))

u<-ggplot(plot_U_long, aes_string(x = colnames(U)[2], y = "Value", color = "Series", group = "Series")) +
  geom_line(size = 0.75, linetype = "solid", show.legend = FALSE) +  
  scale_color_manual(values = custom_colors) +  
  labs(x = "U probability", y = "U fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.8)) +  
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 24, colour = "black")) +
  theme(axis.line = element_line(size = 1),  
        axis.ticks = element_line(size = 1))

m<-ggplot(plot_m6A_long, aes_string(x = colnames(m6A)[2], y = "Value", color = "Series", group = "Series")) +
  geom_line(size = 0.75, linetype = "solid", show.legend = FALSE) +  
  scale_color_manual(values = custom_colors) +  
  labs(x = expression("m"^6*"A probability"), y = expression("m"^6*"A fraction")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.1)) +  
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 24, colour = "black")) +
  theme(axis.line = element_line(size = 1),  
        axis.ticks = element_line(size = 1))

a<-ggplot(plot_A_long, aes_string(x = colnames(A)[2], y = "Value", color = "Series", group = "Series")) +
  geom_line(size = 0.75, linetype = "solid", show.legend = FALSE) +  
  scale_color_manual(values = custom_colors) +  
  labs(x = "A probability", y = "A fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.8)) +  
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 24, colour = "black")) +
  theme(axis.line = element_line(size = 1),  
        axis.ticks = element_line(size = 1))


combined_plot <- (u | a) / (p | m)

combined_plot

pdf(file = "ProfilePlot.pdf", width = 10, height = 10)
combined_plot
dev.off()


p2<-ggplot(plot_pseU_long, aes_string(x = colnames(pseU)[2], y = "Value", color = "Series", group = "Series")) +
  geom_line(size = 0.75, linetype = "solid", show.legend = FALSE) +  
  scale_color_manual(values = custom_colors) +  
  labs(x = expression(psi ~ "probability"), y = expression(psi ~ "fraction")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.05)) +  
  scale_x_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0.9, 1.0)) + 
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 24, colour = "black")) +
  theme(axis.line = element_line(size = 1),  
        axis.ticks = element_line(size = 1))

m2<-ggplot(plot_m6A_long, aes_string(x = colnames(m6A)[2], y = "Value", color = "Series", group = "Series")) +
  geom_line(size = 0.75, linetype = "solid", show.legend = FALSE) +  
  scale_color_manual(values = custom_colors) +  
  labs(x = expression("m"^6*"A probability"), y = expression("m"^6*"A fraction")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.05)) +  
  scale_x_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0.9, 1.0)) + 
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 24, colour = "black")) +
  theme(axis.line = element_line(size = 1),  
        axis.ticks = element_line(size = 1))

combined_plot2 <- (p2 | m2)

combined_plot2

pdf(file = "ProfilePlotCloseup.pdf", width = 10, height = 5)
combined_plot2
dev.off()

# Add a column to distinguish the datasets
plot_pseU_long$Dataset <- "pseU"
plot_m6A_long$Dataset <- "m6A"

# Combine the two datasets by binding rows
combined_data_mods <- rbind(plot_pseU_long, plot_m6A_long)

# Get the unique series names from the combined data
series_listA <- unique(combined_data_mods$Series)

# Custom colors for the two datasets (e.g., A and a)
custom_colorsA <- c("pseU" = "#5979bc", "m6A" = "red")

# Create individual plots for each "Series" in the combined dataset
plots <- list()
for (series in series_listA) {
  # Subset the data for each series
  plot_data <- subset(combined_data_mods, Series == series)
  
  # Create the plot for this series
  p <- ggplot(plot_data, aes_string(x = colnames(A)[2], y = "Value", fill = "Dataset")) +
    geom_bar(stat = "identity", position = "identity", alpha = 0.6, show.legend = FALSE) +
    scale_fill_manual(values = custom_colorsA) +
    labs(x = "modification probability", y = "fraction", title = series) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01), limits = c(0, 0.08)) +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          text = element_text(size = 24, color = "black"),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(size = 1),
          legend.position = "none",
          plot.title = element_text(size = 12, hjust = 0.5)) # Make title text smaller and center it
  
  plots[[series]] <- p  # Store each plot in the list
}

# Combine the individual plots into a 2x4 grid using patchwork
combined_plotA <- wrap_plots(plots, ncol = 2, nrow = 4)
combined_plotA

pdf(file = "SepProfilePlot_pseU_m6A.pdf", width = 10, height = 20)
combined_plotA
dev.off()



### PLOTTING FOR FIG S5A&B

# Add a column to distinguish the datasets
plot_U_long$Dataset <- "U"
plot_A_long$Dataset <- "A"

# Combine the two datasets by binding rows
combined_data_norm <- rbind(plot_U_long, plot_A_long)

combined_data_norm <- combined_data_norm %>%
  arrange(Series, factor(Dataset, levels = c("U", "A")))

# Get the unique series names from the combined data
series_listU <- unique(combined_data_norm$Series)

# Custom colors for the two datasets (e.g., A and a)
custom_colorsU <- c("U" = "black", "A" = "violetred3")

# Create individual plots for each "Series" in the combined dataset
plots <- list()
for (series in series_listU) {
  # Subset the data for each series
  plot_data2 <- subset(combined_data_norm, Series == series)
  
  # Create the plot for this series
  p <- ggplot(plot_data2, aes_string(x = colnames(A)[2], y = "Value", fill = "Dataset")) +
    geom_bar(stat = "identity", position = "identity", alpha = 0.8, show.legend = FALSE) +
    scale_fill_manual(values = custom_colorsU) +
    labs(x = "modification probability", y = "fraction", title = series) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.01), limits = c(0, 0.8)) +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          text = element_text(size = 24, color = "black"),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(size = 1),
          legend.position = "none",
          plot.title = element_text(size = 12, hjust = 0.5))  # Make title text smaller and center it
  
  plots[[series]] <- p  # Store each plot in the list
}

# Combine the individual plots into a 2x4 grid using patchwork
combined_plotU <- wrap_plots(plots, ncol = 2, nrow = 4)
combined_plotU

pdf(file = "SepProfilePlot_U_A.pdf", width = 10, height = 20)
combined_plotU
dev.off()









