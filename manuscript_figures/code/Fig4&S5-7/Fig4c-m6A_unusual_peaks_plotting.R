
library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)

setwd("")

# Import the tab-delimited file
m6A <- read.delim("Fig4c-m6A-unusual-peaks.txt", header = TRUE, sep = "\t")

# Extract only the required columns

plot_m6A <- m6A[, c(2, 4:8)]

# Reshape the data to long format for ggplot2
plot_m6A_long <- melt(plot_m6A, id.vars = colnames(m6A)[2], variable.name = "Series", value.name = "Value")


custom_colors <- c("grey", "orchid1", "darkorange3", "thistle3", "navy")


# Create the plot
p<-ggplot(plot_m6A_long, aes_string(x = colnames(m6A)[2], y = "Value", color = "Series", group = "Series")) +
  geom_line(size = 0.75, linetype = "solid", show.legend = FALSE) +  # Make lines thicker and smoother
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  labs(x = expression("m"^6*"A probability"), y = expression("m"^6*"A fraction")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.1)) +  # Set limits and format y-axis
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 24, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness


pdf(file = "m6a-peaks.pdf", width = 5, height = 5)
p
dev.off()

