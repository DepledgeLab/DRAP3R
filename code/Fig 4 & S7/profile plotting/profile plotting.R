library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)

# Import the tab-delimited file
pseU <- read.delim("pseU.txt", header = TRUE, sep = "\t")
U <- read.delim("U.txt", header = TRUE, sep = "\t")
m6A <- read.delim("m6A.txt", header = TRUE, sep = "\t")
A <- read.delim("A.txt", header = TRUE, sep = "\t")



# Extract only the required columns
plot_pseU <- pseU[, c(2, 4:7)]
plot_U <- U[, c(2, 4:7)]
plot_m6A <- m6A[, c(2, 4:7)]
plot_A <- A[, c(2, 4:7)]

# Reshape the data to long format for ggplot2
plot_pseU_long <- melt(plot_pseU, id.vars = colnames(pseU)[2], variable.name = "Series", value.name = "Value")
plot_U_long <- melt(plot_U, id.vars = colnames(U)[2], variable.name = "Series", value.name = "Value")
plot_m6A_long <- melt(plot_m6A, id.vars = colnames(m6A)[2], variable.name = "Series", value.name = "Value")
plot_A_long <- melt(plot_A, id.vars = colnames(A)[2], variable.name = "Series", value.name = "Value")




custom_colors <- c("deeppink", "deeppink", "deeppink4", "cadetblue")



# Create the plot
p<-ggplot(plot_pseU_long, aes_string(x = colnames(pseU)[2], y = "Value", color = "Series", group = "Series")) +
  geom_line(size = 0.75, linetype = "solid", show.legend = FALSE) +  # Make lines thicker and smoother
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  labs(x = "mod probability", y = "pseU fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.1)) +  # Set limits and format y-axis
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 24, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness

u<-ggplot(plot_U_long, aes_string(x = colnames(U)[2], y = "Value", color = "Series", group = "Series")) +
  geom_line(size = 0.75, linetype = "solid", show.legend = FALSE) +  # Make lines thicker and smoother
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  labs(x = "mod probability", y = "U fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.8)) +  # Set limits and format y-axis
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 24, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness

m<-ggplot(plot_m6A_long, aes_string(x = colnames(m6A)[2], y = "Value", color = "Series", group = "Series")) +
  geom_line(size = 0.75, linetype = "solid", show.legend = FALSE) +  # Make lines thicker and smoother
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  labs(x = "mod probability", y = "m6A fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.1)) +  # Set limits and format y-axis
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 24, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness

a<-ggplot(plot_A_long, aes_string(x = colnames(A)[2], y = "Value", color = "Series", group = "Series")) +
  geom_line(size = 0.75, linetype = "solid", show.legend = FALSE) +  # Make lines thicker and smoother
  scale_color_manual(values = custom_colors) +  # Apply custom colors
  labs(x = "mod probability", y = "A fraction") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001), limits = c(0, 0.8)) +  # Set limits and format y-axis
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(text = element_text(size = 24, colour = "black")) +
  theme(axis.line = element_line(size = 1),  # Increase axis line thickness
        axis.ticks = element_line(size = 1)) # Increase tick line thickness


combined_plot <- (u | a) / (p | m)

combined_plot

pdf(file = "HSV1 ProfilePlot.pdf", width = 10, height = 10)
combined_plot
dev.off()



