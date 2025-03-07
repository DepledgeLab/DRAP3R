library(ggplot2)
library(reshape2)
library(dplyr)

data <- read.delim("PolIIItx_heatmap_data/m6A-over10per-polIIItx-summary.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# Convert data to long format
data_long <- melt(data, id.vars = c("name", "position"), variable.name = "Condition", value.name = "Expression")

# Create the heatmap
m6A<-ggplot(data_long, aes(x = Condition, y = paste(name, position, sep = " "), fill = Expression)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkred", limits = c(0,100)) + # Adjust colors as needed
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16, color = "black"),  # Increase x-axis label size & make black
    axis.text.y = element_text(size = 16, color = "black"),  # Increase y-axis label size & make black
    axis.title.x = element_blank(),  # Remove x-axis title
    axis.title.y = element_blank(),  # Remove y-axis title
    plot.title = element_text(size = 18, face = "bold", color = "black"),  # Increase title size & make black
    legend.text = element_text(size = 16, color = "black"),  # Increase legend text size & make black
    legend.title = element_text(size = 16, face = "bold", color = "black")  # Increase legend title size & make black
  )

pdf(file = "m6A-heatmap.pdf", width = 5, height = 5)
m6A
dev.off()