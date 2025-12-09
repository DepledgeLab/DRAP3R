library(ggplot2)
library(reshape2)
library(dplyr)

#setwd("")

data <- read.delim("pseU-over10per-polIIItx-summary.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# Convert data to long format
data_long <- melt(data, id.vars = c("name", "position"), variable.name = "Condition", value.name = "Expression")

# Create the heatmap
pseU<-ggplot(data_long, aes(x = Condition, y = paste(name, position, sep = " "), fill = Expression)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkblue", limits = c(0,100)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16, color = "black"), 
    axis.text.y = element_text(size = 16, color = "black"), 
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(), 
    plot.title = element_text(size = 18, face = "bold", color = "black"),  
    legend.text = element_text(size = 16, color = "black"),  
    legend.title = element_text(size = 16, face = "bold", color = "black") 
  )

pdf(file = "pseU-heatmap.pdf", width = 6, height = 8)
pseU
dev.off()
