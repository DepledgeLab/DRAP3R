# Load required libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Set working directory
#setwd("")

chromosome_filter <- "Homo_sapiens_tRNA-Glu-TTC-14"

# Function to read and preprocess data
process_variant_file <- function(file_path, chromosome_filter) {
  df <- read.csv(file_path, sep = '\t', header = TRUE)
  df <- df[!grepl("ins", df$Position), ]  # Remove insertions
  df <- df[df$Chromosome == chromosome_filter, ]  # Filter based on chromosome
  return(df)
}

# Process DRAP3R and nano data
DRAP3R <- process_variant_file('ARPE19_uninf-2.sup-pseU_bwa_W13k6T20.SR.merged.primary.splint_tRNAdb.variantCalls.txt',chromosome_filter)
nano <- process_variant_file('ARPE19_nanotRNA_uninf-2.sup-pseU_bwa_W13k6T20.SR.merged.primary.splint_tRNAdb.variantCalls.txt',chromosome_filter)

# Function to calculate max value, accuracy, and additional X column
calculate_max_and_accuracy <- function(df) {
  max_values <- numeric(nrow(df))
  X_values <- numeric(nrow(df))
  
  for (i in 1:nrow(df)) {
    # Calculate sum of ACGTN values
    acgtn_sum <- sum(df[i, c("A", "C", "G", "T", "N")], na.rm = TRUE)
    
    # Check if sum of ACGTN is less than 25 and if Depth is less than 25
    if (acgtn_sum < 25 && df$Depth[i] < 25) {
      X_values[i] <- acgtn_sum
    } else {
      X_values[i] <- 0
    }
    
    # Zero out ACGTN values if X is greater than 0
    if (X_values[i] > 0) {
      df[i, c("A", "C", "G", "T", "N")] <- 0
    }
    
    # Calculate max value and accuracy
    if (df$Consensus[i] == df$Reference[i]) {
      max_values[i] <- max(df[i, c("A", "C", "G", "T", "N")])
    } else {
      exclude_column <- df$Consensus[i]
      columns_to_check <- setdiff(names(df)[8:12], exclude_column)
      values_to_check <- df[i, columns_to_check]
      max_values[i] <- max(values_to_check)
    }
  }
  
  df$MaxValue <- max_values
  df$Accuracy <- df$MaxValue / df$Depth
  df$X <- X_values  # Add the X column to the dataframe
  
  return(df)
}

# Apply the function to calculate MaxValue, Accuracy, and X for DRAP3R and nano
DRAP3R <- calculate_max_and_accuracy(DRAP3R)
nano <- calculate_max_and_accuracy(nano)

DRAP3R$Position <- as.numeric(as.character(DRAP3R$Position))
nano$Position <- as.numeric(as.character(nano$Position))


### PLOTTING
nano_complex<-nano
DRAP3R_complex<-DRAP3R
nano_complex$Source <- "nano-tRNA-Seq"
DRAP3R_complex$Source <- "DRAP3R"

merged_complex <- bind_rows(nano_complex, DRAP3R_complex)

merged_complex <- merged_complex %>% filter(Position >= 24 & Position <= 99)# & grepl("Glu", Chromosome))

merged_complex$Adjusted_Position <- merged_complex$Position - 23

# Calculate medians for each Position and Source
medians_df <- merged_complex %>%
  group_by(Position, Source) %>%
  summarise(Median_Accuracy = median(Accuracy, na.rm = TRUE), .groups = "drop")

medians_df$Adjusted_Position <- medians_df$Position - 23

# Now plot using the 'Adjusted_Position' variable
complex <- ggplot(merged_complex, aes(x = Adjusted_Position, y = Accuracy, color = Source)) +
  # Add individual points
  geom_jitter(width = 0.2, size = 0.1, alpha = 0.6) +
  # Add median lines
  geom_line(data = medians_df, aes(x = Adjusted_Position, y = Median_Accuracy, group = Source, color = Source), size = 1) +
  # Optional: Add points at the median values for clarity
  geom_point(data = medians_df, aes(x = Adjusted_Position, y = Median_Accuracy, color = Source), size = 0, shape = 18) +
  labs(title = "",
       x = "tRNA position",
       y = "Fraction") +
  theme_minimal() +  # Clean theme
  scale_color_manual(values = c("#6FAF60", "#FF9A3D")) +  # Custom colors for each dataset
  scale_x_continuous(breaks = seq(0, max(merged_complex$Adjusted_Position), by = 5))  # Custom x-axis breaks

# Customize the appearance
complex <- complex + theme(
  # Set text color to black
  text = element_text(size = 14, colour = "black"),  # Text size and color
  axis.title = element_text(size = 18, colour = "black"),  # Axis titles
  axis.text = element_text(size = 18, colour = "black"),   # Axis labels
  plot.title = element_text(size = 18, face = "bold", colour = "black"),  # Plot title
  # Remove grid lines
  panel.grid = element_blank(),
  # Add axis lines only on the x and y axes
  axis.line.x = element_line(colour = "black"),  # Bottom axis line
  axis.line.y = element_line(colour = "black"),  # Left axis line
  # Remove axis lines on the top and right
  axis.ticks.length = unit(0.3, "cm"),
  panel.border = element_blank()  # No border around the plot panel
)

# Output to PDF
pdf(file = "Glu-complex.pdf", width = 10, height = 3)
complex
dev.off()
