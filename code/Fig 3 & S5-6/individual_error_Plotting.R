# Load required libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

chromosome_filter <- "Homo_sapiens_tRNA-Glu-TTC-1"

# Function to read and preprocess data
process_variant_file <- function(file_path, chromosome_filter) {
  df <- read.csv(file_path, sep = '\t', header = TRUE)
  df <- df[!grepl("ins", df$Position), ]  # Remove insertions
  df <- df[df$Chromosome == chromosome_filter, ]  # Filter based on chromosome
  return(df)
}

# Process DRAP3R and nano data
DRAP3R <- process_variant_file('./variantFilter/ARPE19_UNINF_24h-4.sup-pseU_bwa_W13k6T20.SR.merged.primary.splint_tRNAdb.variantCalls.txt',chromosome_filter)
nano <- process_variant_file('./variantFilter/ARPE19_UNINF_tRNA-24h-4.sup-pseU_bwa_W13k6T20.SR.merged.primary.splint_tRNAdb.variantCalls.txt',chromosome_filter)

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

# Function to reshape data for plotting
reshape_data_for_plotting <- function(df) {
  reshaped <- data.frame(Position = integer(), Reference = character(), Count = integer(), stringsAsFactors = FALSE)
  for (i in 1:nrow(df)) {
    # Add the reference count for the "Ref" value (match to its column in the dataset)
    reshaped <- rbind(reshaped, data.frame(Position = df$Position[i], Reference = "Ref", Count = df[i, match(df$Reference[i], colnames(df))]))
    
    # Add counts for other bases (A, C, G, T, N) even if they are zero for the given reference
    for (base in c("A", "C", "G", "T", "N")) {
      reshaped <- rbind(reshaped, data.frame(Position = df$Position[i], Reference = base, Count = ifelse(base == df$Reference[i], 0, df[i, match(base, colnames(df))])))
    }
    
    # Handle "X" values properly
    if (!is.na(df$X[i]) && df$X[i] > 0) {
      reshaped <- rbind(reshaped, data.frame(Position = df$Position[i], Reference = "X", Count = df$X[i]))
    }
  }
  return(reshaped)
}



# Reshape both datasets
reshaped_DRAP3R <- reshape_data_for_plotting(DRAP3R)
reshaped_nano <- reshape_data_for_plotting(nano)

# Function to create fraction stacked bar plot
create_fraction_stacked_bar_plot <- function(df, df_name, x_range = NULL, position_to_letter_df = NULL) {
  df$Position <- as.numeric(df$Position)
  
  if (!is.null(x_range)) {
    df <- df %>% filter(Position >= x_range[1] & Position <= x_range[2])
  }
  
  df <- df %>%
    group_by(Position) %>%
    mutate(Fraction = Count / sum(Count))
  
  plot_title <- paste("Fraction Stacked Barplot for", df_name)
  fill_order <- c("N", "T", "G", "C", "A", "Ref", "X")
  
  p <- ggplot(df, aes(x = factor(Position, levels = as.character(unique(df$Position))),
                      y = Fraction, fill = factor(Reference, levels = fill_order),
                      color = factor(Reference, levels = fill_order))) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Ref" = "#d8d8d8", "A" = "#d66060", "C" = "#3eabca", "G" = "#d1cd5d", "T" = "#63ac52", "N" = "grey48", "X" = "white")) +
    scale_color_manual(values = c("Ref" = "black", "A" = "black", "C" = "black", "G" = "black", "T" = "black", "N" = "black", "X" = "black")) +
    labs(title = paste("Fraction Stacked Barplot for", df_name, chromosome_filter), x = "Position", y = "Fraction") +
    theme_minimal() +
    theme(legend.position = "right") +
    guides(fill = guide_legend(reverse = TRUE), color = FALSE) +
    guides(fill = guide_legend(title = NULL, reverse = TRUE))
  
  if (!is.null(position_to_letter_df) && all(c("Position", "Reference") %in% colnames(position_to_letter_df))) {
    p <- p + scale_x_discrete(breaks = position_to_letter_df$Position, labels = position_to_letter_df$Reference)
  }
  
  return(p)
}

# Prepare position-to-letter mapping and apply T -> U conversion
pos2letDRAP3R <- select(DRAP3R, c('Position', 'Reference'))
pos2letnano <- select(nano, c('Position', 'Reference'))

pos2letDRAP3R$Reference <- gsub("T", "U", pos2letDRAP3R$Reference)
pos2letnano$Reference <- gsub("T", "U", pos2letnano$Reference)

# Set x-axis range dynamically based on max position
Max <- max(as.numeric(pos2letDRAP3R$Position), na.rm = TRUE)
x_range_new <- c(25, Max - 33)

# Generate the plots
plot1 <- create_fraction_stacked_bar_plot(reshaped_DRAP3R, "DRAP3R", x_range = x_range_new, position_to_letter_df = pos2letDRAP3R)
plot2 <- create_fraction_stacked_bar_plot(reshaped_nano, "nano", x_range = x_range_new, position_to_letter_df = pos2letnano)

# Customize plot appearance
plot1 <- plot1 + theme(axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 16), legend.text = element_text(size = 16))
plot2 <- plot2 + theme(axis.text.y = element_text(size = 18), axis.title.y = element_text(size = 20), axis.title.x = element_text(size = 18), axis.text.x = element_text(size = 16), legend.text = element_text(size = 16))

# Combine the plots using patchwork
final_plot <- plot1 / plot2

# Generate the filename with the chromosome filter included
pdf_filename <- paste0("DRAP3R-nano_", chromosome_filter, "_fraction_stacked_bar.pdf")

# Save the final plot as a PDF
pdf(file = pdf_filename, width = 30, height = 10)
final_plot
dev.off()
