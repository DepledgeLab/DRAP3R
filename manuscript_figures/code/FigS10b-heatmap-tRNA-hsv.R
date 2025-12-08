library(dplyr)
library(tidyr)

setwd("")

# Set path to your folder of files
files <- list.files(path = "./", pattern = ".*hg38.*\\.txt", full.names = TRUE)

# Now exclude unwanted patterns
files <- files[!grepl("PolIIItx\\.pseudo|noPolIIItx|UNINF", files)]

# Function to read a file and name columns based on filename
read_file <- function(file) {
  # Get filename without path and extension
  name_raw <- tools::file_path_sans_ext(basename(file))
  
  # Clean it up: remove suffix like ".hg38_counts"
  sample_name <- gsub("_bwa.*", "", name_raw)
  
  # Read and rename count column
  df <- read.table(file, header = FALSE, col.names = c("count", "transcript_id"))
  df <- df %>% rename(!!sample_name := count)
}

# Read and merge all files
data_list <- lapply(files, read_file)
merged_df <- Reduce(function(x, y) full_join(x, y, by = "transcript_id"), data_list)

# Replace NAs with 0 (assuming missing values mean 0 count)
merged_df[is.na(merged_df)] <- 0

# Select count columns only (excluding 'transcript_id')
count_cols <- setdiff(names(merged_df), "transcript_id")

# Compute TPM for each sample
tpm_df <- merged_df
tpm_df[count_cols] <- lapply(merged_df[count_cols], function(x) {
  x / sum(x, na.rm = TRUE) * 1e6
})


library(pheatmap)

# Set transcript_id as rownames
rownames(tpm_df) <- tpm_df$transcript_id
tpm_df <- tpm_df[, count_cols]

# Optional: log-transform for better visualization
tpm_log <- log2(tpm_df + 1)

# Optional: row-scaling
tpm_scaled <- t(scale(t(tpm_log)))  # scale across samples

# Plot heatmap
plot<-pheatmap(tpm_log,
               clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               show_rownames = FALSE,
               fontsize_row = 12,
               fontsize_col = 12)


pdf(file = "tRNA-HSV-heatmap.pdf", width = 5, height = 10)
plot
dev.off()


