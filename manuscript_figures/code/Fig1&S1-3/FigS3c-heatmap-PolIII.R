library(dplyr)
library(tidyr)

setwd("")

# Set path to your folder of files
files <- list.files(path = "./", pattern = ".*PolIII.*\\.txt", full.names = TRUE)

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

library(biomaRt)

# Use Ensembl BioMart
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Clean transcript IDs (remove .1, .2, etc.)
rownames(tpm_df) <- tpm_df$transcript_id
tpm_df$transcript_id <- NULL
transcripts <- sub("\\..*", "", rownames(tpm_df))

# Query gene names
transcript2gene <- getBM(
  attributes = c("ensembl_transcript_id", "external_gene_name"),
  filters = "ensembl_transcript_id",
  values = transcripts,
  mart = mart
)

# Map gene names to transcript IDs
mapping_df <- data.frame(
  ensembl_transcript_id = transcripts,
  rownames = rownames(tpm_df)
)
mapping_df <- left_join(mapping_df, transcript2gene, by = "ensembl_transcript_id")

# Fallback if gene name is missing
gene_names <- ifelse(
  is.na(mapping_df$external_gene_name) | mapping_df$external_gene_name == "",
  mapping_df$ensembl_transcript_id,
  mapping_df$external_gene_name
)

# Replace row names in TPM matrix
rownames(tpm_df) <- make.unique(gene_names)
#rownames(tpm_df) <- gene_names


# Add gene names as column
tpm_df$gene_name <- gene_names

# Collapse by gene (sum TPM across transcripts)
collapsed_df <- tpm_df %>%
  group_by(gene_name) %>%
  summarise(across(everything(), sum))

# Clean up
collapsed_df <- as.data.frame(collapsed_df)
rownames(collapsed_df) <- collapsed_df$gene_name
collapsed_df$gene_name <- NULL

tpm_df2<-collapsed_df

library(pheatmap)

# Set transcript_id as rownames
#rownames(tpm_df) <- tpm_df$transcript_id
tpm_df2 <- tpm_df2[, count_cols]

# Optional: log-transform for better visualization
tpm_log <- log2(tpm_df2 + 1)

# Optional: row-scaling
#tpm_scaled <- t(scale(t(tpm_log)))  # scale across samples

rows_to_remove <- c("5S_rRNA", "7SK")
tpm_log <- tpm_log[!(rownames(tpm_log) %in% rows_to_remove), ]

# Plot heatmap
plot<-pheatmap(tpm_log,
               # clustering_distance_rows = "euclidean",
               clustering_distance_cols = "euclidean",
               # clustering_method = "complete",
               cluster_rows = FALSE,
               #  cluster_cols = FALSE,
               show_rownames = TRUE,
               fontsize_row = 12,
               fontsize_col = 12)

pdf(file = "PolIIItx-heatmap.pdf", width = 5, height = 10)
plot
dev.off()
