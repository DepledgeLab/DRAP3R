library(ggplot2)
library(dplyr)
library(DESeq2)
library(EnhancedVolcano)

setwd("")

### Read in count files
UNINF_1 <- read.csv('./Uninf-1.merged.counts.txt', sep=',', header = TRUE, row.names = 1)
UNINF_2 <- read.csv('./Uninf-2.merged.counts.txt', sep=',', header = TRUE, row.names = 1)
INF_1 <- read.csv('./HSV1-12hpi-1.merged.counts.txt', sep=',', header = TRUE, row.names = 1)
INF_2 <- read.csv('./HSV1-12hpi-2.merged.counts.txt', sep=',', header = TRUE, row.names = NULL)

### drop first and last column
UNINF_1 <- UNINF_1[, -ncol(UNINF_1)]
UNINF_2 <- UNINF_2[, -ncol(UNINF_2)]
INF_1 <- INF_1[, -ncol(INF_1)]
INF_2 <- INF_2[, -ncol(INF_2)]

### Set columns names
colnames(UNINF_1)<-c("UNINF_1","ID")
colnames(UNINF_2)<-c("UNINF_2","ID")
colnames(INF_1)<-c("INF_1","ID")
colnames(INF_2)<-c("INF_2","ID")

### merge dataframes into single dataframe & reorder columns
counts <- full_join(UNINF_1, UNINF_2, by = "ID") %>%
  full_join(INF_1, by = "ID") %>%
  full_join(INF_2, by = "ID")
col_order <- c("ID", "UNINF_1", "UNINF_2", "INF_1", "INF_2")
counts_fixed <- counts[, col_order]

### Drop rows with NA
counts_fixed <- na.omit(counts_fixed)

countTable <- counts_fixed
rownames(countTable) <- countTable$ID
countTable <- countTable[, -1]

### Read in experimental setup
samples <- read.table("samples.txt", header=TRUE)

### Make DeSeq dataset
Dataset <- DESeqDataSetFromMatrix(countData = countTable, colData=samples, design=~condition)

# Run DESEQ and generate a simple plot showing the distribution of regulated and unregulated genes
Dataset$condition <- relevel(Dataset$condition, "UNINF")
DatasetProcessed <- DESeq(Dataset)

### Perform contrast analyses to produce lists of differentially regulated genes between conditions (pairwise)
result <- lfcShrink(DatasetProcessed, contrast=c("condition","INF","UNINF"), type="ashr")
baseMean_INF = rowMeans(counts(DatasetProcessed,normalized=TRUE)[,DatasetProcessed$condition == "INF"])
baseMean_UNINF = rowMeans(counts(DatasetProcessed,normalized=TRUE)[,DatasetProcessed$condition == "UNINF"])
result = cbind(as.data.frame(result), baseMean_INF, baseMean_UNINF)
write.csv(result, "deseq2-INF12h-vs-UNINF.output.csv", row.names=TRUE)

# Read the BED file (tab-separated, no header)
pol3_bed <- read.delim("gencode.v45.PolIIItx.primaryOnly.bed", header = FALSE, stringsAsFactors = FALSE)
pol3_ensts <- pol3_bed$V4
pol3_ensts <- intersect(rownames(result), pol3_ensts)

# Step 1: Create base color vector (all black)
# Step 1: Initialize color vector (black by default)
custom_colors <- setNames(rep("lightgrey", nrow(result)), rownames(result))

# Step 2: Define significant tRNAs (padj < 0.05, abs(log2FC) > 1, and name starts with tRNA)
sig_tRNAs <- rownames(result)[
  result$padj < 0.05 &
    abs(result$log2FoldChange) > 1 &
    grepl("^tRNA", rownames(result))
]

# Safety check
sig_tRNAs <- sig_tRNAs[sig_tRNAs %in% rownames(result)]

# Step 3: Define all significant transcripts not starting with tRNA
sig_tx <- rownames(result)[
  result$padj < 0.05 &
    abs(result$log2FoldChange) > 1 &
    !grepl("^tRNA", rownames(result))
]

sig_P3tx <- sig_tx[sig_tx %in% pol3_ensts]    # Pol III
sig_P2tx <- sig_tx[!sig_tx %in% pol3_ensts]   # Not Pol III = assumed Pol II

sig_tRNAs <- intersect(sig_tRNAs, rownames(result))
custom_colors[sig_tRNAs] <- "#002060"

sig_P3tx <- intersect(sig_P3tx, rownames(result))
custom_colors[sig_P3tx] <- "#f3511f"

sig_P2tx <- intersect(sig_P2tx, rownames(result))
custom_colors[sig_P2tx] <- "#00b050"

volc <- EnhancedVolcano(result,
                        x = 'log2FoldChange',
                        y = 'padj',
                        lab = NA,
                        pointSize = 1,
                        FCcutoff = 1,
                        pCutoff = 0.05,
                        colCustom = custom_colors,
                        colAlpha = 1,
                        gridlines.major = TRUE,
                        gridlines.minor = FALSE,
                        title = 'INF vs. UNINF',
                        hline = 0.05,
                        xlim = c(-2.5, 2.5),
                        ylim = c(0, 10)) +
  theme(axis.line = element_line(size = 1.5),
        legend.position = "none") +  # 👈 hides the legend
  xlab("Log2 Fold Change") +
  ylab("-Log10 (adjusted p)")

volc

pdf(file = "volcano-INF12h-vs-UNINF.pdf", width = 7, height = 9)
print(volc)
dev.off()

