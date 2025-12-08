library(ggplot2)
library(dplyr)
library(patchwork)

setwd("")

### Read in count files
INF1_Pol3pri <- read.csv('ARPE19_HSV1_6hpi-1_bwa_W13k6T20.merged.HG38.primary.sorted.PolIIItx.primaryOnly.count.txt', sep='\t', header = FALSE)
INF2_Pol3pri <- read.csv('ARPE19_HSV1_6hpi-2_bwa_W13k6T20.merged.HG38.primary.sorted.PolIIItx.primaryOnly.count.txt', sep='\t', header = FALSE)
INF1_Pol3psu <- read.csv('ARPE19_HSV1_6hpi-1_bwa_W13k6T20.merged.HG38.primary.sorted.PolIIItx.pseudoOnly.count.txt', sep='\t', header = FALSE)
INF2_Pol3psu <- read.csv('ARPE19_HSV1_6hpi-2_bwa_W13k6T20.merged.HG38.primary.sorted.PolIIItx.pseudoOnly.count.txt', sep='\t', header = FALSE)
INF1_tRNA <- read.csv('ARPE19_HSV1_6hpi-1_bwa_W13k6T20.merged.HG38.primary.sorted.hg38-tRNAs.count.txt', sep='\t', header = FALSE)
INF2_tRNA <- read.csv('ARPE19_HSV1_6hpi-2_bwa_W13k6T20.merged.HG38.primary.sorted.hg38-tRNAs.count.txt', sep='\t', header = FALSE)
INF1_SNAR <- read.csv('ARPE19_HSV1_6hpi-1_bwa_W13k6T20.merged.HG38.primary.sorted.SNARs.count.txt', sep='\t', header = FALSE)
INF2_SNAR <- read.csv('ARPE19_HSV1_6hpi-2_bwa_W13k6T20.merged.HG38.primary.sorted.SNARs.count.txt', sep='\t', header = FALSE)
INF1_noPolIII <- read.csv('ARPE19_HSV1_6hpi-1_bwa_W13k6T20.merged.HG38.primary.sorted.noPolIIItx.counts.txt', sep='\t', header = FALSE)
INF2_noPolIII <- read.csv('ARPE19_HSV1_6hpi-2_bwa_W13k6T20.merged.HG38.primary.sorted.noPolIIItx.counts.txt', sep='\t', header = FALSE)

### Set columns names & generate sum + tpm counts
colnames(INF1_Pol3pri)<-c("count","ID")
colnames(INF2_Pol3pri)<-c("count","ID")
colnames(INF1_Pol3psu)<-c("count","ID")
colnames(INF2_Pol3psu)<-c("count","ID")
colnames(INF1_tRNA)<-c("count","ID")
colnames(INF2_tRNA)<-c("count","ID")
colnames(INF1_SNAR)<-c("count","ID")
colnames(INF2_SNAR)<-c("count","ID")
colnames(INF1_noPolIII)<-c("count","ID")
colnames(INF2_noPolIII)<-c("count","ID")

INF1_Pol3pri$category<-"primary"
INF2_Pol3pri$category<-"primary"
INF1_tRNA$category<-"pre-tRNA"
INF2_tRNA$category<-"pre-tRNA"
INF1_Pol3psu$category<-"pseudogene"
INF2_Pol3psu$category<-"pseudogene"
INF1_SNAR$category<-"SNAR"
INF2_SNAR$category<-"SNAR"
INF1_noPolIII$category<-"Pol II"
INF2_noPolIII$category<-"Pol II"

INF1_combined <- rbind(INF1_Pol3pri, INF1_tRNA, INF1_Pol3psu, INF1_SNAR, INF1_noPolIII)
INF2_combined <- rbind(INF2_Pol3pri, INF2_tRNA, INF2_Pol3psu, INF2_SNAR, INF2_noPolIII)

sum_INF1 <- sum(INF1_combined$count)
sum_INF2 <- sum(INF2_combined$count)

norm_INF1 <- INF1_combined %>%
  mutate(INF1_tpm = (INF1_combined$count / sum_INF1) * 1e6)

norm_INF2 <- INF2_combined %>%
  mutate(INF2_tpm = (INF2_combined$count / sum_INF2) * 1e6)

### PAO1 vs PA14wt
# merge dataframes into single dataframe
counts <- full_join(norm_INF1, norm_INF2, by = "ID") 

columns_to_drop <- c("category.y")
counts <- counts[, -which(names(counts) %in% columns_to_drop)]
names(counts)[names(counts) == "count.x"] <- "inf1_count"
names(counts)[names(counts) == "count.y"] <- "inf2_count"
names(counts)[names(counts) == "category.x"] <- "category"

### Reorder columns
col_order <- c("ID", "category", "inf1_count", "inf2_count", "INF1_tpm", "INF2_tpm")
counts_fixed <- counts[, col_order]

### Set NA to 0
counts_fixed <-replace(counts_fixed, is.na(counts_fixed), 0)

counts_filtered <- subset(counts_fixed, category != 0)

counts_filtered$category <- factor(counts_filtered$category, levels = c("Pol II", "SNAR", "pseudogene", "pre-tRNA", "primary"))  # Adjust the order as needed

counts_filtered <- counts_filtered[order(counts_filtered$category), ]

custom_colors <- c("Pol II" = "#00b050", "SNAR" = "#f3511f", "pseudogene" = "#c00000", "pre-tRNA" = "#002060", "primary" = "#f3511f")

# Calculate Pearson correlation coefficient
corr <- cor(counts_filtered$inf1_count, counts_filtered$inf2_count, method = "pearson")

corr2 <- corr^2


# Create a scatter plot
inf_bioreps <- ggplot(counts_filtered, aes(x = INF1_tpm, y = INF2_tpm, color = category)) +
  geom_point(size = 1.5) +
  theme_classic() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "HSV-1 6hpi #1 (tpm)", 
       y = "HSV-1 6hpi #2 (tpm)", 
       title = paste("Pearson r2=", round(corr2, 4))) + 
  scale_color_manual(values = custom_colors) +
  scale_x_log10() +
  scale_y_log10() +
  theme(plot.title = element_text(size = 20), 
        text = element_text(size = 22, colour = "black"),  # Set all text to black
        axis.line = element_line(size = 1.5, color = "black"),  # Thicken axis lines
        axis.text = element_text(size = 20, colour = "black"),  # Black axis text
        axis.ticks = element_line(size = 1, color = "black")) +  # Thicker axis ticks
  theme(legend.position = "none")


pdf(file = "inf_6hpi_bioreps.pdf", width = 6, height = 6)
inf_bioreps
dev.off()

