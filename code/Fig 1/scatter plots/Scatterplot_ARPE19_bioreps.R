library(ggplot2)
library(dplyr)
library(patchwork)

### Read in count files
UNINF1_Pol3pri <- read.csv('./counts/ARPE19-UNINF-24h-4_bwa_W13k6T20.SR.primary-merged.HG38.sorted.PolIIItx.primaryOnly.count.txt', sep='\t', header = FALSE)
UNINF2_Pol3pri <- read.csv('./counts/ARPE19-UNINF-24h-6_bwa_W13k6T20.SR.primary-merged.HG38.sorted.PolIIItx.primaryOnly.count.txt', sep='\t', header = FALSE)
UNINF1_Pol3psu <- read.csv('./counts/ARPE19_UNINF_24h-4_bwa_W13k6T20.SR.primary-merged.HG38.sorted.PolIIItx.pseudoOnly.count.txt', sep='\t', header = FALSE)
UNINF2_Pol3psu <- read.csv('./counts/ARPE19_UNINF_24h-6_bwa_W13k6T20.SR.primary-merged.HG38.sorted.PolIIItx.pseudoOnly.count.txt', sep='\t', header = FALSE)
UNINF1_tRNA <- read.csv('./counts/ARPE19-UNINF-24h-4_bwa_W13k6T20.SR.primary-merged.HG38.sorted.hg38-tRNAs.count.txt', sep='\t', header = FALSE)
UNINF2_tRNA <- read.csv('./counts/ARPE19-UNINF-24h-6_bwa_W13k6T20.SR.primary-merged.HG38.sorted.hg38-tRNAs.count.txt', sep='\t', header = FALSE)
UNINF1_SNAR <- read.csv('./counts/ARPE19_UNINF_24h-4_bwa_W13k6T20.SR.primary-merged.HG38.sorted.SNARs.count.txt', sep='\t', header = FALSE)
UNINF2_SNAR <- read.csv('./counts/ARPE19_UNINF_24h-6_bwa_W13k6T20.SR.primary-merged.HG38.sorted.SNARs.count.txt', sep='\t', header = FALSE)
UNINF1_noPolIII <- read.csv('./counts/ARPE19_UNINF_24h-4_bwa_W13k6T20.SR.primary-merged.HG38.sorted.noPolIIItx.counts.txt', sep='\t', header = FALSE)
UNINF2_noPolIII <- read.csv('./counts/ARPE19_UNINF_24h-6_bwa_W13k6T20.SR.primary-merged.HG38.sorted.noPolIIItx.counts.txt', sep='\t', header = FALSE)

### Set columns names & generate sum + tpm counts
colnames(UNINF1_Pol3pri)<-c("count","ID")
colnames(UNINF2_Pol3pri)<-c("count","ID")
colnames(UNINF1_Pol3psu)<-c("count","ID")
colnames(UNINF2_Pol3psu)<-c("count","ID")
colnames(UNINF1_tRNA)<-c("count","ID")
colnames(UNINF2_tRNA)<-c("count","ID")
colnames(UNINF1_SNAR)<-c("count","ID")
colnames(UNINF2_SNAR)<-c("count","ID")
colnames(UNINF1_noPolIII)<-c("count","ID")
colnames(UNINF2_noPolIII)<-c("count","ID")

UNINF1_Pol3pri$category<-"primary"
UNINF2_Pol3pri$category<-"primary"
UNINF1_tRNA$category<-"pre-tRNA"
UNINF2_tRNA$category<-"pre-tRNA"
UNINF1_Pol3psu$category<-"pseudogene"
UNINF2_Pol3psu$category<-"pseudogene"
UNINF1_SNAR$category<-"SNAR"
UNINF2_SNAR$category<-"SNAR"
UNINF1_noPolIII$category<-"Pol II"
UNINF2_noPolIII$category<-"Pol II"

UNINF1_combined <- rbind(UNINF1_Pol3pri, UNINF1_tRNA, UNINF1_Pol3psu, UNINF1_SNAR, UNINF1_noPolIII)
UNINF2_combined <- rbind(UNINF2_Pol3pri, UNINF2_tRNA, UNINF2_Pol3psu, UNINF2_SNAR, UNINF2_noPolIII)

write.csv(UNINF1_combined, "Uninf-1.merged.counts.txt")
write.csv(UNINF2_combined, "Uninf-2.merged.counts.txt")

sum_UNINF1 <- sum(UNINF1_combined$count)
sum_UNINF2 <- sum(UNINF2_combined$count)

norm_UNINF1 <- UNINF1_combined %>%
  mutate(UNINF1_tpm = (UNINF1_combined$count / sum_UNINF1) * 1e6)

norm_UNINF2 <- UNINF2_combined %>%
  mutate(UNINF2_tpm = (UNINF2_combined$count / sum_UNINF2) * 1e6)

### PAO1 vs PA14wt
# merge dataframes into single dataframe
counts <- full_join(norm_UNINF1, norm_UNINF2, by = "ID") 

columns_to_drop <- c("category.y")
counts <- counts[, -which(names(counts) %in% columns_to_drop)]
names(counts)[names(counts) == "count.x"] <- "uninf1_count"
names(counts)[names(counts) == "count.y"] <- "uninf2_count"
names(counts)[names(counts) == "category.x"] <- "category"

### Reorder columns
col_order <- c("ID", "category", "uninf1_count", "uninf2_count", "UNINF1_tpm", "UNINF2_tpm")
counts_fixed <- counts[, col_order]

### Set NA to 0
counts_fixed <-replace(counts_fixed, is.na(counts_fixed), 0)

counts_filtered <- subset(counts_fixed, category != 0)

counts_filtered$category <- factor(counts_filtered$category, levels = c("Pol II", "SNAR", "pseudogene", "pre-tRNA", "primary"))  # Adjust the order as needed

counts_filtered <- counts_filtered[order(counts_filtered$category), ]

custom_colors <- c("Pol II" = "#00b050", "SNAR" = "#f3511f", "pseudogene" = "#c00000", "pre-tRNA" = "#002060", "primary" = "#f3511f")

# Calculate Spearman correlation coefficient
corr <- cor(counts_filtered$uninf1_count, counts_filtered$uninf2_count, method = "spearman")

corr2 <- corr^2



# Create a scatter plot
uninf_bioreps <- ggplot(counts_filtered, aes(x = UNINF1_tpm, y = UNINF2_tpm, color = category)) +
  geom_point(size = 1.5) +
  theme_classic() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
  labs(x = "ARPE-19 #1 (tpm)", 
       y = "ARPE-19 #2 (tpm)", 
       title = paste("Pearson r?=", round(corr2, 4))) + 
  scale_color_manual(values = custom_colors) +
  scale_x_log10() +
  scale_y_log10() +
  theme(plot.title = element_text(size = 20), 
        text = element_text(size = 22, colour = "black"),  # Set all text to black
        axis.line = element_line(size = 1.5, color = "black"),  # Thicken axis lines
        axis.text = element_text(size = 20, colour = "black"),  # Black axis text
        axis.ticks = element_line(size = 1, color = "black")) +  # Thicker axis ticks
        theme(legend.position = "none")

pdf(file = "arpe19_bioreps.pdf", width = 7, height = 7)
uninf_bioreps
dev.off()

