library(ggplot2)
library(tidyr)
library(dplyr)

data <- read.csv("ARPE19_UNINF_24h-4_bwa_W13k6T20.merged.primary.sorted.HG38.bed", sep="\t", header=FALSE)
data<-data[, c("V2","V3","V4")]
data$aln_length <- data$V3 - data$V2
colnames(data) <- c("start","stop","read_id","aln_length")


SeqSum <- read.table("ARPE19_UNINF_24h-4.hac.notrim.seqSum.txt", sep="\t", header=F, stringsAsFactors=FALSE)
colnames(SeqSum) <- c("read_id","sequence_length")

merged <- merge(data, SeqSum, by = "read_id")

merged$sequence_length <- as.numeric(merged$sequence_length)

# Set a seed for reproducibility (optional)
set.seed(123)

# Subsample 500,000 rows
sample_size <- 500000
if (nrow(merged) < sample_size) {
  stop("The dataframe does not have enough rows to sample.")
}
sample_indices <- sample(seq_len(nrow(merged)), size = sample_size)

# Create the new dataframe with the sampled rows
merged_sampled <- merged[sample_indices, ]


### MAKE SCATTER PLOT
p<-ggplot(merged_sampled, aes(sequence_length, aln_length)) +
  geom_point(size = 0.05) + theme_classic() + ylim(0,1200) + xlim(0,1200) + 
  geom_abline(intercept = 0, slope = 1, color = "blue", linetype = "dashed") +
  geom_abline(intercept = -70, slope = 1, color = "red", linetype = "dashed") +
  labs(title = "read length vs alignment length", x = "read length", y = "alignment length") +
  theme(plot.title = element_text(size = 20), 
        text = element_text(size = 18, colour = "black"),  # Set all text to black
        axis.line = element_line(size = 1, color = "black"),  # Thicken axis lines
        axis.text = element_text(size = 18, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, colour = "black"),# Black axis text
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks = element_line(size = 1, color = "black"))  # Thicker axis ticks


pdf(file = "readLengthVSalnLength.pdf", width = 7, height = 6)
p
dev.off()



