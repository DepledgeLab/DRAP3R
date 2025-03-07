library(data.table)
library(dplyr)
library(tidyr)
library (ggplot2)
library(scales)
library(patchwork)

setwd("C:/Users/depledgd/Dropbox/DRAP3R/Fig 2 - mods-new/")
#setwd("D:/Dropbox/DRAP3R/Fig 2 - mods-new/")

### Set basic filtering parameters
COV <- 20
FREQ <- 5

### Read in  transcriptome summary
base <- fread("hg38-output.txt", header=TRUE)
#base$gene_name <- sub("^mRNA\\.", "", base$gene_name)
base$result <- NA


###### SECTION FOR PROCESSING m6A DATA ######

### Read in Dorado m6A output
m6A_sites <- fread("ARPE19_polyA_UNINF_24h-6.sup-m6A.hg38tx.modFilt.0.98.cov20filt.sorted.bed", header = FALSE)

colnames(m6A_sites) <- c("chrom", "start","end", "mod","score", "strand","start_tmp", "end_tmp","color", "valid_cov", "frac_mod", "N_mod","N_canon", "N_other", "N_del","N_fail", "N_diff","N_nocall")

m6A_sites$chrom <- sub("\\|.*", "", m6A_sites$chrom)

m6A_split <- m6A_sites %>%
  filter(valid_cov >= COV, frac_mod > FREQ, N_diff <= valid_cov) #N_fail <= valid_cov/2

### Perform the join and calculation using dplyr
m6A_split <- m6A_split %>%
  left_join(base, by = c("chrom" = "gene_name")) %>%
  rowwise() %>%
  mutate(result = ifelse(is.na(start), NA, 
                         ifelse(start < utr5_size, (start / utr5_size), 
                                ifelse(start < (utr5_size + cds_size), ((start - utr5_size) / cds_size ) + 1, 
                                       ((start - utr5_size - cds_size) / utr3_size) + 2)))) %>%
  ungroup()  # Ensure we ungroup after row-wise operations

### Print the resulting dataframe m6A_split
#print(m6A_split)

### Fix dataframe and perform filtering based on coverage and frac_mod
m6A_split$frac_mod <- as.numeric(m6A_split$frac_mod)
m6A_split$valid_cov <- as.numeric(m6A_split$valid_cov)
m6A_parsed <- m6A_split %>% filter(valid_cov >= COV & frac_mod >= FREQ)


### rescale 5'UTR and 3'UTR
utr5.SF <- median(m6A_parsed$utr5_size, na.rm = T)/median(m6A_parsed$cds_size, na.rm = T)
utr3.SF <- median(m6A_parsed$utr3_size, na.rm = T)/median(m6A_parsed$cds_size, na.rm = T)

utr5.m6a.dist <- m6A_parsed[m6A_parsed$result < 1, ]
cds.m6a.dist <- m6A_parsed[m6A_parsed$result < 2 & m6A_parsed$result >= 1, ]
utr3.m6a.dist <- m6A_parsed[m6A_parsed$result >= 2, ]

utr5.m6a.dist$result <- rescale(utr5.m6a.dist$result, to = c(1-utr5.SF, 1), from = c(0,1))
utr3.m6a.dist$result <- rescale(utr3.m6a.dist$result, to = c(2, 2+utr3.SF), from = c(2,3))
m6a.metagene.coord <- c(utr5.m6a.dist$result, cds.m6a.dist$result, utr3.m6a.dist$result)

df <- data.frame(m6a.metagene.coord)

m6A<-ggplot(df) + geom_density(aes(x = m6a.metagene.coord), size = 1.5,color="red") + ylim(0, 1.5) + xlim(0,3) +
  geom_vline(xintercept = 1:2, col = "grey") + 
  theme_classic() + 
  theme(
    legend.position = "none", 
    axis.text = element_text(colour = "black", size = 24), 
    plot.title = element_text(size = 24, hjust = 0.5), 
    text = element_text(size = 24, color = "black"), 
    axis.line = element_line(size = 1), 
    axis.title = element_text(size = 24)
  )



###### SECTION FOR PROCESSING m6A DATA ######

### Read in Dorado m6A output
pseU_sites <- fread("ARPE19_polyA_UNINF_24h-6.sup-pseU.hg38tx.modFilt.0.98.cov20filt.sorted.bed", header = FALSE)

colnames(pseU_sites) <- c("chrom", "start","end", "mod","score", "strand","start_tmp", "end_tmp","color", "valid_cov", "frac_mod", "N_mod","N_canon", "N_other", "N_del","N_fail", "N_diff","N_nocall")

pseU_sites$chrom <- sub("\\|.*", "", pseU_sites$chrom)

pseU_split <- pseU_sites %>%
  filter(valid_cov >= COV, frac_mod > FREQ, N_diff <= valid_cov)




### Perform the join and calculation using dplyr
pseU_split <- pseU_split %>%
  left_join(base, by = c("chrom" = "gene_name")) %>%
  rowwise() %>%
  mutate(result = ifelse(is.na(start), NA, 
                         ifelse(start < utr5_size, (start / utr5_size), 
                                ifelse(start < (utr5_size + cds_size), ((start - utr5_size) / cds_size ) + 1, 
                                       ((start - utr5_size - cds_size) / utr3_size) + 2)))) %>%
  ungroup()  # Ensure we ungroup after row-wise operations

### Fix dataframe and perform filtering based on coverage and frac_mod
pseU_split$frac_mod <- as.numeric(pseU_split$frac_mod)
pseU_split$valid_cov <- as.numeric(pseU_split$valid_cov)
pseU_parsed <- pseU_split %>% filter(valid_cov >= COV & frac_mod >= FREQ)

### rescale 5'UTR and 3'UTR
utr5.SF <- median(pseU_parsed$utr5_size, na.rm = T)/median(pseU_parsed$cds_size, na.rm = T)
utr3.SF <- median(pseU_parsed$utr3_size, na.rm = T)/median(pseU_parsed$cds_size, na.rm = T)

utr5.pseU.dist <- pseU_parsed[pseU_parsed$result < 1, ]
cds.pseU.dist <- pseU_parsed[pseU_parsed$result < 2 & pseU_parsed$result >= 1, ]
utr3.pseU.dist <- pseU_parsed[pseU_parsed$result >= 2, ]

utr5.pseU.dist$result <- rescale(utr5.pseU.dist$result, to = c(1-utr5.SF, 1), from = c(0,1))
utr3.pseU.dist$result <- rescale(utr3.pseU.dist$result, to = c(2, 2+utr3.SF), from = c(2,3))
pseU.metagene.coord <- c(utr5.pseU.dist$result, cds.pseU.dist$result, utr3.pseU.dist$result)

df2 <- data.frame(pseU.metagene.coord)

pseU<-ggplot(df2) + geom_density(aes(x = pseU.metagene.coord),size = 1.5,color="blue") + ylim(0,1.5) + xlim(0,3) +
  geom_vline(xintercept = 1:2, col = "grey") + 
  theme_classic() + 
  theme(
    legend.position = "none", 
    axis.text = element_text(colour = "black", size = 24), 
    plot.title = element_text(size = 24, hjust = 0.5), 
    text = element_text(size = 24, color = "black"), 
    axis.line = element_line(size = 1), 
    axis.title = element_text(size = 24)
  )



### FINAL PLOTS + EXPORT


final <- (m6A | pseU)

pdf(file = "metaplot-polyA.pdf", width = 12, height = 6) 
final
dev.off()