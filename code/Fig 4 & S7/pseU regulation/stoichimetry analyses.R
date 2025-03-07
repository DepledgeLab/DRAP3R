library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

### Set basic filtering parameters

COV <- 20
FREQ <- 1

# Step 1: Group file names based on their types
### tRNAs
file_groups <- list(
  pseU_sites = list(
    files = list(
      "ARPE19_HSV1_6hpi-2.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed",
      "ARPE19_HSV1_6hpi-3.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed",
      "ARPE19-HSV1-12hpi-1.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed",
      "../../Fig 2/stochiometry histograms/IntersectMods-pseU/ARPE19_UNINF_24h-4.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed",
      "../../Fig 2/stochiometry histograms/IntersectMods-pseU/ARPE19_UNINF_24h-6.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed",
      "ARPE19-IFNag-24h-1.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed"
    ),
    variable_names = list(
      "HSV1_6hpi_1_tRNA_pseU_sites",
      "HSV1_6hpi_2_tRNA_pseU_sites",
      "HSV1_12hpi_1_tRNA_pseU_sites",
      "UNINF_1_tRNA_pseU_sites",
      "UNINF_2_tRNA_pseU_sites",
      "IFN_tRNA_pseU_sites"
    )
  ),
 PolIIItx_pseU_sites = list(
    files = list(
      "ARPE19_HSV1_6hpi-2.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "ARPE19_HSV1_6hpi-3.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "ARPE19-HSV1-12hpi-1.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "../../Fig 2/stochiometry histograms/IntersectMods-pseU/ARPE19_UNINF_24h-4.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "../../Fig 2/stochiometry histograms/IntersectMods-pseU/ARPE19_UNINF_24h-6.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "ARPE19-IFNag-24h-1.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed"
    ),
    variable_names = list(
      "HSV1_6hpi_1_PolIIItx_pseU_sites",
      "HSV1_6hpi_2_PolIIItx_pseU_sites",
      "HSV1_12hpi_1_PolIIItx_pseU_sites",
      "UNINF_1_PolIIItx_pseU_sites",
      "UNINF_2_PolIIItx_pseU_sites",
      "IFN_PolIIItx_pseU_sites"
    )
  )
)

# Define column names
column_names <- list(
  "chrom", "start", "end", "mod", "score", "strand", "start_tmp", "end_tmp", 
  "color", "valid_cov", "frac_mod", "N_mod", "N_canon", "N_other", "N_del", 
  "N_fail", "N_diff", "N_nocall", "chrom2", "start2", "stop2", "tRNA", 
  "na1", "na2", "na3", "na4", "na5", "exons", "blockSizes", "blockStarts", "na6"
)

# Initialize an empty list to store the resulting dataframes
dataframes <- list()

# Process each group of files
for (group in names(file_groups)) {
  files <- file_groups[[group]]$files
  variable_names <- file_groups[[group]]$variable_names
  
  for (i in seq_along(files)) {
    # Read the data from the file with mixed delimiters using read.table
    mod_sites <- read.table(files[[i]], header = FALSE, sep = "", fill = TRUE)
    
    # Set column names
    colnames(mod_sites) <- unlist(column_names)
    
    # Filter rows where valid_cov is below the specified threshold
    mod_sites <- mod_sites %>%
      filter(valid_cov >= COV, frac_mod > FREQ, N_diff <= valid_cov, N_fail <= valid_cov/2)
    
    # Ensure result is initialized
    mod_sites <- mod_sites %>% mutate(result = NA)
    
    # Apply the transformation conditionally
    mod_sites <- mod_sites %>%
      mutate(result = if_else(
        strand == "+",
        (start - start2) / (stop2 - start2),
        if_else(strand == "-",
                (stop2 - end) / (stop2 - start2),
                NA_real_)
      ))
    
    # Store the dataframe in the list with the corresponding variable name
    dataframes[[variable_names[[i]]]] <- mod_sites
    
    # Optionally, assign the dataframe to a variable with the constructed name
    assign(variable_names[[i]], mod_sites)
  }
}

### INIDIVIDUAL PLOTS

h6_1_tRNA<-ggplot() +
  geom_histogram(data = HSV1_6hpi_1_tRNA_pseU_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "#5979bc", alpha = 0.9, color = "black", boundary = 0) +  # First histogram
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Modification stoichiometry", y = "Count") + #, title = "ARPE19 rep2 pre-tRNA") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 24),  # Axis text size
    axis.title = element_text(size = 24)  # Axis title size
  ) +
  ylim(0, 500)

h6_1_pIII<-ggplot() +
  geom_histogram(data = HSV1_6hpi_1_PolIIItx_pseU_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "#5979bc", alpha = 0.9, color = "black", boundary = 0) +  # First histogram
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Modification stoichiometry", y = "Count") + #, title = "ARPE19 rep1 Pol III tx") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 24),  # Axis text size
    axis.title = element_text(size = 24)  # Axis title size
  ) +
  ylim(0, 50)


h6_2_tRNA<-ggplot() +
  geom_histogram(data = HSV1_6hpi_2_tRNA_pseU_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "#5979bc", alpha = 0.9, color = "black", boundary = 0) +  # First histogram
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Modification stoichiometry", y = "Count") + #, title = "ARPE19 rep2 pre-tRNA") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 24),  # Axis text size
    axis.title = element_text(size = 24)  # Axis title size
  ) +
  ylim(0, 500)

h6_2_pIII<-ggplot() +
  geom_histogram(data = HSV1_6hpi_2_PolIIItx_pseU_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "#5979bc", alpha = 0.9, color = "black", boundary = 0) +  # First histogram
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Modification stoichiometry", y = "Count") + #, title = "ARPE19 rep1 Pol III tx") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 24),  # Axis text size
    axis.title = element_text(size = 24)  # Axis title size
  ) +
  ylim(0, 50)


h12_1_tRNA<-ggplot() +
  geom_histogram(data = HSV1_12hpi_1_tRNA_pseU_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "#5979bc", alpha = 0.9, color = "black", boundary = 0) +  # First histogram
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Modification stoichiometry", y = "Count") + #, title = "ARPE19 rep2 pre-tRNA") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 24),  # Axis text size
    axis.title = element_text(size = 24)  # Axis title size
  ) +
 ylim(0, 500)

h12_1_pIII<-ggplot() +
  geom_histogram(data = HSV1_12hpi_1_PolIIItx_pseU_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "#5979bc", alpha = 0.9, color = "black", boundary = 0) +  # First histogram
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Modification stoichiometry", y = "Count") + #, title = "ARPE19 rep1 Pol III tx") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 24),  # Axis text size
    axis.title = element_text(size = 24)  # Axis title size
  ) +
  ylim(0, 50)


IFN_1_tRNA<-ggplot() +
  geom_histogram(data = IFN_tRNA_pseU_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "#5979bc", alpha = 0.9, color = "black", boundary = 0) +  # First histogram
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Modification stoichiometry", y = "Count") + #, title = "ARPE19 rep2 pre-tRNA") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 24),  # Axis text size
    axis.title = element_text(size = 24)  # Axis title size
  ) +
  ylim(0, 500)

IFN_1_pIII<-ggplot() +
  geom_histogram(data = IFN_PolIIItx_pseU_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "#5979bc", alpha = 0.9, color = "black", boundary = 0) +  # First histogram
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Modification stoichiometry", y = "Count") + #, title = "ARPE19 rep1 Pol III tx") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 24),  # Axis text size
    axis.title = element_text(size = 24)  # Axis title size
  ) +
  ylim(0, 50)






combined_plot2 <- h6_1_tRNA + h6_2_tRNA + h12_1_tRNA + IFN_1_tRNA + h6_2_pIII + h6_2_pIII + h12_1_pIII + IFN_1_pIII + plot_layout(ncol = 4)


pdf(file = "HSV1 FreqHistograms.pdf", width = 15, height = 10)
combined_plot2
dev.off()




### THE FUN STUFF PART 1 - tRNA stoichiometry analyses

### comparison of stoichiometries in uninfected
uninf_comp <- merge(UNINF_1_tRNA_pseU_sites, UNINF_2_tRNA_pseU_sites, by = c("start", "end", "tRNA"), all = TRUE)
uninf_comp <- na.omit(uninf_comp)
uninf_comp <- uninf_comp[, c("start", "end", "tRNA", "frac_mod.x", "frac_mod.y")]
uninf_R <- cor(uninf_comp$frac_mod.x, uninf_comp$frac_mod.y, method="spearman")
uninf_R2 <- uninf_R^2

un<-ggplot(uninf_comp, aes(x = frac_mod.x, y = frac_mod.y)) +
  geom_point(color = "#5979bc", alpha = 0.5) +  # Scatter plot with blue points
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +  # x = y line (y = x)
  labs(x = "Uninfected ARPE19 #1", y = "Uninfected ARPE19 #2", 
    title = paste("Pearson r2=", round(uninf_R, 4))) + 
  theme_classic() +
  theme(plot.title = element_text(size = 20), 
        text = element_text(size = 22, colour = "black"),  # Set all text to black
        axis.line = element_line(size = 1.5, color = "black"),  # Thicken axis lines
        axis.text = element_text(size = 20, colour = "black"))  # Black axis text

### comparison of stoichiometries in infected
inf_comp <- merge(HSV1_6hpi_1_tRNA_pseU_sites, HSV1_6hpi_2_tRNA_pseU_sites, by = c("start", "end", "tRNA"), all = TRUE)
inf_comp <- na.omit(inf_comp)
inf_comp <- inf_comp[, c("start", "end", "tRNA", "frac_mod.x", "frac_mod.y")]
inf_R <- cor(inf_comp$frac_mod.x, inf_comp$frac_mod.y, method="spearman")
inf_R2 <- inf_R^2

inf<-ggplot(inf_comp, aes(x = frac_mod.x, y = frac_mod.y)) +
  geom_point(color = "#5979bc", alpha = 0.5) +  # Scatter plot with blue points
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +  # x = y line (y = x)
  labs(x = "HSV-1 6 hpi #1", y = "HSV-1 6 hpi #2", 
       title = paste("Pearson r2=", round(inf_R, 4))) + 
  theme_classic() +
  theme(plot.title = element_text(size = 20), 
        text = element_text(size = 22, colour = "black"),  # Set all text to black
        axis.line = element_line(size = 1.5, color = "black"),  # Thicken axis lines
        axis.text = element_text(size = 20, colour = "black"))  # Black axis text


### comparison of stoichiometries: uninf vs 6hpi
comp6 <- merge(UNINF_1_tRNA_pseU_sites, HSV1_6hpi_1_tRNA_pseU_sites, by = c("start", "end", "tRNA"), all = TRUE)
comp6<- na.omit(comp6)
comp6 <- comp6[, c("start", "end", "tRNA", "frac_mod.x", "frac_mod.y")]
comp6_R <- cor(comp6$frac_mod.x, comp6$frac_mod.y, method="spearman")
comp6_R2 <- comp6_R^2

c6<-ggplot(comp6, aes(x = frac_mod.x, y = frac_mod.y)) +
  geom_point(color = "#5979bc", alpha = 0.5) +  # Scatter plot with blue points
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +  # x = y line (y = x)
  labs(x = "Uninfected ARPE19 #1", y = "HSV-1 6 hpi #1", 
       title = paste("Pearson r2=", round(comp6_R, 4))) + 
  theme_classic() +
  theme(plot.title = element_text(size = 20), 
        text = element_text(size = 22, colour = "black"),  # Set all text to black
        axis.line = element_line(size = 1.5, color = "black"),  # Thicken axis lines
        axis.text = element_text(size = 20, colour = "black"))  # Black axis text

### comparison of stoichiometries: uninf vs 6hpi
comp12 <- merge(UNINF_1_tRNA_pseU_sites, HSV1_12hpi_1_tRNA_pseU_sites, by = c("start", "end", "tRNA"), all = TRUE)
comp12<- na.omit(comp12)
comp12 <- comp12[, c("start", "end", "tRNA", "frac_mod.x", "frac_mod.y")]
comp12_R <- cor(comp12$frac_mod.x, comp12$frac_mod.y, method="spearman")
comp12_R2 <- comp12_R^2

c12<-ggplot(comp12, aes(x = frac_mod.x, y = frac_mod.y)) +
  geom_point(color = "#5979bc", alpha = 0.5) +  # Scatter plot with blue points
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +  # x = y line (y = x)
  labs(x = "Uninfected ARPE19 #1", y = "HSV-1 12 hpi #1", 
       title = paste("Pearson r2=", round(comp12_R, 4))) + 
  theme_classic() +
  theme(plot.title = element_text(size = 20), 
        text = element_text(size = 22, colour = "black"),  # Set all text to black
        axis.line = element_line(size = 1.5, color = "black"),  # Thicken axis lines
        axis.text = element_text(size = 20, colour = "black"))  # Black axis text

compIFN <- merge(UNINF_1_tRNA_pseU_sites, IFN_tRNA_pseU_sites, by = c("start", "end", "tRNA"), all = TRUE)
compIFN<- na.omit(compIFN)
compIFN <- compIFN[, c("start", "end", "tRNA", "frac_mod.x", "frac_mod.y")]
compIFN_R <- cor(compIFN$frac_mod.x, compIFN$frac_mod.y, method="spearman")
compIFN_R2 <- compIFN_R^2

IFN<-ggplot(compIFN, aes(x = frac_mod.x, y = frac_mod.y)) +
  geom_point(color = "#5979bc", alpha = 0.5) +  # Scatter plot with blue points
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1) +  # x = y line (y = x)
  labs(x = "Uninfected ARPE19 #1", y = "HSV-1 12 hpi #1", 
       title = paste("Pearson r2=", round(comp12_R, 4))) + 
  theme_classic() +
  theme(plot.title = element_text(size = 20), 
        text = element_text(size = 22, colour = "black"),  # Set all text to black
        axis.line = element_line(size = 1.5, color = "black"),  # Thicken axis lines
        axis.text = element_text(size = 20, colour = "black"))  # Black axis text


combined_plot4 <- (un | inf) / (c6 | c12)

combined_plot4

pdf(file = "tRNA Stoichiometry changes freq0.01.pdf", width = 10, height = 10)
combined_plot4
dev.off()

combined_plot4

pdf(file = "IFN tRNA Stoichiometry changes freq0.01.pdf", width = 5, height = 5)
IFN
dev.off()





### 10% filtering analysis

IFN_PolIIItx_pseU_sites <- IFN_PolIIItx_pseU_sites %>% mutate(source = "IFNag")
UNINF_1_PolIIItx_pseU_sites <- UNINF_1_PolIIItx_pseU_sites %>% mutate(source = "ARPE 1")
UNINF_2_PolIIItx_pseU_sites <-  UNINF_2_PolIIItx_pseU_sites %>% mutate(source = "ARPE 2")
HSV1_6hpi_1_PolIIItx_pseU_sites <- HSV1_6hpi_1_PolIIItx_pseU_sites %>% mutate(source = "HSV-6-1")
HSV1_6hpi_2_PolIIItx_pseU_sites <- HSV1_6hpi_2_PolIIItx_pseU_sites %>% mutate(source = "HSV-6-2")
HSV1_12hpi_1_PolIIItx_pseU_sites <- HSV1_12hpi_1_PolIIItx_pseU_sites %>% mutate(source = "HSV-12")

merged_df <- bind_rows(IFN_PolIIItx_pseU_sites, UNINF_1_PolIIItx_pseU_sites, UNINF_2_PolIIItx_pseU_sites, HSV1_6hpi_1_PolIIItx_pseU_sites, HSV1_6hpi_2_PolIIItx_pseU_sites, HSV1_12hpi_1_PolIIItx_pseU_sites)


write.csv(merged_df, "HSV1-pseU-over10per-polIIItx.csv", row.names = FALSE)
