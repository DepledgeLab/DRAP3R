library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

### Set basic filtering parameters

COV <- 20
FREQ <- 10

# Step 1: Group file names based on their types
### tRNAs
file_groups <- list(
  m6A_sites = list(
    files = list(
      "IntersectMods-m6A/ARPE19_UNINF_24h-4.sup-m6A_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed",
      "IntersectMods-m6A/ARPE19_UNINF_24h-6.sup-m6A_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed",
      "IntersectMods-m6A/CRO-AP5-1.sup-m6A_bwa_W13k6T20.SR.primary-merged.HG38_EBV_KSHV.modFiltU.0.98.hg38-tRNAs.intersect.bed",
      "IntersectMods-m6A/NHDF-uninf-1.sup-m6A_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed"
    ),
    variable_names = list(
      "ARPE19_UNINF_24h_rep4_tRNA_m6A_sites",
      "ARPE19_UNINF_24h_rep6_tRNA_m6A_sites",
      "CRO_AP5_1_rep1_tRNA_m6A_sites",
      "NHDF_uninf_1_rep1_tRNA_m6A_sites"
    )
  ),
  pseU_sites = list(
    files = list(
      "IntersectMods-pseU/ARPE19_UNINF_24h-4.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed",
      "IntersectMods-pseU/ARPE19_UNINF_24h-6.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed",
      "IntersectMods-pseU/CRO-AP5-1.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38_EBV_KSHV.modFiltU.0.98.hg38-tRNAs.intersect.bed",
      "IntersectMods-pseU/NHDF-uninf-1.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed"
    ),
    variable_names = list(
      "ARPE19_UNINF_24h_rep4_tRNA_pseU_sites",
      "ARPE19_UNINF_24h_rep6_tRNA_pseU_sites",
      "CRO_AP5_1_rep1_tRNA_pseU_sites",
      "NHDF_uninf_1_rep1_tRNA_pseU_sites"
    )
  ),
  PolIIItx_sites = list(
    files = list(
      "IntersectMods-m6A/ARPE19_UNINF_24h-4.sup-m6A_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "IntersectMods-m6A/ARPE19_UNINF_24h-6.sup-m6A_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "IntersectMods-m6A/CRO-AP5-1.sup-m6A_bwa_W13k6T20.SR.primary-merged.HG38_EBV_KSHV.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "IntersectMods-m6A/NHDF-uninf-1.sup-m6A_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "IntersectMods-m6A/RN7SK-IVT-1.sup-m6A_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed"
  ),
    variable_names = list(
      "ARPE19_UNINF_24h_rep4_PolIIItx_m6A_sites",
      "ARPE19_UNINF_24h_rep6_PolIIItxA_m6A_sites",
      "CRO_AP5_1_rep1_PolIIItx_m6A_sites",
      "NHDF_uninf_1_rep1_PolIIItx_m6A_sites",
      "RN7SK-IVT-1_PolIIItx_m6A_sites"
    )
  ),
  PolIIItx_pseU_sites = list(
    files = list(
      "IntersectMods-pseU/ARPE19_UNINF_24h-4.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "IntersectMods-pseU/ARPE19_UNINF_24h-6.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "IntersectMods-pseU/CRO-AP5-1.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38_EBV_KSHV.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "IntersectMods-pseU/NHDF-uninf-1.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "IntersectMods-pseU/RN7SK-IVT-1.sup-pseU_bwa_W13k6T20.SR.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed"
      
    ),
    variable_names = list(
      "ARPE19_UNINF_24h_rep4_PolIIItx_pseU_sites",
      "ARPE19_UNINF_24h_rep6_PolIIItx_pseU_sites",
      "CRO_AP5_1_rep1_PolIIItx_pseU_sites",
      "NHDF_uninf_1_rep1_PolIIItx_pseU_sites",
      "RN7SK-IVT-1_PolIIItx_pseU_sites"
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
#      filter(valid_cov >= COV, frac_mod > FREQ, N_diff <= valid_cov, N_fail <= valid_cov/2)
      filter(valid_cov >= COV, frac_mod > FREQ, N_diff <= valid_cov)#, N_fail <= valid_cov*0.9)
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

### REPRESENTATIVE PLOTS USING UNINF2

rep1<-ggplot() +
#  geom_histogram(data = ARPE19_UNINF_24h_rep4_tRNA_pseU_sites, aes(x = frac_mod), 
#                 binwidth = 2, fill = "#5979bc", alpha = 0.9, color = "black", boundary = 0, linewidth = 0.75, linetype = "solid") +  # First histogram
  geom_histogram(data = ARPE19_UNINF_24h_rep4_tRNA_m6A_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "red", alpha = 0.9, color = "black", boundary = 0, linewidth = 0.75, linetype = "solid") +  # Second histogram
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Modification stoichiometry", y = "Count") + #, title = "ARPE19 rep1 pre-tRNA") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),  # Title size and centering
    text = element_text(size = 24, color = "black"),
    axis.text = element_text(size = 24),  # Axis text size
    axis.line = element_line(size = 1),
    axis.title = element_text(size = 24)  # Axis title size
  ) +
  ylim(0, 600) # 600 for main plot

rep2<-ggplot() +
#  geom_histogram(data = ARPE19_UNINF_24h_rep4_PolIIItx_pseU_sites, aes(x = frac_mod), 
#                 binwidth = 2, fill = "#5979bc", alpha = 0.9, color = "black", boundary = 0, linewidth = 0.75, linetype = "solid") +  # First histogram
  geom_histogram(data = ARPE19_UNINF_24h_rep4_PolIIItx_m6A_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "red", alpha = 0.9, color = "black", boundary = 0, linewidth = 0.75, linetype = "solid") +  # Second histogram
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Modification stoichiometry", y = "Count") + #, title = "ARPE19 rep1 Pol III tx") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),  # Title size and centering
    text = element_text(size = 24, color = "black"),
    axis.text = element_text(size = 24),  # Axis text size
    axis.line = element_line(size = 1),
    axis.title = element_text(size = 24)  # Axis title size
  ) +
  ylim(0, 50) #80 for main plot

combined_plot <- rep1 + rep2 + plot_layout(ncol = 1)

pdf(file = "FreqHistogramMain.pdf", width = 5, height = 10)
combined_plot
dev.off()


### INIDIVIDUAL PLOTS

#UNINF1

uninf1m<-ggplot() +
  geom_histogram(data = ARPE19_UNINF_24h_rep6_tRNA_pseU_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "#5979bc", alpha = 0.9, color = "black", boundary = 0) +  # First histogram
  geom_histogram(data = ARPE19_UNINF_24h_rep6_tRNA_m6A_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "red", alpha = 0.9, color = "black", boundary = 0) +  # Second histogram
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Modification stoichiometry", y = "Count") + #, title = "ARPE19 rep2 pre-tRNA") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 24),  # Axis text size
    axis.title = element_text(size = 24)  # Axis title size
  ) +
  ylim(0, 600)

uninf1p<-ggplot() +
  geom_histogram(data = ARPE19_UNINF_24h_rep4_PolIIItx_pseU_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "#5979bc", alpha = 0.9, color = "black", boundary = 0) +  # First histogram
  geom_histogram(data = ARPE19_UNINF_24h_rep4_PolIIItx_m6A_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "red", alpha = 0.9, color = "black", boundary = 0) +  # Second histogram
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

#UNINF3  FIX THE TWO BELOW

uninf3m<-ggplot() +
  geom_histogram(data = NHDF_uninf_1_rep1_tRNA_pseU_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "#5979bc", alpha = 0.9, color = "black", boundary = 0) +  # First histogram
  geom_histogram(data = NHDF_uninf_1_rep1_tRNA_m6A_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "red", alpha = 0.9, color = "black", boundary = 0) +  # Second histogram
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Modification stoichiometry", y = "Count") + #, title = "NHDF pre-tRNA") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 24),  # Axis text size
    axis.title = element_text(size = 24)  # Axis title size
  ) +
  ylim(0, 600)

uninf3p<-ggplot() +
  geom_histogram(data = NHDF_uninf_1_rep1_PolIIItx_pseU_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "#5979bc", alpha = 0.9, color = "black", boundary = 0) +  # First histogram
  geom_histogram(data = NHDF_uninf_1_rep1_PolIIItx_m6A_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "red", alpha = 0.9, color = "black", boundary = 0) +  # Second histogram
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Modification stoichiometry", y = "Count") + #, title = "NHDF Pol III tx") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 24),  # Axis text size
    axis.title = element_text(size = 24)  # Axis title size
  ) +
  ylim(0, 50)

#CRO
inf2m<-ggplot() +
  geom_histogram(data = CRO_AP5_1_rep1_tRNA_pseU_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "#5979bc", alpha = 0.9, color = "black", boundary = 0) +  # First histogram
  geom_histogram(data = CRO_AP5_1_rep1_tRNA_m6A_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "red", alpha = 0.9, color = "black", boundary = 0) +  # Second histogram
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Modification stoichiometry", y = "Count") + #, title = "CRO-AP5 pre-tRNA") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 24),  # Axis text size
    axis.title = element_text(size = 24)  # Axis title size
  ) +
    ylim(0, 600)

inf2p<-ggplot() +
  geom_histogram(data = CRO_AP5_1_rep1_PolIIItx_pseU_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "#5979bc", alpha = 0.9, color = "black", boundary = 0) +  # First histogram
  geom_histogram(data = CRO_AP5_1_rep1_PolIIItx_m6A_sites, aes(x = frac_mod), 
                 binwidth = 2, fill = "red", alpha = 0.9, color = "black", boundary = 0) +  # Second histogram
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Modification stoichiometry", y = "Count") + #, title = "CRO-AP5 Pol III tx") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 24),  # Axis text size
    axis.title = element_text(size = 24)  # Axis title size
  ) +
  ylim(0, 50)

#IVT
IVT2p<-ggplot() +
  geom_histogram(data = `RN7SK-IVT-1_PolIIItx_pseU_sites`, aes(x = frac_mod), 
                 binwidth = 2, fill = "#5979bc", alpha = 0.9, color = "black", boundary = 0) +  # First histogram
  geom_histogram(data = `RN7SK-IVT-1_PolIIItx_m6A_sites`, aes(x = frac_mod), 
                 binwidth = 2, fill = "red", alpha = 0.9, color = "black", boundary = 0) +  # Second histogram
  scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  labs(x = "Modification stoichiometry", y = "Count") + #, title = "CRO-AP5 Pol III tx") +  # Axis labels and title
  theme_classic() +
  theme(axis.text=element_text(colour="black")) +
  theme(
    plot.title = element_text(size = 24, hjust = 0.5),  # Title size and centering
    axis.text = element_text(size = 24),  # Axis text size
    axis.title = element_text(size = 24)  # Axis title size
  ) +
  ylim(0, 20)



combined_plot2 <- rep1 + uninf1m + uninf3m + inf2m + rep2 + uninf1p + uninf3p + inf2p + plot_layout(ncol = 5)
  
pdf(file = "FreqHistogramSupplementary.pdf", width = 25, height = 10)
combined_plot2
dev.off()

pdf(file = "FreqHistogramSupplementaryIVT.pdf", width = 5, height = 5)
IVT2p
dev.off()

