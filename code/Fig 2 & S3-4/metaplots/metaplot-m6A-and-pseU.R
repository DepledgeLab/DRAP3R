library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(patchwork)

setwd("C:/Users/depledgd/Desktop/code/Fig 2/metaplots")


### Set basic filtering parameters

COV <- 20
FREQ <- 1 #10

# Step 1: Group file names based on their types
file_groups <- list(
  m6A_sites = list(
    files = list(
      "../stochiometry histograms/IntersectMods-m6A/ARPE19_UNINF_24h-4.sup-m6A_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed",
      "../stochiometry histograms/IntersectMods-m6A/ARPE19_UNINF_24h-6.sup-m6A_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed",
      "../stochiometry histograms/IntersectMods-m6A/CRO-AP5-1.sup-m6A_bwa_W13k6T20.SR.primary-merged.HG38_EBV_KSHV.modFiltU.0.98.hg38-tRNAs.intersect.bed",
      "../stochiometry histograms/IntersectMods-m6A/NHDF-uninf-1.sup-m6A_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed"
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
      "../stochiometry histograms/IntersectMods-pseU/ARPE19_UNINF_24h-4.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed",
      "../stochiometry histograms/IntersectMods-pseU/ARPE19_UNINF_24h-6.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed",
      "../stochiometry histograms/IntersectMods-pseU/CRO-AP5-1.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38_EBV_KSHV.modFiltU.0.98.hg38-tRNAs.intersect.bed",
      "../stochiometry histograms/IntersectMods-pseU/NHDF-uninf-1.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed"
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
      "../stochiometry histograms/IntersectMods-m6A/ARPE19_UNINF_24h-4.sup-m6A_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "../stochiometry histograms/IntersectMods-m6A/ARPE19_UNINF_24h-6.sup-m6A_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "../stochiometry histograms/IntersectMods-m6A/CRO-AP5-1.sup-m6A_bwa_W13k6T20.SR.primary-merged.HG38_EBV_KSHV.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "../stochiometry histograms/IntersectMods-m6A/NHDF-uninf-1.sup-m6A_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "../stochiometry histograms/IntersectMods-m6A/RN7SK-IVT-1.sup-m6A_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed"
    ),
    variable_names = list(
      "ARPE19_UNINF_24h_rep4_PolIIItx_m6A_sites",
      "ARPE19_UNINF_24h_rep6_PolIIItx_m6A_sites",
      "CRO_AP5_1_rep1_PolIIItx_m6A_sites",
      "NHDF_uninf_1_rep1_PolIIItx_m6A_sites",
      "RN7SK-IVT-1_PolIIItx_m6A_sites"
    )
  ),
  PolIIItx_pseU_sites = list(
    files = list(
      "../stochiometry histograms/IntersectMods-pseU/ARPE19_UNINF_24h-4.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "../stochiometry histograms/IntersectMods-pseU/ARPE19_UNINF_24h-6.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "../stochiometry histograms/IntersectMods-pseU/CRO-AP5-1.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38_EBV_KSHV.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed",
      "../stochiometry histograms/IntersectMods-pseU/NHDF-uninf-1.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.PolIIItx.primaryOnly.intersect.bed"
      
    ),
    variable_names = list(
      "ARPE19_UNINF_24h_rep4_PolIIItx_pseU_sites",
      "ARPE19_UNINF_24h_rep6_PolIIItx_pseU_sites",
      "CRO_AP5_1_rep1_PolIIItx_pseU_sites",
      "NHDF_uninf_1_rep1_PolIIItx_pseU_sites"
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

#### Combine and plot tRNA m6A modification data
combined.metagene.coord <- c(ARPE19_UNINF_24h_rep4_tRNA_m6A_sites$result, ARPE19_UNINF_24h_rep6_tRNA_m6A_sites$result, NHDF_uninf_1_rep1_tRNA_m6A_sites$result, CRO_AP5_1_rep1_tRNA_m6A_sites$result)
mod <- c(rep("tRNA_arpe19_rep1", length(ARPE19_UNINF_24h_rep4_tRNA_m6A_sites$result)), 
         rep("tRNA_arpe19_rep2", length(ARPE19_UNINF_24h_rep6_tRNA_m6A_sites$result)),
         rep("tRNA_nhdf_rep1", length(NHDF_uninf_1_rep1_tRNA_m6A_sites$result)),
         rep("tRNA_croap5_rep1", length(CRO_AP5_1_rep1_tRNA_m6A_sites$result))) 
tRNA_m6A <- data.frame(combined.metagene.coord, mod)

plot_tRNA_m6A<-ggplot(tRNA_m6A) + 
  geom_density(aes(x = combined.metagene.coord, colour = mod), show.legend = FALSE, size=1.5) + 
  xlim(0, 1) + ylim(0,2.5) + 
  theme_classic() + 
  scale_color_manual(values = c(
    "tRNA_arpe19_rep1" = "cadetblue", 
    "tRNA_arpe19_rep2" = "cadetblue", 
    "tRNA_nhdf_rep1" = "darkorange2", 
    "tRNA_croap5_rep1" = "goldenrod1"
  )) + 
  theme(
    legend.position = "none", 
    axis.text = element_text(colour = "black", size = 24), 
    plot.title = element_text(size = 24, hjust = 0.5), 
    text = element_text(size = 24, color = "black"), 
    axis.line = element_line(size = 1), 
    axis.title = element_text(size = 24)
  )


#### Combine and plot tRNA pseU modification data
combined.metagene.coord <- c(ARPE19_UNINF_24h_rep4_tRNA_pseU_sites$result, ARPE19_UNINF_24h_rep6_tRNA_pseU_sites$result, NHDF_uninf_1_rep1_tRNA_pseU_sites$result, CRO_AP5_1_rep1_tRNA_pseU_sites$result)
mod <- c(rep("tRNA_arpe19_rep1", length(ARPE19_UNINF_24h_rep4_tRNA_pseU_sites$result)), 
         rep("tRNA_arpe19_rep2", length(ARPE19_UNINF_24h_rep6_tRNA_pseU_sites$result)),
         rep("tRNA_nhdf_rep1", length(NHDF_uninf_1_rep1_tRNA_pseU_sites$result)),
         rep("tRNA_croap5_rep1", length(CRO_AP5_1_rep1_tRNA_pseU_sites$result))) 
tRNA_pseU <- data.frame(combined.metagene.coord, mod)

plot_tRNA_pseU<-ggplot(tRNA_pseU) + 
  geom_density(aes(x = combined.metagene.coord, colour = mod), show.legend = FALSE, size=1.5) + 
  xlim(0, 1) + ylim(0,2.5) + 
  theme_classic() + 
  scale_color_manual(values = c(
    "tRNA_arpe19_rep1" = "cadetblue", 
    "tRNA_arpe19_rep2" = "cadetblue",
    "tRNA_nhdf_rep1" = "darkorange2", 
    "tRNA_croap5_rep1" = "goldenrod1"
  )) + 
  theme(
    legend.position = "none", 
    axis.text = element_text(colour = "black", size = 24), 
    plot.title = element_text(size = 24, hjust = 0.5), 
    text = element_text(size = 24, color = "black"), 
    axis.line = element_line(size = 1), 
    axis.title = element_text(size = 24)
  )

#### Combine and plot PolIIItx m6A modification data
combined.metagene.coord <- c(ARPE19_UNINF_24h_rep4_PolIIItx_m6A_sites$result, ARPE19_UNINF_24h_rep6_PolIIItx_m6A_sites$result, NHDF_uninf_1_rep1_PolIIItx_m6A_sites$result, CRO_AP5_1_rep1_PolIIItx_m6A_sites$result)
mod <- c(rep("PolIIItx_arpe19_rep1", length(ARPE19_UNINF_24h_rep4_PolIIItx_m6A_sites$result)), 
         rep("PolIIItx_arpe19_rep2", length(ARPE19_UNINF_24h_rep6_PolIIItx_m6A_sites$result)),
         rep("PolIIItx_nhdf_rep1", length(NHDF_uninf_1_rep1_PolIIItx_m6A_sites$result)),
         rep("PolIIItx_croap5_rep1", length(CRO_AP5_1_rep1_PolIIItx_m6A_sites$result))) 
PolIIItx_m6A <- data.frame(combined.metagene.coord, mod)

plot_PolIIItx_m6A<-ggplot(PolIIItx_m6A) + 
  geom_density(aes(x = combined.metagene.coord, colour = mod), show.legend = FALSE, size=1.5) + 
  xlim(0, 1) + ylim(0,2.5) + 
  theme_classic() + 
  scale_color_manual(values = c(
    "PolIIItx_arpe19_rep1" = "cadetblue", 
    "PolIIItx_arpe19_rep2" = "cadetblue", 
    "PolIIItx_arpe19_ifn_rep1" = "cadetblue2", 
    "PolIIItx_nhdf_rep1" = "darkorange2", 
    "PolIIItx_croap5_rep1" = "goldenrod1"
  )) + 
  theme(
    legend.position = "none", 
    axis.text = element_text(colour = "black", size = 24), 
    plot.title = element_text(size = 24, hjust = 0.5), 
    text = element_text(size = 24, color = "black"), 
    axis.line = element_line(size = 1), 
    axis.title = element_text(size = 24)
  )


#### Combine and plot PolIIItx pseU modification data
combined.metagene.coord <- c(ARPE19_UNINF_24h_rep4_PolIIItx_pseU_sites$result, ARPE19_UNINF_24h_rep6_PolIIItx_pseU_sites$result, NHDF_uninf_1_rep1_PolIIItx_pseU_sites$result, CRO_AP5_1_rep1_PolIIItx_pseU_sites$result)
mod <- c(rep("PolIIItx_arpe19_rep1", length(ARPE19_UNINF_24h_rep4_PolIIItx_pseU_sites$result)), 
         rep("PolIIItx_arpe19_rep2", length(ARPE19_UNINF_24h_rep6_PolIIItx_pseU_sites$result)),
         rep("PolIIItx_nhdf_rep1", length(NHDF_uninf_1_rep1_PolIIItx_pseU_sites$result)),
         rep("PolIIItx_croap5_rep1", length(CRO_AP5_1_rep1_PolIIItx_pseU_sites$result))) 
PolIIItx_pseU <- data.frame(combined.metagene.coord, mod)

plot_PolIIItx_pseU<-ggplot(PolIIItx_pseU) + 
  geom_density(aes(x = combined.metagene.coord, colour = mod), show.legend = FALSE, size=1.5) + 
  xlim(0, 1) + ylim(0,2.5) + 
  theme_classic() + 
  scale_color_manual(values = c(
    "PolIIItx_arpe19_rep1" = "cadetblue", 
    "PolIIItx_arpe19_rep2" = "cadetblue",
    "PolIIItx_nhdf_rep1" = "darkorange2", 
    "PolIIItx_croap5_rep1" = "goldenrod1"
  )) + 
  theme(
    legend.position = "none", 
    axis.text = element_text(colour = "black", size = 24), 
    plot.title = element_text(size = 24, hjust = 0.5), 
    text = element_text(size = 24, color = "black"), 
    axis.line = element_line(size = 1), 
    axis.title = element_text(size = 24)
  )


final <- (plot_tRNA_m6A | plot_PolIIItx_m6A) / (plot_tRNA_pseU | plot_PolIIItx_pseU)

print(final)

pdf(file = "metaplot-freq1.pdf", width = 10, height = 10) #5
final
dev.off()

#write.csv(ARPE19_UNINF_24h_rep4_tRNA_pseU_sites, "ARPE19_UNINF_24h_rep4_tRNA_pseU_sites.csv")
#write.csv(ARPE19_UNINF_24h_rep4_tRNA_m6A_sites, "ARPE19_UNINF_24h_rep4_tRNA_m6A_sites.csv")
#write.csv(ARPE19_UNINF_24h_rep4_PolIIItx_pseU_sites, "ARPE19_UNINF_24h_rep4_PolIIItx_pseU_sites.csv")
#write.csv(ARPE19_UNINF_24h_rep4_PolIIItx_m6A_sites, "ARPE19_UNINF_24h_rep4_PolIIItx_m6A_sites.csv")

### ANALYSIS OF RN7SK


