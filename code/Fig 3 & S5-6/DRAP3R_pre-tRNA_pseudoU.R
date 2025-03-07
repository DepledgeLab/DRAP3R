# Load necessary libraries
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)

bed <- read.table("pre-tRNA/ARPE19_UNINF_24h-4.sup-pseU_bwa_W13k6T20.SR.primary-merged.HG38.modFiltU.0.98.hg38-tRNAs.intersect.bed", 
                  header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "") %>%
    filter(V22 != ".") %>%
    rename(coverage = V10, mod_freq = V11, Nmod = V12, 
    Ncanonical = V13, Nother_mod = V14, Ndelete = V15, 
    Nfail = V16, Ndiff = V17, Nno_call = V18, V12 = V20, V13 = V21, V14 = V22) %>%
    mutate(across(c(coverage, Nmod, Ncanonical, Nother_mod, Ndiff), as.numeric))

# Filter and calculate modification frequency
bed_filtered <- bed %>%
  filter(coverage > 100, mod_freq > 5) %>%
  mutate(across(c(V29, V30), ~ str_remove(.x, ",$")))  # Remove trailing commas

# Adjust transcript IDs by removing "tRNA" prefix
bed_filtered$V14 <- str_remove(str_sub(bed_filtered$V14, start = 5), '^-')

# Process and adjust positions
bed_filtered <- bed_filtered %>%
  separate(V29, into = c("V21.1", "V21.2"), sep = ",") %>%
  separate(V30, into = c("V22.1", "V22.2"), sep = ",") %>%
  mutate(across(c(V13, V21.1, V21.2, V22.1, V22.2, V12), as.numeric),
         V2 = V2 + 1, V12 = V12 + 1, V22.1 = V22.1 + 1, V22.2 = V22.2 + 1)

# Calculate modification positions
bed_position <- bed_filtered %>%
  mutate(Mod_position = ifelse(V6 == "-", V13 - V2 + 1, V2 - V12 + 1),
         Length = V13 - V12 + 1)

# Handle introns
bed_introns <- bed_position %>%
  filter(V28 == 2) %>%
  mutate(exon1_end = V21.1, intron_length = V22.2 - V21.1 - 1)

bed_introns_pos <- bed_introns %>%
  mutate(tRNA_length_spliced = Length - intron_length)


# Adjust modification positions based on splicing
bed_position <- bed_position %>%
  left_join(bed_introns_pos %>% select(V14, tRNA_length_spliced, intron_length), by = "V14") %>%
  mutate(tRNA_length_spliced = coalesce(tRNA_length_spliced, Length),
         intron_length = coalesce(intron_length, 0)) %>%
  mutate(Mod_position = ifelse(V28 == 2 & Mod_position > (V22.2 - 1), Mod_position - intron_length, Mod_position))


# Create matrix for visualization
df1 <- bed_position %>%
  select(V14, Mod_position, tRNA_length_spliced, mod_freq, V6) %>%
  distinct() %>%
  group_by(V14, tRNA_length_spliced, V6) %>%
  summarise(position = list(1:tRNA_length_spliced), .groups = "drop") %>%
  unnest(cols = c(position)) %>%
  mutate(ids = paste0(V14, "(", V6, ")"),
         Score = 0,
         Mod_freq = 0)

# Assign modification frequencies
df1 <- df1 %>%
  left_join(bed_position %>% select(V14, Mod_position, mod_freq), 
            by = c("V14" = "V14", "position" = "Mod_position")) %>%
  mutate(mod_freq = as.numeric(mod_freq),  # Convert to numeric
         Score = ifelse(!is.na(mod_freq), 1, 0),
         Mod_freq = coalesce(mod_freq, 0))

# Generate heatmap
pdf("pseU_ARPE19_4-1_tRNAs_spliced-alt.pdf", width = 8, height = 8)
ggplot(df1, aes(x = position, y = factor(ids, levels = rev(unique(ids))), fill = Mod_freq)) +
  geom_tile(color = "white", lwd = 0, linetype = 1) +
  scale_fill_gradient(low = "white", high = "darkblue", limits = c(0, 100), name = "Modification Frequency") +
  coord_fixed() +
  ggtitle("pseU ARPE-19 #1 tRNAs") +
  xlab("tRNA position") +
  scale_x_continuous(breaks = seq(0, max(df1$position, na.rm = TRUE), by = 5)) +
  ylab("pre_tRNA ID") +
  theme(axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        axis.text.y = element_text(margin = margin(r = -20), angle = 0, size = 5),
        axis.text.x = element_text(margin = margin(t = 0), angle = 90, hjust = 1, size = 5))
dev.off()


Glu <- df1[grep("Glu", df1$ids),]
Glu$position <- Glu$position# + 1

# Generate Glu specific heatmant
pdf("pseU_ARPE19_4-1_GLUonly_tRNAs_spliced.pdf", width = 12, height = 6)
ggplot(Glu, aes(x = position, y = factor(ids, levels = rev(unique(ids))), fill = Mod_freq)) +
  geom_tile(color = "white", lwd = 0, linetype = 1) +
  scale_fill_gradient(low = "white", high = "darkblue", limits = c(0, 100), name = "Modification Frequency") +
  coord_fixed() +
  ggtitle("pseU ARPE-19 #1 Glu tRNAs") +
  xlab("position") +
  scale_x_continuous(breaks = seq(5, max(Glu$position, na.rm = TRUE), by = 5)) +
  ylab("pre_tRNA ID") +
  theme(axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        axis.text.y = element_text(margin = margin(r = -20), angle = 0, size = 10),
        axis.text.x = element_text(margin = margin(t = 0), angle = 90, hjust = 1, size = 14))
dev.off()




p13 <- df1[grep("13", df1$position),]
p20 <- df1[grep("21", df1$position),]
p27 <- df1[grep("27|28", df1$position),]
p31 <- df1[grep("31|32", df1$position),]
p35 <- df1[grep("35|36", df1$position),]
p55 <- df1[grep("54|55|56", df1$position),]
p65 <- df1[grep("64|65|67", df1$position),]

# Add a category column to each subset
p13$category <- "p 13"
p20$category <- "p 20"
p27$category <- "p 27"
p31$category <- "p 31"
p35$category <- "p 35"
p55$category <- "p 55"
p65$category <- "p 65"

# Combine all subsets into a single dataframe
combined_df <- rbind(p13, p20, p27, p31, p35, p55, p65)
combined_df <- combined_df[!is.na(combined_df$mod_freq) & combined_df$mod_freq >= 5, ]

# Boxplot
modplot<-ggplot(combined_df, aes(x = category, y = mod_freq)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, color = "darkblue") +  # Remove fill, add black outline
  geom_jitter(aes(fill = mod_freq), width = 0.2, size = 2, alpha = 0.7, 
              shape = 21, color = "black", stroke = 0.8) + 
  scale_fill_gradient(low = "white", high = "darkblue") +  # Color points by value
  theme_minimal() +
  labs(x = "", 
       y = "stoichiometry", 
       color = "stoichiometry") +  
  theme(
    panel.grid.major = element_blank(),   # Remove major grid lines
    panel.grid.minor = element_blank(),   # Remove minor grid lines
    panel.border = element_blank(),       # Remove box around plot
    axis.line = element_line(color = "black", linewidth = 0.75),  # Thicker axis lines
    axis.title = element_text(size = 18),  # Increase axis title text size
    axis.text = element_text(color = "black", size = 18),   # Increase axis tick labels text size
    plot.title = element_text(size = 18, face = "bold"),  # Increase plot title size
    legend.position = "none"              # Remove the legend
  )


pdf("pseu_modlevel_pretRNA_DRAP3R_ARPE19-4.pdf", width = 8, height = 4)
modplot
dev.off()




