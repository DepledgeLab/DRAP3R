# Load necessary libraries
library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(purrr)

# Set working directory
#setwd("")

#bed <- read.table("tRNA-splintdb/ARPE19_nanotRNA-uninf-2.sup-pseU_bwa_W13k6T20.SR.merged.primary.splint_tRNAdb.modFilt.0.98.sorted.bed", 
bed <- read.table("tRNA-splintdb/ARPE19_uninf-2.sup-pseU_bwa_W13k6T20.SR.merged.primary.splint_tRNAdb.modFilt.0.98.sorted.bed", 
                                    
                                    header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "") %>%
  rename(coverage = V10, mod_freq = V11, Nmod = V12, 
         Ncanonical = V13, Nother_mod = V14, Ndelete = V15, 
         Nfail = V16, Ndiff = V17, Nno_call = V18) %>%
  mutate(across(c(coverage, Nmod, Ncanonical, Nother_mod, Ndiff), as.numeric),
         V1 = as.character(trimws(V1)))  # Ensure V1 is character

info <- read.table("tRNA-splintdb/mature-tRNA-isodecoders-wsplintadapt-new.list.txt", 
                   header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "") %>%
  mutate(V1 = as.character(trimws(V1)))  # Ensure V1 is character

# Merge, keeping only the 3rd column from info
bed_filtered <- bed %>%
  left_join(info %>% select(V1, V3), by = "V1") %>%
  filter(coverage > 100, mod_freq > 5)

bed_filtered$V1 <- str_remove(str_sub(bed_filtered$V1, start = 19), '^-')

# Sanity check for missing matches
#setdiff(bed$V1, info$V1) 

# Calculate modification positions
bed_position <- bed_filtered %>%
  mutate(Mod_position = V3.x, Length = V3.y)

# Create matrix for visualization
df1 <- bed_position %>%
  select(V1, Mod_position, Length) %>%
  distinct() %>%
  mutate(position = map(Length, ~seq(pmax(24, 1), pmax(1, .x - 33)))) %>%  #remove splint adapter sequences
  unnest(cols = c(position)) %>%
  mutate(ids = V1, 
         Score = 0, 
         Mod_freq = 0) 


# Assign modification frequencies
df1 <- df1 %>%
  left_join(bed_position %>% select(V1, Mod_position, mod_freq), 
            by = c("V1" = "V1", "position" = "Mod_position")) %>%
  mutate(mod_freq = as.numeric(mod_freq),  # Convert to numeric
         Score = ifelse(!is.na(mod_freq), 1, 0),
         Mod_freq = coalesce(mod_freq, 0),
         renumbered_position = dense_rank(position))  # Renumber positions from 1 to max

# Generate heatmap
pdf("nanotRNA-ARPE19_UNINF_24h-4-splint-tRNADB.pdf", width = 11.7, height = 18)
#pdf("DRAP3R-ARPE19_UNINF_24h-4-splint-tRNADB.pdf", width = 11.7, height = 18)
ggplot(df1, aes(x = renumbered_position, y = factor(ids, levels = rev(unique(ids))), fill = Mod_freq)) +
  geom_tile(color = "white", lwd = 0, linetype = 1) +
  scale_fill_gradient(low = "white", high = "darkblue", limits = c(0, 100), name = "Modification Frequency") +
  coord_fixed() +
  ggtitle("nanotRNA-ARPE19_UNINF_24h-4-splint-tRNADB") +
#  ggtitle("DRAP3R-ARPE19_UNINF_24h-4-splint-tRNADB") +
  xlab("position") +
  scale_x_continuous(breaks = seq(0, max(df1$position, na.rm = TRUE), by = 5)) +
  ylab("tRNA isodecoder ID") +
  theme(axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        axis.text.y = element_text(margin = margin(r = -20), angle = 0, size = 5),
        axis.text.x = element_text(margin = margin(t = 0), angle = 90, hjust = 1, size = 5))
dev.off()


Glu <- df1[grep("Glu", df1$ids),]
Glu <- Glu[grep("Glu-TTC-3", Glu$ids,invert=TRUE),]
Glu$position <- Glu$position - 24 #+ 1

# Generate Glu specific heatmant
#pdf("nano-ARPE19_UNINF_24h-4-splint-tRNADB-GluOnly.pdf", width = 12, height = 6)
pdf("DRAP3R-ARPE19_UNINF_24h-4-splint-tRNADB-GluOnly.pdf", width = 12, height = 6)
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

