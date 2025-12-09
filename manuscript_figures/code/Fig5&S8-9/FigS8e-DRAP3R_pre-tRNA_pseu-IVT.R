library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(purrr)

# Function to create a simplified sample name
make_sample_name <- function(filepath) {
  fname <- basename(filepath)
  # Take everything before the first "."
  base <- str_split(fname, "\\.", n = 2)[[1]][1]
  # Check for "0.98" or "nofilt"
  tag <- case_when(
    str_detect(fname, "0.98") ~ "0.98",
    str_detect(fname, "nofilt") ~ "nofilt",
    TRUE ~ ""
  )
  if (tag != "") {
    return(paste0(base, "_", tag))
  } else {
    return(base)
  }
}

# Function to process a single bed file and return per-position data
process_bed <- function(filepath, keep_ids) {
  bed <- read.table(filepath, 
                    header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "") %>%
    filter(V22 != ".") %>%
    rename(coverage = V10, mod_freq = V11, Nmod = V12, 
           Ncanonical = V13, Nother_mod = V14, Ndelete = V15, 
           Nfail = V16, Ndiff = V17, Nno_call = V18, 
           V12 = V20, V13 = V21, V14 = V22) %>%
    mutate(across(c(coverage, Nmod, Ncanonical, Nother_mod, Ndiff), as.numeric)) %>%
    filter(V14 %in% keep_ids)
  
  bed_filtered <- bed %>%
    filter(coverage > 100, mod_freq > 5) %>%
    mutate(across(c(V29, V30), ~ str_remove(.x, ",$")))
  
  bed_filtered$V14 <- str_remove(str_sub(bed_filtered$V14, start = 5), '^-')
  
  bed_filtered <- bed_filtered %>%
    separate(V29, into = c("V21.1", "V21.2"), sep = ",") %>%
    separate(V30, into = c("V22.1", "V22.2"), sep = ",") %>%
    mutate(across(c(V13, V21.1, V21.2, V22.1, V22.2, V12), as.numeric),
           V2 = V2 + 1, V12 = V12 + 1, V22.1 = V22.1 + 1, V22.2 = V22.2 + 1)
  
  bed_position <- bed_filtered %>%
    mutate(Mod_position = ifelse(V6 == "-", V13 - V2 + 1, V2 - V12 + 1),
           Length = V13 - V12 + 1)
  
  bed_introns <- bed_position %>%
    filter(V28 == 2) %>%
    mutate(exon1_end = V21.1, intron_length = V22.2 - V21.1 - 1)
  
  bed_introns_pos <- bed_introns %>%
    mutate(tRNA_length_spliced = Length - intron_length)
  
  bed_position <- bed_position %>%
    left_join(bed_introns_pos %>% select(V14, tRNA_length_spliced, intron_length), by = "V14") %>%
    mutate(tRNA_length_spliced = coalesce(tRNA_length_spliced, Length),
           intron_length = coalesce(intron_length, 0)) %>%
    mutate(Mod_position = ifelse(V28 == 2 & Mod_position > (V22.2 - 1), 
                                 Mod_position - intron_length, Mod_position))
  
  # Collapse to one row per sample
  df1 <- bed_position %>%
    select(Mod_position, mod_freq) %>%
    mutate(Mod_position = as.integer(Mod_position)) %>%
    group_by(Mod_position) %>%
    summarise(Mod_freq = mean(as.numeric(mod_freq), na.rm = TRUE), .groups = "drop") %>%
    mutate(sample = make_sample_name(filepath))   # <-- smart naming
  
  return(df1)
}

# === Usage ===
keep_ids <- c("tRNA-Glu-CTC-1-1")
#keep_ids <- c("tRNA-Glu-TTC-2-1")
#keep_ids <- c("tRNA-Arg-ACG-1-1")

# List all bed files
files <- list.files("pre-tRNA/", pattern = "*.bed", full.names = TRUE)

# Process all files and combine
all_data <- map_dfr(files, ~ process_bed(.x, keep_ids))

unique(all_data$sample)

# === Custom order for plotting ===
desired_order <- c(
  "CRO-AP5-1_0.98",
  "NHDF-uninf-1_0.98",  
  "ARPE19_uninf-1_0.98",
  "ARPE19_uninf-2_0.98",  
  "tRNA-IVTmix-1_0.98",  
  "tRNA-IVTmix-1"
)


all_data_complete <- all_data %>%
  group_by(sample) %>%
  complete(Mod_position = full_seq(1:max(Mod_position, na.rm = TRUE), 1), fill = list(Mod_freq = 0)) %>%
  ungroup() %>%
  mutate(sample = factor(sample, levels = desired_order))   # <-- enforce order



# === Heatmap (one row per dataset) ===
pdf(paste0("pseU_",keep_ids,"_allSamples_oneRowPerFile.pdf"), width = 10, height = 6)
ggplot(all_data_complete, aes(x = Mod_position, 
                              y = sample,   # <- just use the factor column
                              fill = Mod_freq)) +
  geom_tile(color = "white", lwd = 0, linetype = 1) +
  scale_fill_gradient(low = "white", high = "darkblue", limits = c(0, 100), name = "Modification Frequency") +
  coord_fixed() +
  ggtitle(paste0("pseU_",keep_ids)) +
  xlab("tRNA position") +
  ylab("Sample") +
  theme(axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "transparent"),
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
dev.off()
