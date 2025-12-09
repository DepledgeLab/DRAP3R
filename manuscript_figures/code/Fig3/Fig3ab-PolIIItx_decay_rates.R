library(tidyverse)

#setwd("")

# Read decay files
control <- read_csv("ARPE19_uninf-1_PolIIItx_primary_decay_minReads100.csv") %>%
  rename(k_control = "decay_rate_k")

control2 <- read_csv("ARPE19_uninf-2_PolIIItx_primary_decay_minReads100.csv") %>%
  rename(k_control2 = "decay_rate_k")

knockdown <- read_csv("ARPE19_siLa-1_PolIIItx_primary_decay_minReads100.csv") %>%
  rename(k_knockdown = "decay_rate_k")

# Merge by tRNA ID
mergedCont <- inner_join(control, control2, by = "tRNA_id")

mergedTest <- inner_join(control2, knockdown, by = "tRNA_id")

mergedtmp <- full_join(control, control2, by = "tRNA_id")
mergedfull <- full_join(mergedtmp, knockdown, by = "tRNA_id")


write.csv(mergedfull, "PolIII_k_values.csv")


p3_control<-ggplot(mergedfull, aes(x = k_control, y = k_control2)) +
  geom_point(size = 2, color = "#f35120") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Pol III Decay Rate Comparison",
    x = "ARPE-19 #1 Decay Rate (k)",
    y = "ARPE-19 #2 Decay Rate (k)"
  ) +
  xlim(-0.05, 0.1) +
  ylim(-0.05, 0.1) +
  theme_classic() +
  theme(axis.text=element_text(colour="black")) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  theme(text = element_text(size = 22, colour = "black"), axis.line=element_line(size=0.75), legend.position = "none") 




p3_comp<-ggplot(mergedfull, aes(x = k_control2, y = k_knockdown)) +
  geom_point(size = 2, color = "#f35120") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Pol III Decay Rate Comparison",
    x = "ARPE-19 #1 Decay Rate (k)",
    y = "siLA #1 Decay Rate (k)"
  ) +
  xlim(-0.05, 0.1) +
  ylim(-0.05, 0.1) +
  theme_classic() +
  theme(axis.text=element_text(colour="black")) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  theme(text = element_text(size = 22, colour = "black"), axis.line=element_line(size=0.75), legend.position = "none") 


wilcox.test(mergedTest$k_knockdown, mergedTest$k_control2, paired = TRUE)

wilcox.test(mergedCont$k_control, mergedCont$k_control2, paired = TRUE)

pdf(file = "PolIII_control_decay.pdf", width = 5, height = 5)
p3_control
dev.off()

pdf(file = "PolIII_comp_decay.pdf", width = 5, height = 5)
p3_comp
dev.off()






