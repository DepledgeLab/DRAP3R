library(ggplot2)
library(dplyr)
library(patchwork)

setwd("")

bwasw <- read.table("ARPE19_uninf-1_bwa_bwasw.readlengths.txt", sep="\t", header=F, stringsAsFactors=FALSE)
bwasw_opts <- read.table("ARPE19_uninf-1_bwa_bwasw_opts.readlengths.txt", sep="\t", header=F, stringsAsFactors=FALSE)
W13k6ont2d <- read.table("ARPE19_uninf-1_bwa_W13k6ont2d.readlengths.txt", sep="\t", header=F, stringsAsFactors=FALSE)
W13k6T10 <- read.table("ARPE19_uninf-1_bwa_W13k6T10.readlengths.txt", sep="\t", header=F, stringsAsFactors=FALSE)
W13k6T20 <- read.table("ARPE19_uninf-1_bwa_W13k6T20.readlengths.txt", sep="\t", header=F, stringsAsFactors=FALSE)
W9k5T10 <- read.table("ARPE19_uninf-1_bwa_W9k5T10.readlengths.txt", sep="\t", header=F, stringsAsFactors=FALSE)
map_ont <- read.table("ARPE19_uninf-1_minimap2_map-ont.readlengths.txt", sep="\t", header=F, stringsAsFactors=FALSE)
splice_k15 <- read.table("ARPE19_uninf-1_minimap2_splice_k15.readlengths.txt", sep="\t", header=F, stringsAsFactors=FALSE)

bwasw$dataset <- "opt 1"
bwasw_opts$dataset <- "opt 2"
W13k6ont2d$dataset <- "-W13 -k6"
W13k6T10$dataset <- "-W13 -k6 -T10"
W13k6T20$dataset <- "-W13 -k6 -T20"
W9k5T10$dataset <- "-W9 -k5 -T10"
map_ont$dataset <- "-ax map-ont -k15"
splice_k15$dataset <- "-ax splice -k15"

tmp1<-rbind(bwasw, bwasw_opts)
tmp2<-rbind(tmp1,W13k6ont2d)
tmp3<-rbind(tmp2,W13k6T10)
tmp4<-rbind(tmp3,W13k6T20)
tmp5<-rbind(tmp4,W9k5T10)
tmp6<-rbind(tmp5,map_ont)
merged<-rbind(tmp6,splice_k15)

### Set up colours and plotting order for violin plot
conditions<-c("opt 1","opt 2","-W13 -k6","-W13 -k6 -T10","-W13 -k6 -T20","-W9 -k5 -T10","-ax map-ont -k15","-ax splice -k15")
condition_colours<-c("darkslategray4","darkslategray4","darkseagreen2","darkseagreen2","darkseagreen2","darkseagreen2","goldenrod2","goldenrod2")
merged$dataset <-factor(merged$dataset, levels=conditions)
merged$V1<-as.numeric(merged$V1)

### Generate violin plot
p<-ggplot(merged, aes(x=dataset, y=V1, fill=dataset)) + geom_violin(linewidth=1) + ylim(0, 500) + theme_classic() +
    labs(y = "read length") + scale_fill_manual(values = condition_colours) + theme(legend.position = "none") +
  theme(plot.title = element_text(size = 14), 
        text = element_text(size = 18, colour = "black"),  # Set all text to black
        axis.line = element_line(size = 1, color = "black"),  # Thicken axis lines
        axis.text = element_text(size = 18, colour = "black"),  # Black axis text
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        axis.ticks = element_line(size = 1, color = "black"))  # Thicker axis ticks

print(p)

pdf(file = "parameterSweepReadLengthViolinPlot.pdf", width = 7, height = 6)
p
dev.off()



