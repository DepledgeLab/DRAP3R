library(ggplot2)
library(dplyr)
library(patchwork)

#setwd("I:/Dropbox/DRAP3R/final_figures/Fig 4 & S7/read count & violin plot/")
setwd("C:/Users/depledgd/Dropbox/DRAP3R/final_figures/Fig 4 & S7/read count & violin plot/")

INF1 <- read.table("ARPE19_HSV1_6hpi-2.sup-pseU.notrim.dorado.0.7.0.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)
INF2 <- read.table("ARPE19_HSV1_6hpi-3.sup-pseU.notrim.dorado.0.7.0.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)
INF3 <- read.table("ARPE19-HSV1-12hpi-1.sup-pseU.notrim.dorado.0.7.0.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)
INF4 <- read.table("ARPE19-HSV1-12hpi-3.sup-pseU.notrim.dorado.0.7.0.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)
IFN1 <- read.table("ARPE19-IFNag-24h-1.sup-pseU.notrim.dorado.0.7.0.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)

INF1$dataset <- "HSV1 6hpi #1"
INF2$dataset <- "HSV1 6hpi #2"
INF3$dataset <- "HSV1 12hpi #1"
INF4$dataset <- "HSV1 12hpi #2"
IFN1$dataset <- "Mock + IFNag #1"

INF1 <- tail(INF1, n = -1)    # remove header prior to merging
INF2 <- tail(INF2, n = -1)
INF3 <- tail(INF3, n = -1)
INF4 <- tail(INF4, n = -1)
IFN1 <- tail(IFN1, n = -1) 

tmp1<-rbind(INF1, INF2)
tmp2<-rbind(INF3,tmp1)
tmp3<-rbind(INF4,tmp2)
merged<-rbind(tmp3,IFN1)

### Set up colours and plotting order for violin plot
conditions<-c("HSV1 6hpi #1","HSV1 6hpi #2","HSV1 12hpi #1","HSV1 12hpi #2","Mock + IFNag #1")
condition_colours<-c("deeppink","deeppink","deeppink4","deeppink4","cadetblue")
merged$dataset <-factor(merged$dataset, levels=conditions)
merged$sequence_length_template<-as.numeric(merged$sequence_length_template)

### Generate violin plot
p<-ggplot(merged, aes(x=dataset, y=sequence_length_template, fill=dataset)) + geom_violin(linewidth=0.75) + ylim(0, 500) + theme_classic() +
  ylab("read length") + theme(axis.title.x = element_blank()) + scale_fill_manual(values = condition_colours) + theme(axis.text=element_text(colour="black")) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  theme(text = element_text(size = 22, colour = "black"), axis.line=element_line(size=0.75), legend.position = "none") 

#labs(y = "read length")
# face = "bold"

### Generate total read counts for each dataset
count_merged <- merged %>%
  group_by(dataset) %>%
  summarise(count = n())

### Generate barplot
b<-ggplot(count_merged, aes(x = dataset, y = count, fill = dataset)) + 
  geom_bar(stat = "identity", position = position_dodge(), color = "black", width = 0.65, linewidth=0.75) + 
  theme_classic() +
  ylab("# reads") + 
  theme(axis.title.x = element_blank()) + 
  scale_fill_manual(values = condition_colours) + theme(legend.position = "none") +
  scale_y_continuous(labels = scales::number_format(scale = 1e-3, suffix = "K"), limits = c(0,4000000)) + 
  theme(axis.text=element_text(colour="black")) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  theme(text = element_text(size = 22, colour = "black"), axis.line=element_line(size=0.75), legend.position = "none") 



### Combine barplot (left) and violinplot (right)
combined_plot <- b + p + plot_layout(ncol = 2, widths = c(0.5, 0.5))

print(combined_plot)

pdf(file = "HSV1 readLengthCount&DistributionBySampleNew.pdf", width = 10, height = 7)
combined_plot
dev.off()


