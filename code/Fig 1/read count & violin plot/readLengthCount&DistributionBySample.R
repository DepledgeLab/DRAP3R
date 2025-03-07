library(ggplot2)
library(dplyr)
library(patchwork)

ARPE1 <- read.table("ARPE19_UNINF_24h-4.sup-pseU.notrim.dorado.0.7.0.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)
ARPE2 <- read.table("ARPE19_UNINF_24h-6.sup-pseU.notrim.dorado.0.7.0.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)
CROAP5 <- read.table("CRO-AP5-1.sup-pseU.notrim.dorado.0.7.0.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)
ARPEpA <- read.table("ARPE19_polyA_UNINF_24h-6.sup-pseU.notrim.dorado.0.7.0.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)
NHDF <- read.table("NHDF-uninf-1.sup-pseU.notrim.dorado.0.7.0.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)

ARPE1$dataset <- "ARPE-19 #1"
ARPE2$dataset <- "ARPE-19 #2"
CROAP5$dataset <- "CRO-AP5 #1"
ARPEpA$dataset <- "ARPE-19 poly(A) #1"
NHDF$dataset <- "NHDF #1"

### remove headers prior to merging
ARPE1 <- tail(ARPE1, n = -1)
ARPE2 <- tail(ARPE2, n = -1)
CROAP5 <- tail(CROAP5, n = -1)
ARPEpA <- tail(ARPEpA, n = -1)
NHDF <- tail(NHDF, n = -1)

tmp1<-rbind(ARPE1, ARPE2)
tmp2<-rbind(tmp1,CROAP5)
tmp3<-rbind(tmp2,ARPEpA)
merged<-rbind(tmp3,NHDF)

### Set up colours and plotting order for violin plot
conditions<-c("ARPE-19 poly(A) #1","ARPE-19 #1","ARPE-19 #2","NHDF #1","CRO-AP5 #1")
condition_colours<-c("slategray1","cadetblue","cadetblue","cadetblue","cadetblue")
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
  scale_y_continuous(labels = scales::number_format(scale = 1e-3, suffix = "K"), limits = c(0,3000000)) + 
  theme(axis.text=element_text(colour="black")) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  theme(text = element_text(size = 22, colour = "black"), axis.line=element_line(size=0.75), legend.position = "none") 


### Combine barplot (left) and violinplot (right)
combined_plot <- b + p + plot_layout(ncol = 2, widths = c(0.5, 0.5))

print(combined_plot)

pdf(file = "readLengthCount&DistributionBySample-v2.pdf", width = 12, height = 7)
combined_plot
dev.off()


