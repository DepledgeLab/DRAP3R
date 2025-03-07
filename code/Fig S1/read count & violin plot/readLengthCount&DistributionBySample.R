library(ggplot2)
library(dplyr)
library(patchwork)

IVT_7SK <- read.table("RN7SK-IVT-1.sup-pseU.notrim.dorado.0.7.0.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)
IVT_EBER <- read.table("EBER2-IVT-3.sup-pseU.notrim.dorado.0.7.0.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)

IVT_7SK$dataset <- "IVT RN7SK"
IVT_EBER$dataset <- "IVT EBER2"

IVT_7SK <- tail(IVT_7SK, n = -1)    # remove header prior to merging
IVT_EBER <- tail(IVT_EBER, n = -1)

merged<-rbind(IVT_7SK, IVT_EBER)

### Set up colours and plotting order for violin plot
conditions<-c("IVT RN7SK","IVT EBER2")
condition_colours<-c("grey","grey")
merged$dataset <-factor(merged$dataset, levels=conditions)
merged$sequence_length_template<-as.numeric(merged$sequence_length_template)

### Generate violin plot
p<-ggplot(merged, aes(x=dataset, y=sequence_length_template, fill=dataset)) + geom_violin(linewidth=0.75) + ylim(0, 500) + theme_classic() +
  ylab("read length") + theme(axis.title.x = element_blank()) + scale_fill_manual(values = condition_colours) + theme(axis.text=element_text(colour="black")) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  theme(text = element_text(size = 22, colour = "black"), axis.line=element_line(size=0.75), legend.position = "none") 

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
  scale_y_continuous(labels = scales::number_format(scale = 1e-3, suffix = "K"), limits = c(0,150000)) + 
  theme(axis.text=element_text(colour="black")) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
  theme(text = element_text(size = 22, colour = "black"), axis.line=element_line(size=0.75), legend.position = "none") 

### Combine barplot (left) and violinplot (right)
combined_plot <- b + p + plot_layout(ncol = 2, widths = c(0.5, 0.5))

print(combined_plot)

pdf(file = "readLengthCount&DistributionBySample.pdf", width = 7, height = 7)
combined_plot
dev.off()

