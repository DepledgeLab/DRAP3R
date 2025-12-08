library(ggplot2)
library(dplyr)

setwd("")

TRIM <- read.table("ARPE19_uninf-1.sup.trimAdapters.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)
NOTRIM <- read.table("ARPE19_uninf-1.sup.notrim.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)

#TRIM <- read.table("ARPE19_uninf-2.sup.trimAdapters.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)
#NOTRIM <- read.table("ARPE19_uninf-2.sup.notrim.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)

#TRIM <- read.table("CRO-AP5-1.sup.trimAdapters.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)
#NOTRIM <- read.table("CRO-AP5-1.sup.notrim.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)

#TRIM <- read.table("NHDF-uninf-1.sup.trimAdapters.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)
#NOTRIM <- read.table("NHDF-uninf-1.sup.notrim.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)

#TRIM <- read.table("ARPE19_polyA_uninf-1.sup.trimAdapters.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)
#NOTRIM <- read.table("ARPE19_polyA_uninf-1.sup.notrim.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)

#TRIM <- read.table("ARPE19_siLa-1.hac.trim.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)
#NOTRIM <- read.table("ARPE19_siLa-1.hac.notrim.id.length.txt", sep="\t", header=T, stringsAsFactors=FALSE)


names(TRIM)[names(TRIM) == "sequence_length_template"] <- "seq_length_trim"
names(NOTRIM)[names(NOTRIM) == "sequence_length_template"] <- "seq_length_notrim"

merged<-merge(TRIM, NOTRIM, by="read_id")

merged$diff<-merged$seq_length_notrim - merged$seq_length_trim

filtered_merged <- merged %>%
  filter(diff <= 500 & diff >= -50)

filt_test  <- merged %>%
  filter(diff <= 10 & diff >= -10)

filt_test2  <- merged %>%
  filter(diff <= 500 & diff >= 10)


# Specify bin edges with intervals of 10
bin_width <- 10
custom_bins <- seq(0, 100, by = bin_width)

options(scipen = 100000)

# Create a histogram with custom bin sizes
pdf(file = "ARPE19-1.pdf", width = 7, height = 6)
#pdf(file = "ARPE19-2.pdf", width = 7, height = 6)
#pdf(file = "CRO-AP5-1.pdf", width = 7, height = 6)
#pdf(file = "NHDF-1.pdf", width = 7, height = 6)
#pdf(file = "ARPE19-1-poly(A).pdf", width = 7, height = 6)
#pdf(file = "ARPE19-siLa.pdf", width = 7, height = 6)

par(lwd = 1,         # Set line width for axis
    cex.axis = 1.5,  # Set size of axis labels (tick marks)
    cex.lab = 1.5,   # Set size of axis titles
    cex.main = 1.5)  # Set size of the plot title
hist(filtered_merged$diff, 
     breaks = 250, 
     xlim=c(-50,250), 
     ylim=c(0,250000),
     xlab = "trim length (estimated adapter length)", 
     main = "",
#     col = "slategray1",
     col = "cadetblue",

     lwd = 2)
dev.off()





### Determine numbers of reads falling into each trim category

untrimmed<-sum(merged$diff < 20)

trimmed<-sum(merged$diff >= 20)

underTrim<-sum(merged$diff >= 20 & merged$diff < 55)

correctTrim<-sum(merged$diff >= 55 & merged$diff < 85)

overTrim<-sum(merged$diff >= 85)

median<-median(test$diff)




