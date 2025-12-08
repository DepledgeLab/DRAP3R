
library(data.table)
library(Gviz)
library(GenomicFeatures)

setwd("")

### Change co-ordinates and chromosome below according to novel_tx of interest (see bottom of script for details)
#novel_tx_1
myChr = "chr17"
myStart = 17960200
myEnd = 17960380

### READ IN GENE MODELS ###
gtrack<-GenomeAxisTrack(col="black") ##Adds genome axis

### Change filenames below according to novel_tx of interest
file1 <- fread('./Fig2a-novel_tx_1.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
file2 <- fread('./Fig2a-novel_tx_1.fwd.TSS.txt', col.names = c('chromosome', 'start', 'end', 'value'))
file3 <- fread('./Fig2a-novel_tx_1.fwd.CPAS.txt', col.names = c('chromosome', 'start', 'end', 'value'))
Regions<-read.table("Fig2a-novel_tx_1.regions.txt", header=T, sep="\t")

file1<-file1[file1$start>myStart]
file1<-file1[file1$end<myEnd,]
max1<-max(file1$value)

dataTrack1 <- DataTrack(range = file1, type = "a", chromosome=myChr, genome = 'HG38', fill = "#ffc000", col = "goldenrod1", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,max1))
dataTrack2 <- DataTrack(range = file2, type = "a", chromosome=myChr, genome = 'HG38', fill = "red", col = "red", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,max1))
dataTrack3 <- DataTrack(range = file3, type = "a", chromosome=myChr, genome = 'HG38', fill = "black", col = "black", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,max1))

### GENERATE OVERLAY TRACKS ###
displayPars(dataTrack1) <- list(groups = factor("sample 1"))
OVERLAY <- OverlayTrack(trackList=list(dataTrack1,dataTrack2,dataTrack3), col.axis="black", background.title = "transparent")

rtrackRegions <- GeneRegionTrack(Regions, genome = "HG38", chromosome=myChr, name = "Gene Model", col="black", fill="white", shape="box", background.title = "transparent", stacking="dense")


### FORWARD STRAND
plotTracks(list(OVERLAY,rtrackRegions,gtrack), from = myStart, to = myEnd, sizes=c(0.80,0.1,0.1),cex.axis = 4, type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)

###export 20 x 5


##novel_tx_2
#myChr = "chr1"
#myStart = 160925900
#myEnd = 160926100

##novel_tx_3
#myChr = "chr7"
#myStart = 18038500
#myEnd = 18038750

##novel_tx_4
#myChr = "chr10"
#myStart = 49387350
#myEnd = 49387550

##novel_tx_5
#myChr = "chr7"
#myStart = 28405550
#myEnd = 28405750

##novel_tx_6
#myChr = "chr8"
#myStart = 106864020
#myEnd = 106864180

##novel_tx_7
#myChr = "chr16"
#myStart = 22298350
#myEnd = 22298650

##novel_tx_8
#myChr = "chr12"
#myStart = 93243550
#myEnd = 93243750










