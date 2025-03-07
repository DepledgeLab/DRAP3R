
library(data.table)
library(Gviz)
library(GenomicFeatures)

myChr = "chr16"
myStart = 22298350
myEnd = 22298650

file1 <- fread('./data/target13.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
file2 <- fread('./data/target13.rev.TSS.txt', col.names = c('chromosome', 'start', 'end', 'value'))
file3 <- fread('./data/target13.rev.CPAS.txt', col.names = c('chromosome', 'start', 'end', 'value'))

file1<-file1[file1$start>myStart]
file1<-file1[file1$end<myEnd,]
max1<-max(file1$value)

dataTrack1 <- DataTrack(range = file1, type = "a", chromosome=myChr, genome = 'HG38', fill = "goldenrod1", col = "goldenrod1", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,400))
dataTrack2 <- DataTrack(range = file2, type = "a", chromosome=myChr, genome = 'HG38', fill = "red", col = "red", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,400))
dataTrack3 <- DataTrack(range = file3, type = "a", chromosome=myChr, genome = 'HG38', fill = "black", col = "black", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,400))

### GENERATE OVERLAY TRACKS ###
displayPars(dataTrack1) <- list(groups = factor("sample 1"))
OVERLAY <- OverlayTrack(trackList=list(dataTrack1,dataTrack2,dataTrack3), col.axis="black", background.title = "transparent")

Regions<-read.table("data/novel_tx_7.regions.txt", header=T, sep="\t")
rtrackRegions <- GeneRegionTrack(Regions, genome = "HG38", chromosome=myChr, name = "Gene Model", col="black", fill="white", shape="box", background.title = "transparent", stacking="dense")

### FORWARD STRAND
plotTracks(list(OVERLAY,rtrackRegions,gtrack), from = myStart, to = myEnd, sizes=c(0.80,0.1,0.1),cex.axis = 4, type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)


