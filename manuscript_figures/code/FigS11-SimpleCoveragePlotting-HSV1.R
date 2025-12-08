
library(data.table)
library(Gviz)
library(GenomicFeatures)

setwd("")

myChr = "HSV1-Kos" # name in first column in bedgraph
myStart = 1
myEnd = 152000 #151974

file1 <- fread('ARPE19-HSV1-6hpi-1.w13k6t20.HSV1.sorted.fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
file2 <- fread('ARPE19-HSV1-6hpi-1.w13k6t20.HSV1.sorted.rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
file3 <- fread('ARPE19-HSV1-6hpi-2.w13k6t20.HSV1.sorted.fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
file4 <- fread('ARPE19-HSV1-6hpi-2.w13k6t20.HSV1.sorted.rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
file5 <- fread('ARPE19-HSV1-12hpi-1.sup-pseU.primary-merged.HSV1-Kos.sorted.fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
file6 <- fread('ARPE19-HSV1-12hpi-1.sup-pseU.primary-merged.HSV1-Kos.sorted.rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
file9 <- fread('ARPE19-HSV1-12hpi-2.sup-pseU.primary-merged.HSV1-Kos.sorted.fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
file10 <- fread('ARPE19-HSV1-12hpi-2.sup-pseU.primary-merged.HSV1-Kos.sorted.rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))


#Add more off these if more tracks are needed

dataTrack1 <- DataTrack(range = file1, type = "a", chromosome=myChr, genome = 'HSV1-Kos', fill = "#ff1493", col = "#ff1493", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,100))
dataTrack2 <- DataTrack(range = file2, type = "a", chromosome=myChr, genome = 'HSV1-Kos', fill = "#ff1493", col = "#ff1493", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(100,0))
dataTrack3 <- DataTrack(range = file3, type = "a", chromosome=myChr, genome = 'HSV1-Kos', fill = "#ff1493", col = "#ff1493", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,100))
dataTrack4 <- DataTrack(range = file4, type = "a", chromosome=myChr, genome = 'HSV1-Kos', fill = "#ff1493", col = "#ff1493", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(100,0))
dataTrack5 <- DataTrack(range = file5, type = "a", chromosome=myChr, genome = 'HSV1-Kos', fill = "#8b0a50", col = "#8b0a50", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,3000))
dataTrack6 <- DataTrack(range = file6, type = "a", chromosome=myChr, genome = 'HSV1-Kos', fill = "#8b0a50", col = "#8b0a50", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(3000,0))
dataTrack9 <- DataTrack(range = file9, type = "a", chromosome=myChr, genome = 'HSV1-Kos', fill = "#8b0a50", col = "#8b0a50", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,3000))
dataTrack10 <- DataTrack(range = file10, type = "a", chromosome=myChr, genome = 'HSV1-Kos', fill = "#8b0a50", col = "#8b0a50", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(3000,0))

### READ IN GENE MODELS ###
gtrack<-GenomeAxisTrack(col="black") ##Adds genome axis

modelsPos<-makeTxDbFromGFF("HSV1-Kos-forward-v1.2.gff3")
modelsNeg<-makeTxDbFromGFF("HSV1-Kos-reverse-v1.2.gff3")

rtrackFor <- GeneRegionTrack(modelsPos, genome = "HSV1-Kos", chromosome = myChr, name = "Gene Model", col="black", fill="grey", stacking="squish", shape="smallArrow", background.title = "transparent") #squish #dense
rtrackRev <- GeneRegionTrack(modelsNeg, genome = "HSV1-Kos", chromosome = myChr, name = "Gene Model", col="black", fill="grey", stacking="squish", shape="smallArrow", background.title = "transparent") #squish #dense

### GENERATE PLOT ###

plotTracks(list(dataTrack1,dataTrack3,dataTrack5,dataTrack9,rtrackFor,gtrack,rtrackRev,dataTrack10,dataTrack6,dataTrack4,dataTrack2), from = myStart, to = myEnd, sizes=c(0.1,0.1,0.1,0.1,0.06,0.1,0.06,0.1,0.1,0.1,0.1), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)



