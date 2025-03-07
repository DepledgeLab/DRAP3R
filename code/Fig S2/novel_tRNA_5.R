library(data.table)
library(Gviz)
library(GenomicFeatures)
library(txdbmaker)

myChr = "chr11"
myStart = 68459950
myEnd = 68460250

file1 <- fread('./data/target12.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
file2 <- fread('./data/target12.rev.TSS.txt', col.names = c('chromosome', 'start', 'end', 'value'))
file3 <- fread('./data/target12.rev.CPAS.txt', col.names = c('chromosome', 'start', 'end', 'value'))

file1<-file1[file1$start>myStart]
file1<-file1[file1$end<myEnd,]
max1<-max(file1$value)

dataTrack1 <- DataTrack(range = file1, type = "a", chromosome=myChr, genome = 'HG38', fill = "darkblue", col = "darkblue", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,600))
dataTrack2 <- DataTrack(range = file2, type = "a", chromosome=myChr, genome = 'HG38', fill = "red", col = "red", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,600))
dataTrack3 <- DataTrack(range = file3, type = "a", chromosome=myChr, genome = 'HG38', fill = "grey", col = "grey", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,600))

### GENERATE OVERLAY TRACKS ###
displayPars(dataTrack1) <- list(groups = factor("sample 1"))
OVERLAY <- OverlayTrack(trackList=list(dataTrack1,dataTrack2,dataTrack3), col.axis="black", background.title = "transparent")

models<-makeTxDbFromGFF("gencode.v45.tRNAs.gtf")
rtrackmodels <- GeneRegionTrack(models, genome = "Hsapiens", chromosome = myChr, name = "Gene Model", col="black", fill="black", stacking="squish", shape="smallArrow", background.title = "transparent") #squish #dense

Regions<-read.table("data/novel_tRNA_5.regions.txt", header=T, sep="\t")
rtrackRegions <- GeneRegionTrack(Regions, genome = "HG38", chromosome=myChr, name = "Gene Model", col="black", fill="white", shape="box", background.title = "transparent", stacking="dense")

gtrack<-GenomeAxisTrack(col="black") ##Adds genome axis

### FORWARD STRAND
plotTracks(list(OVERLAY,gtrack,rtrackRegions), from = myStart, to = myEnd, sizes=c(0.8,0.1,0.05), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)

#export as 5 x 20
