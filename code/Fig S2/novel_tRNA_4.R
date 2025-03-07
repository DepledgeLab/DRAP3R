library(data.table)
library(Gviz)
library(GenomicFeatures)
library(txdbmaker)

myChr = "chr16"
myStart = 3152850
myEnd = 3153200

file1 <- fread('./data/target17.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
file2 <- fread('./data/target17.fwd.TSS.txt', col.names = c('chromosome', 'start', 'end', 'value'))
file3 <- fread('./data/target17.fwd.CPAS.txt', col.names = c('chromosome', 'start', 'end', 'value'))

file1<-file1[file1$start>myStart]
file1<-file1[file1$end<myEnd,]
max1<-max(file1$value)

dataTrack1 <- DataTrack(range = file1, type = "a", chromosome=myChr, genome = 'HG38', fill = "darkblue", col = "darkblue", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,100))
dataTrack2 <- DataTrack(range = file2, type = "a", chromosome=myChr, genome = 'HG38', fill = "red", col = "red", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,100))
dataTrack3 <- DataTrack(range = file3, type = "a", chromosome=myChr, genome = 'HG38', fill = "grey", col = "grey", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,100))

### GENERATE OVERLAY TRACKS ###
displayPars(dataTrack1) <- list(groups = factor("sample 1"))
OVERLAY <- OverlayTrack(trackList=list(dataTrack1,dataTrack2,dataTrack3), col.axis="black", background.title = "transparent")

models<-makeTxDbFromGFF("gencode.v45.tRNAs.gtf")
rtrackmodels <- GeneRegionTrack(models, genome = "Hsapiens", chromosome = myChr, name = "Gene Model", col="black", fill="black", stacking="squish", shape="smallArrow", background.title = "transparent") #squish #dense

Regions<-read.table("data/novel_tRNA_4.regions.txt", header=T, sep="\t")
rtrackRegions <- GeneRegionTrack(Regions, genome = "HG38", chromosome=myChr, name = "Gene Model", col="black", fill="white", shape="box", background.title = "transparent", stacking="dense")

gtrack<-GenomeAxisTrack(col="black") ##Adds genome axis

### FORWARD STRAND
plotTracks(list(OVERLAY,rtrackRegions,rtrackmodels,gtrack), from = myStart, to = myEnd, sizes=c(0.75,0.05,0.1,0.1), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)

#export as 5 x 20

