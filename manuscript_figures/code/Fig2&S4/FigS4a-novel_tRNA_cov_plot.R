library(data.table)
library(Gviz)
library(GenomicFeatures)
library(txdbmaker)

setwd("")

### Change co-ordinates and chromosome below according to novel_tx of interest (see bottom of script for details)
#novel_tRNA_1
myChr = "chr1"
myStart = 248874200
myEnd = 248874500

### Change filenames below according to novel_tx of interest
file1 <- fread('./FigS4a-novel_tRNA_1.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
file2 <- fread('./FigS4a-novel_tRNA_1.fwd.TSS.txt', col.names = c('chromosome', 'start', 'end', 'value'))
file3 <- fread('./FigS4a-novel_tRNA_1.fwd.CPAS.txt', col.names = c('chromosome', 'start', 'end', 'value'))
Regions<-read.table("FigS4a-novel_tRNA_1.regions.txt", header=T, sep="\t")

file1<-file1[file1$start>myStart]
file1<-file1[file1$end<myEnd,]
max1<-max(file1$value)

dataTrack1 <- DataTrack(range = file1, type = "a", chromosome=myChr, genome = 'SVV', fill = "darkblue", col = "darkblue", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,max1))
dataTrack2 <- DataTrack(range = file2, type = "a", chromosome=myChr, genome = 'SVV', fill = "red", col = "red", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,max1))
dataTrack3 <- DataTrack(range = file3, type = "a", chromosome=myChr, genome = 'SVV', fill = "grey", col = "grey", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,max1))

### GENERATE OVERLAY TRACKS ###
displayPars(dataTrack1) <- list(groups = factor("sample 1"))
OVERLAY <- OverlayTrack(trackList=list(dataTrack1,dataTrack2,dataTrack3), col.axis="black", background.title = "transparent")

models<-makeTxDbFromGFF("gencode.v45.tRNAs.gtf")
rtrackmodels <- GeneRegionTrack(models, genome = "Hsapiens", chromosome = myChr, name = "Gene Model", col="black", fill="black", stacking="squish", shape="smallArrow", background.title = "transparent") #squish #dense

rtrackRegions <- GeneRegionTrack(Regions, genome = "HG38", chromosome=myChr, name = "Gene Model", col="black", fill="white", shape="box", background.title = "transparent", stacking="dense")

gtrack<-GenomeAxisTrack(col="black") ##Adds genome axis

### PLOT
plotTracks(list(OVERLAY,rtrackRegions,rtrackmodels,gtrack), from = myStart, to = myEnd, sizes=c(0.75,0.05,0.1,0.1), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
#export as 3 x 10


##novel_tRNA_2
#myChr = "chr17"
#myStart = 31550000
#myEnd = 31550300

##novel_tRNA_3
#myChr = "chr14"
#myStart = 20632950
#myEnd = 20633350

##novel_tRNA_4
#myChr = "chr16"
#myStart = 3152850
#myEnd = 3153200

##novel_tRNA_5
#myChr = "chr11"
#myStart = 68459950
#myEnd = 68460250

##novel_tRNA_6
#myChr = "chr1"
#myStart = 22255620
#myEnd = 22255800




