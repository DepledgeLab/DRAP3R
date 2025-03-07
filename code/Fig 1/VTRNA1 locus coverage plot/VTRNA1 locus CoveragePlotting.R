library(data.table)
library(Gviz)
library(GenomicFeatures)

myChr = "chr5"

myStart = 140710000
myEnd = 140730000

bedFile1 <- fread('./ARPE19_UNINF_24h-4_dorado.0.6.0_bwa_W13k6T20.merged.VTRNA1.sorted.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = 'Hsapiens', name = "Seq. Depth", fill = "#f3511f", col = "na", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,50000))

### READ IN GENE MODELS ###
gtrack<-GenomeAxisTrack(col="black")
models<-makeTxDbFromGFF("gencode.v45.PolIIItx.gtf")
rtrackmodels <- GeneRegionTrack(models, genome = "Hsapiens", chromosome = "chr5", name = "Gene Model", col="black", fill="black", stacking="squish", shape="smallArrow", background.title = "transparent") #squish #dense

### GENERATE PLOT ###

#export as 15 x 3
plotTracks(list(dataTrack1,rtrackmodels,gtrack), from = myStart, to = myEnd, sizes=c(0.85,0.05,0.1), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)

#export as 5 x 3
VT1Start = 140711250
VT1End = 140711400
plotTracks(list(dataTrack1,rtrackmodels,gtrack), from = VT1Start, to = VT1End, sizes=c(0.85,0.05,0.1), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)

VT2Start = 140718900
VT2End = 140719060
dataTrack2 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = 'Hsapiens', name = "Seq. Depth", fill = "#f3511f", col = "na", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,10000))
plotTracks(list(dataTrack2,rtrackmodels,gtrack), from = VT2Start, to = VT2End, sizes=c(0.85,0.05,0.1), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)

VT3Start = 140726100 
VT3End = 140726400
dataTrack3 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = 'Hsapiens', name = "Seq. Depth", fill = "#f3511f", col = "na", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,3000))
plotTracks(list(dataTrack3,rtrackmodels,gtrack), from = VT3Start, to = VT3End, sizes=c(0.85,0.05,0.1), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
