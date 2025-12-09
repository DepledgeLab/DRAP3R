
library(data.table)
library(Gviz)
library(GenomicFeatures)

setwd("")

myChr = "AJ507799.2"
myStart = 1 #1
myEnd = 172000 #172000

file1 <- fread('CRO-AP5-1_bwa_W13k6T20.SR.primary-merged.EBV.sorted.fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
file2 <- fread('EBER2-IVT-3.sup-pseU_bwa_W13k6T20.SR.HG38_EBV_KSHV.sorted.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

file1<-file1[file1$start>myStart]
file1<-file1[file1$end<myEnd,]
max1<-max(file1$value)

file2<-file2[file2$start>myStart]
file2<-file2[file2$end<myEnd,]
max2<-max(file2$value)

dataTrack1 <- DataTrack(range = file1, type = "a", chromosome=myChr, genome = 'KSHV', fill = "deeppink", col = "deeppink", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,max1))
dataTrack2 <- DataTrack(range = file2, type = "a", chromosome=myChr, genome = 'KSHV', fill = "deeppink", col = "deeppink", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,max2))

### READ IN GENE MODELS ###
gtrack<-GenomeAxisTrack(col="black") ##Adds genome axis
modelsPos<-makeTxDbFromGFF("../../EBV/annotation/EBV.AJ507799.forward-v1.0.gff3")
rtrackFor <- GeneRegionTrack(modelsPos, genome = "EBV", chromosome = myChr, name = "Gene Model", col="black", fill="grey", stacking="squish", shape="smallArrow", background.title = "transparent") #squish #dense

### GENERATE PLOT ###
plotTracks(list(dataTrack1,dataTrack2,rtrackFor,gtrack), from = myStart, to = myEnd, sizes=c(0.4,0.4,0.05,0.1), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)


