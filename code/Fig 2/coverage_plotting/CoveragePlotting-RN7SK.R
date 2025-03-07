
library(data.table)
library(Gviz)
library(GenomicFeatures)

#setwd("C:/Users/depledgd/Dropbox/EBV/")
setwd("D:/Dropbox/DRAP3R/Fig 2 - mods-new/")

myChr = "chr6"
myStart = 52995500 #
myEnd = 52996100 

file1 <- fread('./bedgraphs/CRO-AP5-1.sup-m6A_bwa_W13k6T20.SR.primary-merged.RN7SK.sorted.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
file3 <- fread('./bedgraphs/RN7SK-IVT-1.sup-m6A_bwa_W13k6T20.SR.RN7SK.sorted.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

file1<-file1[file1$start>myStart]
file1<-file1[file1$end<myEnd,]
max1<-max(file1$value)

file3<-file3[file3$start>myStart]
file3<-file3[file3$end<myEnd,]
max3<-max(file3$value)


dataTrack1 <- DataTrack(range = file1, type = "a", chromosome=myChr, genome = 'EBV', fill = "#f3511f", col = "#f3511f", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,8000))
dataTrack3 <- DataTrack(range = file3, type = "a", chromosome=myChr, genome = 'EBV', fill = "grey", col = "grey", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,8000))

### READ IN GENE MODELS ###
gtrack<-GenomeAxisTrack(col="black") ##Adds genome axis

modelsPos<-makeTxDbFromGFF("gencode.v45.RN7SK.gtf")

rtrackFor <- GeneRegionTrack(modelsPos, genome = "EBV", chromosome = myChr, name = "Gene Model", col="black", fill="grey", stacking="squish", shape="smallArrow", background.title = "transparent") #squish #dense


### GENERATE PLOT ###

### BOTH STRANDS

plotTracks(list(dataTrack1,dataTrack3,rtrackFor,gtrack), from = myStart, to = myEnd, sizes=c(0.4,0.4,0.1,0.1), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)

#export landscape 10 x 7.5
