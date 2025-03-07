
library(data.table)
library(Gviz)
library(GenomicFeatures)

#setwd("C:/Users/depledgd/Dropbox/EBV/")
setwd("D:/Dropbox/DRAP3R/Fig 2 - mods-new/")

myChr = "AJ507799.2"
myStart = 1 #1
myEnd = 172000 #172000

file1 <- fread('./bedgraphs/CRO-AP5-1_bwa_W13k6T20.SR.primary-merged.EBV.sorted.fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
file2 <- fread('./bedgraphs/CRO-AP5-1_bwa_W13k6T20.SR.primary-merged.EBV.sorted.rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
file3 <- fread('./bedgraphs/EBER2-IVT-3.sup-pseU_bwa_W13k6T20.SR.HG38_EBV_KSHV.sorted.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
#file4 <- fread('./bedgraphs/KSHV-iSLK-72h-1.GK18.reverse.TSS.txt', col.names = c('chromosome', 'start', 'end', 'value'))
#file5 <- fread('./support/KSHV-iSLK-72h-1.GK18.forward.CPAS.txt', col.names = c('chromosome', 'start', 'end', 'value'))
#file6 <- fread('./support/KSHV-iSLK-72h-1.GK18.reverse.CPAS.txt', col.names = c('chromosome', 'start', 'end', 'value'))


file1<-file1[file1$start>myStart]
file1<-file1[file1$end<myEnd,]
max1<-max(file1$value)

file2<-file2[file2$start>myStart]
file2<-file2[file2$end<myEnd,]
max2<-max(file2$value)

file3<-file3[file3$start>myStart]
file3<-file3[file3$end<myEnd,]
max3<-max(file3$value)

#file4<-file4[file4$start>myStart]
#file4<-file4[file4$end<myEnd,]
#max4<-max(file4$value)

#file5<-file5[file5$start>myStart]
#file5<-file5[file5$end<myEnd,]
#max5<-max(file5$value)

#file6<-file6[file6$start>myStart]
#file6<-file6[file6$end<myEnd,]
#max6<-max(file6$value)


dataTrack1 <- DataTrack(range = file1, type = "a", chromosome=myChr, genome = 'KSHV', fill = "deeppink", col = "deeppink", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,1500000))
dataTrack2 <- DataTrack(range = file2, type = "a", chromosome=myChr, genome = 'KSHV', fill = "deeppink", col = "deeppink", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(max2,0))
dataTrack3 <- DataTrack(range = file3, type = "a", chromosome=myChr, genome = 'KSHV', fill = "grey", col = "grey", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,120000))
#dataTrack4 <- DataTrack(range = file4, type = "a", chromosome=myChr, genome = 'KSHV', fill = "red", col = "red", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(max4,0))
#dataTrack5 <- DataTrack(range = file5, type = "a", chromosome=myChr, genome = 'KSHV', fill = "black", col = "black", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(0,max5))
#dataTrack6 <- DataTrack(range = file6, type = "a", chromosome=myChr, genome = 'KSHV', fill = "black", col = "black", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent", ylim=c(max6,0))


### READ IN GENE MODELS ###
gtrack<-GenomeAxisTrack(col="black") ##Adds genome axis

modelsPos<-makeTxDbFromGFF("../../EBV/annotation/EBV.AJ507799.forward-v1.0.gff3")
modelsNeg<-makeTxDbFromGFF("../../EBV/annotation/EBV.AJ507799.reverse-v1.0.gff3")

rtrackFor <- GeneRegionTrack(modelsPos, genome = "EBV", chromosome = myChr, name = "Gene Model", col="black", fill="grey", stacking="squish", shape="smallArrow", background.title = "transparent") #squish #dense
rtrackRev <- GeneRegionTrack(modelsNeg, genome = "EBV", chromosome = myChr, name = "Gene Model", col="black", fill="grey", stacking="squish", shape="smallArrow", background.title = "transparent") #squish #dense


### GENERATE PLOT ###

### BOTH STRANDS

plotTracks(list(dataTrack1,dataTrack3,rtrackFor,gtrack,rtrackRev), from = myStart, to = myEnd, sizes=c(0.4,0.4,0.1,0.1,0.1), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)

plotTracks(list(dataTrack1,dataTrack3,rtrackFor,gtrack), from = myStart, to = myEnd, sizes=c(0.4,0.4,0.05,0.1), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)


