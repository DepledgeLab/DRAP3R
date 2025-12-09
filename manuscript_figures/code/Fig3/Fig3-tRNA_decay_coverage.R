library(data.table)
library(Gviz)
library(GenomicFeatures)
library(txdbmaker)

setwd("")

#chr17:38869000-38869500 tRNA-Cys-GCA-4-1
#chr17:75035000-75035500 tRNA-Arg-TCG-3-1
#chr2:27050700-27050950 tRNA-Tyr-GTA-2-1
#chr17:8119000-8119500  tRNA-Lys-TTT-3-5
#chr17:19508180-19508252 tRNA-Trp-CCA-2-1
#chr1:120844261-120844335 tRNA-Asn-GTT-7-1 - no decay
#chr9:35657750-35658019 RMRP - no decay

myChr = "chr9"
myStart = 35657720
myEnd = 35658050 

file1 <- fread('ARPE19_uninf-1_select_tRNAs.sorted.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
file2 <- fread('ARPE19_siLa-1_select_tRNAs.sorted.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Use for tRNA plots
#dataTrack1 <- DataTrack(range = file1, type = "hist", chromosome=myChr, genome = 'SVV', fill = "#002060", col = "#002060", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent")#, ylim=c(0,300))
#dataTrack2 <- DataTrack(range = file2, type = "hist", chromosome=myChr, genome = 'SVV', fill = "#597ec8", col = "#597ec8", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent")#, ylim=c(0,300))

### Use for Pol III plots
dataTrack1 <- DataTrack(range = file1, type = "hist", chromosome=myChr, genome = 'SVV', fill = "#f35120", col = "#f35120", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent")#, ylim=c(0,300))
dataTrack2 <- DataTrack(range = file2, type = "hist", chromosome=myChr, genome = 'SVV', fill = "#f79e84", col = "#f79e84", options(ucscChromosomeNames=FALSE),col.axis="black", background.title = "transparent")#, ylim=c(0,300))

### GENERATE OVERLAY TRACKS ###
displayPars(dataTrack1) <- list(groups = factor("sample 1"))
OVERLAY <- OverlayTrack(trackList=list(dataTrack1,dataTrack2), col.axis="black", background.title = "transparent", type="l")

### Choose appropriate model
#models<-makeTxDbFromGFF("gencode.v45.tRNAs.gtf")
models<-makeTxDbFromGFF("gencode.v45.PolIIItx.gtf")
rtrackmodels <- GeneRegionTrack(models, genome = "Hsapiens", chromosome = myChr, name = "Gene Model", col="black", fill="black", stacking="squish", shape="smallArrow", background.title = "transparent") #squish #dense

gtrack<-GenomeAxisTrack(col="black") ##Adds genome axis

### Plotting
plotTracks(list(OVERLAY,rtrackmodels,gtrack), from = myStart, to = myEnd, sizes=c(0.8,0.04,0.016), lwd = 3, type="l", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)

# export 3 x 8

