#Load libraries
library(Seurat)
library(dplyr)
library(stringi)
library(stringr)
library(ggplot2)
library(colorspace)
library(RColorBrewer)

#Define main directory
dir <- '/project2/gilad/ghousman/skeletal-human-chimp/'

#Load batch info
batch <- read.csv(file=paste0(dir,'human-chimp-skeletal-scRNA/data/scrna-batch.csv'), header=TRUE, sep=",")

#Read in files
data <- readRDS(paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.rds"))

#Separate human and chimp data within each collection into separate objects
objects <- list()
idx <- 1
for (i in 1:length(batch$Sample_Name_at_Core))
{
  obj <- subset(data[[i]], species=='Human')
  objects[[idx]] <- CreateSeuratObject(obj@assays$RNA@counts, meta.data=obj@meta.data)
  idx <- idx+1
  obj <- subset(data[[i]], species=='Chimp')
  objects[[idx]] <- CreateSeuratObject(obj@assays$RNA@counts, meta.data=obj@meta.data)
  idx <- idx+1
}
rm(data,obj)

#Merge all cell types in humans and all cell types in chimps
h.ipsc  <- merge(objects[[1]], y=c(objects[[7]],objects[[13]],objects[[19]],objects[[25]],objects[[31]],objects[[37]]),
               add.cell.ids=c("H1.I","H1r.I","H2.I","H3.I","H4.I","H5.I","H6.I"))
h.msc   <- merge(objects[[3]], y=c(objects[[9]],objects[[15]],objects[[21]],objects[[27]],objects[[33]],objects[[39]]),
               add.cell.ids=c("H1.M","H1r.M","H2.M","H3.M","H4.M","H5.M","H6.M"))
h.osteo <- merge(objects[[5]], y=c(objects[[11]],objects[[17]],objects[[23]],objects[[29]],objects[[35]],objects[[41]]),
               add.cell.ids=c("H1.O","H1r.O","H2.O","H3.O","H4.O","H5.O","H6.O"))
c.ipsc  <- merge(objects[[2]], y=c(objects[[8]],objects[[14]],objects[[20]],objects[[26]],objects[[32]],objects[[38]]),
               add.cell.ids=c("C1.I","C1r.I","C2.I","C3.I","C4.I","C5.I","C6.I"))
c.msc   <- merge(objects[[4]], y=c(objects[[10]],objects[[16]],objects[[22]],objects[[28]],objects[[34]],objects[[40]]),
               add.cell.ids=c("C1.M","C1r.M","C2.M","C3.M","C4.M","C5.M","C6.M"))
c.osteo <- merge(objects[[6]], y=c(objects[[12]],objects[[18]],objects[[24]],objects[[30]],objects[[36]],objects[[42]]),
               add.cell.ids=c("C1.O","C1r.O","C2.O","C3.O","C4.O","C5.O","C6.O"))

objects <- list(h.ipsc,h.msc,h.osteo,c.ipsc,c.msc,c.osteo)

#Perform SCTransform normalization on each collection
for (i in 1:length(objects)) {
  print(i)
  objects[[i]] <- SCTransform(objects[[i]], conserve.memory=TRUE)
}

#Save data
saveRDS(objects, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.spps-cell.sct.rds"))

