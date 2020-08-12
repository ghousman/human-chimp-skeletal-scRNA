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

#Perform SCTransform normalization on each individual and cell type
for (i in 1:length(objects)) {
  print(i)
  objects[[i]] <- SCTransform(objects[[i]], conserve.memory=TRUE)
}

#Save data
saveRDS(objects, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.indv-cell.sct.rds"))

