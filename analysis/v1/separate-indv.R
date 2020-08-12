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

#Perform log normalization
#Reduce dimensionality
#Remove excess files
#Save datasets

#Perform log normalization on each collection
for (i in 1:length(objects)) {
  print(i)
  objects[[i]] <- NormalizeData(objects[[i]], verbose=FALSE)
  objects[[i]] <- FindVariableFeatures(objects[[i]], selection.method="vst", nfeatures=3000, verbose=FALSE)
  objects[[i]] <- ScaleData(objects[[i]], verbose=FALSE)
  objects[[i]] <- RunPCA(object=objects[[i]], npcs=100, verbose=FALSE)
  pva <- objects[[i]]@reductions$pca@stdev^2/objects[[i]]@reductions$pca@misc$total.variance
  ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
  ndim #
  objects[[i]] <- RunUMAP(objects[[i]], dims=1:ndim)
}

saveRDS(objects, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.separate-indv.rds"))

