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

#Perform log normalization
#Reduce dimensionality
#Remove excess files
#Save datasets

#Perform log normalization on each collection
for (i in 1:length(data)) {
  print(i)
  data[[i]] <- NormalizeData(data[[i]], verbose=FALSE)
  data[[i]] <- FindVariableFeatures(data[[i]], selection.method="vst", nfeatures=3000, verbose=FALSE)
  data[[i]] <- ScaleData(data[[i]], verbose=FALSE)
  data[[i]] <- RunPCA(object=data[[i]], npcs=100, verbose=FALSE)
  pva <- data[[i]]@reductions$pca@stdev^2/data[[i]]@reductions$pca@misc$total.variance
  ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
  ndim #
  data[[i]] <- RunUMAP(data[[i]], dims=1:ndim)
}

saveRDS(data, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.separate.rds"))
