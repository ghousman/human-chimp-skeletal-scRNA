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

#Perform SCTransform normalization on each collection
for (i in 1:length(data)) {
  print(i)
  data[[i]] <- NormalizeData(data[[i]], verbose=FALSE)
  data[[i]] <- FindVariableFeatures(data[[i]], selection.method="vst", nfeatures=3000, verbose=FALSE)
}

#Save data
saveRDS(data, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.pair-cell.log.rds"))

