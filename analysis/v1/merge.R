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

#Merge datasets
#Perform SCTransform normalization
#Reduce dimensionality
#Remove excess files
#Save merged datasets

#combine ipscs from all samples
combo.ipsc    <- merge(data[[1]], y=c(data[[4]],data[[7]],data[[10]],data[[13]],data[[16]],data[[19]]),
                      add.cell.ids=c("H1C1.I","H1C1-r2.I","H2C2.I","H3C3.I","H4C4.I","H5C5.I","H6C6.I"))
combo.ipsc <- SCTransform(combo.ipsc)
combo.ipsc <- RunPCA(object=combo.ipsc, npcs=100, verbose=FALSE)
pva <- combo.ipsc@reductions$pca@stdev^2/combo.ipsc@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
ndim #68
combo.ipsc <- RunUMAP(combo.ipsc, dims=1:68)
saveRDS(combo.ipsc, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.merge.ipsc.rds"))
rm(combo.ipsc)

#combine mscs from all samples
combo.msc    <- merge(data[[2]], y=c(data[[5]],data[[8]],data[[11]],data[[14]],data[[17]],data[[20]]),
                      add.cell.ids=c("H1C1.M","H1C1-r2.M","H2C2.M","H3C3.M","H4C4.M","H5C5.M","H6C6.M"))
combo.msc <- SCTransform(combo.msc)
combo.msc <- RunPCA(object=combo.msc, npcs=100, verbose=FALSE)
pva <- combo.msc@reductions$pca@stdev^2/combo.msc@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
ndim #68
combo.msc <- RunUMAP(combo.msc, dims=1:68)
saveRDS(combo.msc, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.merge.msc.rds"))
rm(combo.msc)

#combine osteoblasts from all samples
combo.osteo    <- merge(data[[3]], y=c(data[[6]],data[[9]],data[[12]],data[[15]],data[[18]],data[[21]]),
                      add.cell.ids=c("H1C1.O","H1C1-r2.O","H2C2.O","H3C3.O","H4C4.O","H5C5.O","H6C6.O"))
combo.osteo <- SCTransform(combo.osteo)
combo.osteo <- RunPCA(object=combo.osteo, npcs=100, verbose=FALSE)
pva <- combo.osteo@reductions$pca@stdev^2/combo.osteo@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
ndim #68
combo.osteo <- RunUMAP(combo.osteo, dims=1:68)
saveRDS(combo.osteo, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.merge.osteo.rds"))
rm(combo.osteoblasts)

#combine all collections
combo.total  <- merge(data[[1]], y=c(data[[2]],data[[3]],data[[4]],data[[5]],data[[6]],data[[7]],data[[8]],data[[9]],data[[10]],data[[11]],
                                     data[[12]],data[[13]],data[[14]],data[[15]],data[[16]],data[[17]],data[[18]],data[[19]],data[[20]],data[[21]]),
                      add.cell.ids=c("H1C1.I","H1C1.M","H1C1.O","H1C1-r2.I","H1C1-r2.M","H1C1-r2.O","H2C2.I","H2C2.M","H2C2.O","H3C3.I","H3C3.M","H3C3.O",
                                     "H4C4.I","H4C4.M","H4C4.O","H5C5.I","H5C5.M","H5C5.O","H6C6.I","H6C6.M","H6C6.O"))
combo.total <- SCTransform(combo.total, variable.features.n=1000, conserve.memory=TRUE)
combo.total <- RunPCA(object=combo.total, npcs=100, verbose=FALSE)
pva <- combo.total@reductions$pca@stdev^2/combo.total@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
ndim #91
combo.total <- RunUMAP(combo.total, dims=1:91)
saveRDS(combo.total, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.merge.rds"))
rm(combo.total)

