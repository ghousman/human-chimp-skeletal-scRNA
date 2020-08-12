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

#Define head directory containing scRNA data
dir_scrna <- paste0(dir,'scRNA/cellranger-data-full')

#Load objects
objects <- readRDS(paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.spps.sct.rds"))

#Select features for downstream integration (keeping at 1000 for now)
#Calculate all necessary Pearson residuals
#Identify anchors (used references + RPCA reduction method) #ALL OTHER METHODS CRASH
objects.features <- SelectIntegrationFeatures(object.list=objects, nfeatures=1000)
objects <- PrepSCTIntegration(object.list=objects, anchor.features=objects.features)
reference_dataset <- 1 #make humans the reference
objects <- lapply(X=objects, FUN=RunPCA, verbose=FALSE, features=objects.features)
objects.anchors <- FindIntegrationAnchors(object.list=objects, normalization.method="SCT", anchor.features=objects.features, reference=reference_dataset, reduction="rpca")

#Integrate datasets
integrate <- IntegrateData(anchorset=objects.anchors, normalization.method="SCT")

#Reduce dimensionality
integrate <- RunPCA(object=integrate, npcs=100, verbose=FALSE)
pva <- integrate@reductions$pca@stdev^2/integrate@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
print(ndim)
integrate <- RunUMAP(integrate, dims=1:100)

#Save data
saveRDS(integrate, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.spps.sct.int.rds"))
rm(objects,integrate)

