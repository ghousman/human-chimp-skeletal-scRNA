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

#Load objects
objects <- readRDS(paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.pair.log.rds"))

#Select features for downstream integration (keeping at 1000 for now)
#Identify anchors (used references + RPCA reduction method) #ALL OTHER METHODS CRASH
objects.features <- SelectIntegrationFeatures(object.list=objects, nfeatures=3000)
reference_dataset <- 2 #make H1C1r2 the references
objects <- lapply(X=objects, FUN=function(x) {
    x <- ScaleData(x, features=objects.features, verbose=FALSE)
    x <- RunPCA(x, features=objects.features, verbose=FALSE)
})
objects.anchors <- FindIntegrationAnchors(object.list=objects, normalization.method="LogNormalize", anchor.features=objects.features,
                                          reference=reference_dataset, reduction="rpca")

#Integrate datasets
integrate <- IntegrateData(anchorset=objects.anchors, normalization.method="LogNormalize")
integrate <- ScaleData(integrate, verbose=FALSE)

#Reduce dimensionality
integrate <- RunPCA(object=integrate, npcs=100, verbose=FALSE)
pva <- integrate@reductions$pca@stdev^2/integrate@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
print(ndim)
integrate <- RunUMAP(integrate, dims=1:91)

#Save data
saveRDS(integrate, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.pair.log.3k.int.rds"))
rm(objects,integrate)

