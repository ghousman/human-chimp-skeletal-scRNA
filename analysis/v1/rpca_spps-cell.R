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
objects <- readRDS(paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.spps-cell.sct.rds"))

#Define list subset
#Select features for downstream integration (keeping at 1000 for now)
#Calculate all necessary Pearson residuals
#Identify anchors (used references + RPCA reduction method) #ALL OTHER METHODS CRASH
#Integrate datasets
#Reduce dimensionality
#Save data

#combine ipscs from all samples
obj.ipsc  <- list(objects[[1]],objects[[4]])
obj.ipsc.features <- SelectIntegrationFeatures(object.list=obj.ipsc, nfeatures=1000)
obj.ipsc <- PrepSCTIntegration(object.list=obj.ipsc, anchor.features=obj.ipsc.features)
reference_dataset <- 1 #make humans the reference
obj.ipsc <- lapply(X=obj.ipsc, FUN=RunPCA, verbose=FALSE, features=obj.ipsc.features)
obj.ipsc.anchors <- FindIntegrationAnchors(object.list=obj.ipsc, normalization.method="SCT", anchor.features=obj.ipsc.features,
                                           reference=reference_dataset, reduction="rpca")
int.ipsc <- IntegrateData(anchorset=obj.ipsc.anchors, normalization.method="SCT")
int.ipsc <- RunPCA(object=int.ipsc, npcs=100, verbose=FALSE)
pva <- int.ipsc@reductions$pca@stdev^2/int.ipsc@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
print(ndim)
int.ipsc <- RunUMAP(int.ipsc, dims=1:100)
saveRDS(int.ipsc, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.spps-cell.sct.int.ipsc.rds"))
rm(obj.ipsc,int.ipsc)

#combine mscs from all samples
obj.msc   <- list(objects[[2]],objects[[5]])
obj.msc.features <- SelectIntegrationFeatures(object.list=obj.msc, nfeatures=1000)
obj.msc <- PrepSCTIntegration(object.list=obj.msc, anchor.features=obj.msc.features)
reference_dataset <- 1 #make humans the reference
obj.msc <- lapply(X=obj.msc, FUN=RunPCA, verbose=FALSE, features=obj.msc.features)
obj.msc.anchors <- FindIntegrationAnchors(object.list=obj.msc, normalization.method="SCT", anchor.features=obj.msc.features,
                                           reference=reference_dataset, reduction="rpca")
int.msc <- IntegrateData(anchorset=obj.msc.anchors, normalization.method="SCT")
int.msc <- RunPCA(object=int.msc, npcs=100, verbose=FALSE)
pva <- int.msc@reductions$pca@stdev^2/int.msc@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
print(ndim)
int.msc <- RunUMAP(int.msc, dims=1:100)
saveRDS(int.msc, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.spps-cell.sct.int.msc.rds"))
rm(obj.msc,int.msc)

#combine osteoblasts from all samples
obj.osteo <- list(objects[[3]],objects[[6]])
obj.osteo.features <- SelectIntegrationFeatures(object.list=obj.osteo, nfeatures=1000)
obj.osteo <- PrepSCTIntegration(object.list=obj.osteo, anchor.features=obj.osteo.features)
reference_dataset <- 1 #make humans the reference
obj.osteo <- lapply(X=obj.osteo, FUN=RunPCA, verbose=FALSE, features=obj.osteo.features)
obj.osteo.anchors <- FindIntegrationAnchors(object.list=obj.osteo, normalization.method="SCT", anchor.features=obj.osteo.features,
                                           reference=reference_dataset, reduction="rpca")
int.osteo <- IntegrateData(anchorset=obj.osteo.anchors, normalization.method="SCT")
int.osteo <- RunPCA(object=int.osteo, npcs=100, verbose=FALSE)
pva <- int.osteo@reductions$pca@stdev^2/int.osteo@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
print(ndim)
int.osteo <- RunUMAP(int.osteo, dims=1:100)
saveRDS(int.osteo, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.spps-cell.sct.int.osteo.rds"))
rm(obj.osteo,int.osteo)

