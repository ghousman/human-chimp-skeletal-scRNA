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
objects <- readRDS(paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.filterMT.indv-cell.log.rds"))

#Define list subset
#Select features for downstream integration (keeping at 10000 for now)
#Identify anchors (used references + RPCA reduction method) #ALL OTHER METHODS CRASH
#Integrate datasets
#Reduce dimensionality
#Save data

#combine ipscs from all samples
obj.ipsc  <- list(objects[[1]],objects[[2]],objects[[7]],objects[[8]],objects[[13]],objects[[14]],objects[[19]],objects[[20]],
                  objects[[25]],objects[[26]],objects[[31]],objects[[32]],objects[[37]],objects[[38]])
obj.ipsc.features <- SelectIntegrationFeatures(object.list=obj.ipsc, nfeatures=10000)
reference_dataset <- c(3:4) #H1C1-r2 differentiated the best
obj.ipsc <- lapply(X=obj.ipsc, FUN=function(x) {
    x <- ScaleData(x, features=obj.ipsc.features, verbose=FALSE)
    x <- RunPCA(x, features=obj.ipsc.features, verbose=FALSE)
})
obj.ipsc.anchors <- FindIntegrationAnchors(object.list=obj.ipsc, normalization.method="LogNormalize", anchor.features=obj.ipsc.features,
                                           reference=reference_dataset, reduction="rpca")
int.ipsc <- IntegrateData(anchorset=obj.ipsc.anchors, normalization.method="LogNormalize")
int.ipsc <- ScaleData(int.ipsc, verbose=FALSE, vars.to.regress=c("nCount_RNA","nFeature_RNA","percent.mt"))
int.ipsc <- RunPCA(object=int.ipsc, npcs=100, verbose=FALSE)
pva <- int.ipsc@reductions$pca@stdev^2/int.ipsc@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
print(ndim)
int.ipsc <- RunUMAP(int.ipsc, dims=1:ndim)
saveRDS(int.ipsc, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.filterMT.regMT.indv-cell.log.10k.int.ipsc.rds"))
rm(obj.ipsc,int.ipsc)

#combine mscs from all samples
obj.msc   <- list(objects[[3]],objects[[4]],objects[[9]],objects[[10]],objects[[15]],objects[[16]],objects[[21]],objects[[22]],
                  objects[[27]],objects[[28]],objects[[33]],objects[[34]],objects[[39]],objects[[40]])
obj.msc.features <- SelectIntegrationFeatures(object.list=obj.msc, nfeatures=10000)
reference_dataset <- c(3:4) #H1C1-r2 differentiated the best
obj.msc <- lapply(X=obj.msc, FUN=function(x) {
    x <- ScaleData(x, features=obj.msc.features, verbose=FALSE)
    x <- RunPCA(x, features=obj.msc.features, verbose=FALSE)
})
obj.msc.anchors <- FindIntegrationAnchors(object.list=obj.msc, normalization.method="LogNormalize", anchor.features=obj.msc.features,
                                           reference=reference_dataset, reduction="rpca")
int.msc <- IntegrateData(anchorset=obj.msc.anchors, normalization.method="LogNormalize")
int.msc <- ScaleData(int.msc, verbose=FALSE, vars.to.regress=c("nCount_RNA","nFeature_RNA","percent.mt"))
int.msc <- RunPCA(object=int.msc, npcs=100, verbose=FALSE)
pva <- int.msc@reductions$pca@stdev^2/int.msc@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
print(ndim)
int.msc <- RunUMAP(int.msc, dims=1:ndim)
saveRDS(int.msc, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.filterMT.regMT.indv-cell.log.10k.int.msc.rds"))
rm(obj.msc,int.msc)


#combine osteoblasts from all samples
obj.osteo <- list(objects[[5]],objects[[6]],objects[[11]],objects[[12]],objects[[17]],objects[[18]],objects[[23]],objects[[24]],
                  objects[[29]],objects[[30]],objects[[35]],objects[[36]],objects[[41]],objects[[42]])
obj.osteo.features <- SelectIntegrationFeatures(object.list=obj.osteo, nfeatures=10000)
reference_dataset <- c(3:4) #H1C1-r2 differentiated the best
obj.osteo <- lapply(X=obj.osteo, FUN=function(x) {
    x <- ScaleData(x, features=obj.osteo.features, verbose=FALSE)
    x <- RunPCA(x, features=obj.osteo.features, verbose=FALSE)
})
obj.osteo.anchors <- FindIntegrationAnchors(object.list=obj.osteo, normalization.method="LogNormalize", anchor.features=obj.osteo.features,
                                           reference=reference_dataset, reduction="rpca")
int.osteo <- IntegrateData(anchorset=obj.osteo.anchors, normalization.method="LogNormalize")
int.osteo <- ScaleData(int.osteo, verbose=FALSE, vars.to.regress=c("nCount_RNA","nFeature_RNA","percent.mt"))
int.osteo <- RunPCA(object=int.osteo, npcs=100, verbose=FALSE)
pva <- int.osteo@reductions$pca@stdev^2/int.osteo@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
print(ndim)
int.osteo <- RunUMAP(int.osteo, dims=1:ndim)
saveRDS(int.osteo, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.filterMT.regMT.indv-cell.log.10k.int.osteo.rds"))
rm(obj.osteo,int.osteo)

