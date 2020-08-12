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

