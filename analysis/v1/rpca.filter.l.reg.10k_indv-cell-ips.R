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
objects <- readRDS(paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.filter.indv-cell.log.rds"))

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
int.ipsc <- ScaleData(int.ipsc, verbose=FALSE, vars.to.regress=c("nCount_RNA","nFeature_RNA"))
int.ipsc <- RunPCA(object=int.ipsc, npcs=100, verbose=FALSE)
pva <- int.ipsc@reductions$pca@stdev^2/int.ipsc@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
print(ndim)
int.ipsc <- RunUMAP(int.ipsc, dims=1:ndim)
saveRDS(int.ipsc, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.filter.reg.indv-cell.log.10k.int.ipsc.rds"))
rm(obj.ipsc,int.ipsc)
