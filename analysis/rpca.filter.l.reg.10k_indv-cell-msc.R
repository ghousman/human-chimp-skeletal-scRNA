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
int.msc <- ScaleData(int.msc, verbose=FALSE, vars.to.regress=c("nCount_RNA","nFeature_RNA"))
int.msc <- RunPCA(object=int.msc, npcs=100, verbose=FALSE)
pva <- int.msc@reductions$pca@stdev^2/int.msc@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
print(ndim)
int.msc <- RunUMAP(int.msc, dims=1:ndim)
saveRDS(int.msc, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.filter.reg.indv-cell.log.10k.int.msc.rds"))
rm(obj.msc,int.msc)
