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
data <- readRDS(paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.filterMT.rds"))

#Separate human and chimp data within each collection into separate objects
objects <- list()
idx <- 1
for (i in 1:length(batch$Sample_Name_at_Core))
{
  obj <- subset(data[[i]], species=='Human')
  objects[[idx]] <- CreateSeuratObject(obj@assays$RNA@counts, meta.data=obj@meta.data)
  idx <- idx+1
  obj <- subset(data[[i]], species=='Chimp')
  objects[[idx]] <- CreateSeuratObject(obj@assays$RNA@counts, meta.data=obj@meta.data)
  idx <- idx+1
}
rm(data,obj)

#Merge datasets
#Perform log normalization
#Reduce dimensionality
#Remove excess files
#Save merged datasets

#combine ipscs from all samples
combo.ipsc    <- merge(objects[[1]], y=c(objects[[7]],objects[[13]],objects[[19]],objects[[25]],objects[[31]],
                                         objects[[37]],objects[[2]],objects[[8]],objects[[14]],objects[[20]],
                                         objects[[26]],objects[[32]],objects[[38]]),
                       add.cell.ids=c("H1.I","H1r2.I","H2.I","H3.I","H4.I","H5.I","H6.I",
				      "C1.I","C1r2.I","C2.I","C3.I","C4.I","C5.I","C6.I"))
combo.ipsc <- NormalizeData(combo.ipsc, verbose=FALSE)
combo.ipsc <- FindVariableFeatures(combo.ipsc, selection.method="vst", nfeatures=3000, verbose=FALSE)
combo.ipsc <- ScaleData(combo.ipsc, verbose=FALSE)
combo.ipsc <- RunPCA(object=combo.ipsc, npcs=100, verbose=FALSE)
pva <- combo.ipsc@reductions$pca@stdev^2/combo.ipsc@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
ndim
combo.ipsc <- RunUMAP(combo.ipsc, dims=1:ndim)
saveRDS(combo.ipsc, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.filterMT.merge.l.ipsc.rds"))
rm(combo.ipsc)

#combine mscs from all samples
combo.msc    <- merge(objects[[3]], y=c(objects[[9]],objects[[15]],objects[[21]],objects[[27]],objects[[33]],
                                        objects[[39]],objects[[4]],objects[[10]],objects[[16]],objects[[22]],
                                        objects[[28]],objects[[34]],objects[[40]]),
                       add.cell.ids=c("H1.M","H1r2.M","H2.M","H3.M","H4.M","H5.M","H6.M",
				      "C1.M","C1r2.M","C2.M","C3.M","C4.M","C5.M","C6.M"))
combo.msc <- NormalizeData(combo.msc, verbose=FALSE)
combo.msc <- FindVariableFeatures(combo.msc, selection.method="vst", nfeatures=3000, verbose=FALSE)
combo.msc <- ScaleData(combo.msc, verbose=FALSE)
combo.msc <- RunPCA(object=combo.msc, npcs=100, verbose=FALSE)
pva <- combo.msc@reductions$pca@stdev^2/combo.msc@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
ndim
combo.msc <- RunUMAP(combo.msc, dims=1:ndim)
saveRDS(combo.msc, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.filterMT.merge.l.msc.rds"))
rm(combo.msc)

#combine osteoblasts from all samples
combo.osteo    <- merge(objects[[5]], y=c(objects[[11]],objects[[17]],objects[[23]],objects[[29]],objects[[35]],
                                          objects[[41]],objects[[6]],objects[[12]],objects[[18]],objects[[24]],
                                          objects[[30]],objects[[36]],objects[[42]]),
                       add.cell.ids=c("H1.O","H1r2.O","H2.O","H3.O","H4.O","H5.O","H6.O",
			              "C1.O","C1r2.O","C2.O","C3.O","C4.O","C5.O","C6.O"))
combo.osteo <- NormalizeData(combo.osteo, verbose=FALSE)
combo.osteo <- FindVariableFeatures(combo.osteo, selection.method="vst", nfeatures=3000, verbose=FALSE)
combo.osteo <- ScaleData(combo.osteo, verbose=FALSE)
combo.osteo <- RunPCA(object=combo.osteo, npcs=100, verbose=FALSE)
pva <- combo.osteo@reductions$pca@stdev^2/combo.osteo@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
ndim
combo.osteo <- RunUMAP(combo.osteo, dims=1:ndim)
saveRDS(combo.osteo, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.filterMT.merge.l.osteo.rds"))
rm(combo.osteoblasts)

#combine all collections
combo.total  <- merge(objects[[1]], y=c(objects[[3]],objects[[5]],
		                        objects[[7]],objects[[9]],objects[[11]],
		                        objects[[13]],objects[[15]],objects[[17]],
		                        objects[[19]],objects[[21]],objects[[23]],
		                        objects[[25]],objects[[27]],objects[[29]],
		                        objects[[31]],objects[[33]],objects[[35]],
		                        objects[[37]],objects[[39]],objects[[41]],
		                        objects[[2]],objects[[4]],objects[[6]],
		                        objects[[8]],objects[[10]],objects[[12]],
		                        objects[[14]],objects[[16]],objects[[18]],
		                        objects[[20]],objects[[22]],objects[[24]],
		                        objects[[26]],objects[[28]],objects[[30]],
			    	        objects[[32]],objects[[34]],objects[[36]],
				        objects[[38]],objects[[40]],objects[[42]]),
                      add.cell.ids=c("H1.I","H1.M","H1.O","H1r2.I","H1r2.M","H1r2.O","H2.I","H2.M","H2.O",
				     "H3.I","H3.M","H3.O","H4.I","H4.M","H4.O","H5.I","H5.M","H5.O","H6.I","H6.M","H6.O",
				     "C1.I","C1.M","C1.O","C1r2.I","C1r2.M","C1r2.O","C2.I","C2.M","C2.O",
				     "C3.I","C3.M","C3.O","C4.I","C4.M","C4.O","C5.I","C5.M","C5.O","C6.I","C6.M","C6.O"))
combo.total <- NormalizeData(combo.total, verbose=FALSE)
combo.total <- FindVariableFeatures(combo.total, selection.method="vst", nfeatures=3000, verbose=FALSE)
combo.total <- ScaleData(combo.total, verbose=FALSE)
combo.total <- RunPCA(object=combo.total, npcs=100, verbose=FALSE)
pva <- combo.total@reductions$pca@stdev^2/combo.total@reductions$pca@misc$total.variance
ndim <- length(which(pva>=0.001)) #keep all dims that explaim more than 0.1% of variance
ndim
combo.total <- RunUMAP(combo.total, dims=1:ndim)
saveRDS(combo.total, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.filterMT.merge.l.rds"))
rm(combo.total)
