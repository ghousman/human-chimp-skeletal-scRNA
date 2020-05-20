#Load libraries
library(Seurat)
library(dplyr)
library(stringi)
library(stringr)
library(ggplot2)
library(grid)
library(gbm)
library(colorspace)
library(RColorBrewer)
library(tidyr)
library(UpSetR)
library(reshape2)
library(SC3)
library(scater)
library(SingleCellExperiment)

#Define main directory
dir <- '/project2/gilad/ghousman/skeletal-human-chimp/'

#Load data
indv <- readRDS(paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.indv.log.10k.int.rds"))            #integrate across individuals - total

#Create a SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(logcounts = as.matrix(indv@assays$RNA@data)), 
                            colData = indv@meta.data)
#Define feature names in feature_symbol column
rowData(sce)$feature_symbol <- rownames(sce)
#Remove features with duplicated names
sce <- sce[!duplicated(rowData(sce)$feature_symbol), ]
#Run SC3
sce <- sc3(sce, ks=3:7, gene_filter=FALSE, biology=TRUE, n_cores=8)
#Save to a file
saveRDS(sce, file=paste0(dir,"human-chimp-skeletal-scRNA/data/cellranger-data-full/sc3.indv.log.10k.int.rds"))
