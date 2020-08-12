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
library(CountClust)
library(classtpx)

#Define main directory
dir <- '/project2/gilad/ghousman/skeletal-human-chimp/'

#Load data
indv <- readRDS(paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.indv.log.10k.int.rds")) #integrate across individuals - total

#Load cell atlas data matrix (currated by GAH)
hpca <- read.table("/project2/gilad/ghousman/skeletal-human-chimp/human-chimp-skeletal-scRNA/data/HumanPrimaryCellAtlasData-SingleR")
hpca.names <- read.csv("/project2/gilad/ghousman/skeletal-human-chimp/human-chimp-skeletal-scRNA/data/HumanPrimaryCellAtlasLabels-SingleR.csv", row.names=1)
hpca.names$label <- paste0(hpca.names$label.curated,"_",hpca.names$geo.curated)
i=1
while(i <= length(colnames(hpca))) {
  colnames(hpca)[i] <- hpca.names$label[rownames(hpca.names)==colnames(hpca)[i]]
  i=i+1
}

#Subset cell atlas data matrix
hpca.sub <- hpca[,c(239:261,416:419,470:474,496,499:515,548:551,592:599,608:616,651:665)]
#hpca.sub <- hpca[,unique(c(grep("iPSC", colnames(hpca)),grep("MSC", colnames(hpca)),grep("Osteoblast", colnames(hpca)),
#                           grep("Chondrocyte", colnames(hpca)),grep("Adipocyte", colnames(hpca))))]

#Merge with cell atlas data (integrated data and renormalizing data like ComputeCorMat() does not work)
gene.id <- rownames(indv@assays$RNA@counts)
gene.ref <- rownames(hpca.sub)
common <- intersect(gene.id, gene.ref)
merged.data <- merge(GetAssayData(indv, slot="data", assay="RNA")[common,], hpca.sub[common,], by="row.names")
#merged.data <- merge(GetAssayData(indv, slot="data", assay="RNA")[common,which(indv@meta.data$Sample=="H1-I")], hpca.sub[common,], by="row.names")
rownames(merged.data) <- merged.data$Row.names
merged.data <- merged.data[,c(2:length(colnames(merged.data)))]
labels <- sapply(strsplit(colnames(merged.data),"_"), `[`, 1)

#Fit GoM model (maybe change tol=0.1 at some point)
Combined_GoM <- CountClust::FitGoM(t(as.matrix(merged.data)), K=c(3:7), tol=1, control=list(tmax=100))

#Save data
save(Combined_GoM, file=paste0(dir,"human-chimp-skeletal-scRNA/data/cellranger-data-full/gom.ref.3to7.indv.log.10k.int.rds"))

