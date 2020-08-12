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

#Subset cell atlas data matrix to several cell types
hpca.sub <- hpca[,c(239:261,416:419,470:474,496,499:515,548:551,592:599,608:616,651:665,
                    grep("Endothelial.cells", colnames(hpca)),grep("Smooth.muscle.cells", colnames(hpca)),
                    grep("Epithelial.cells", colnames(hpca)),grep("Neurons", colnames(hpca)),
                    grep("Fibroblasts", colnames(hpca)),grep("B.cell", colnames(hpca)),grep("T.cells", colnames(hpca)),
                    grep("Monocyte", colnames(hpca)),grep("Macrophage", colnames(hpca)),grep("Neutrophils", colnames(hpca)))]

#Merge with cell atlas data (integrated data and renormalizing data like ComputeCorMat() does not work)
gene.id <- rownames(indv@assays$RNA@counts)
gene.ref <- rownames(hpca.sub)
common <- intersect(gene.id, gene.ref)
merged.data <- merge(GetAssayData(indv, slot="data", assay="RNA")[common,], hpca.sub[common,], by="row.names")
#merged.data <- merge(GetAssayData(indv, slot="data", assay="RNA")[common,which(indv@meta.data$Sample=="H1-I")], hpca.sub[common,], by="row.names")
rownames(merged.data) <- merged.data$Row.names
merged.data <- merged.data[,c(2:length(colnames(merged.data)))]
labels <- sapply(strsplit(colnames(merged.data),"_"), `[`, 1)

#Semi-supervised topic modeling (define some known clusters prior to fitting topic model)

#labels for "known" clusters
idx.ips <- which(labels=="iPSC")
idx.msc <- which(labels=="MSC")
idx.ost <- which(labels=="Osteoblast" | labels=="MSC-Osteoblast.D1" | labels=="MSC-Osteoblast.D3" | labels=="MSC-Osteoblast.D7" | labels=="MSC-Osteoblast")
idx.end <- which(labels=="Endothelial.cells")
idx.smo <- which(labels=="Smooth.muscle.cells")
idx.epi <- which(labels=="Epithelial.cells")
idx.nrn <- which(labels=="Neurons")
idx.fib <- which(labels=="Fibroblasts")
idx.bce <- which(labels=="B.cell")
idx.tce <- which(labels=="T.cells")
idx.mon <- which(labels=="Monocyte")
idx.mac <- which(labels=="Macrophage")
idx.neu <- which(labels=="Neutrophils")
kwn.smp <- c(idx.ips,idx.msc,idx.ost,idx.end,idx.smo,idx.epi,idx.nrn,idx.fib,idx.bce,idx.tce,idx.mon,idx.mac,idx.neu)

class.labs <- c(rep("iPSC", length(idx.ips)),
                rep("MSC", length(idx.msc)),
                rep("Osteoblast", length(idx.ost)),
                rep("Endothelial.cells", length(idx.end)),
                rep("Smooth.muscle.cells", length(idx.smo)),
                rep("Epithelial.cells", length(idx.epi)),
                rep("Neurons", length(idx.nrn)),
                rep("Fibroblasts", length(idx.fib)),
                rep("B.cell", length(idx.bce)),
                rep("T.cells", length(idx.tce)),
                rep("Monocyte", length(idx.mon)),
                rep("Macrophage", length(idx.mac)),
                rep("Neutrophils", length(idx.neu)))

#Perform topic modeling k=13:14 (maybe change tol=0.1 at some point)
#Save GoM data
for(k in c(13:14)) {
  print(k)
  tpx.clust <- classtpx::class_topics(
    t(merged.data),
    K=k,
    known_samples=kwn.smp,
    class_labs=class.labs,
    method="omega.fix",
    tol=1)
  save(tpx.clust, file=paste0(dir,"human-chimp-skeletal-scRNA/data/cellranger-data-full/gom.sup.",k,".indv.log.10k.int.rds"))
  print(paste0("Finished k = ",k))
}

