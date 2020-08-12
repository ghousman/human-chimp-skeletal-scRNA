#Load libraries
library(Seurat)
library(dplyr)
library(stringi)
library(stringr)
library(ggplot2)
library(colorspace)
library(RColorBrewer)
library(tidyr)

#Compute the correlation matrix use a seurat object and cell line matrix as input
ComputeCorMat <- function(object, cell.line, CorMethod="pearson"){

    data=GetAssayData(object, slot="data", assay="RNA") #additions - added slot and assay info because SCT assay 
    data_norm <- apply(data, 2, function(x) x/sum(x))
    data_libsize <- apply(data_norm, 2, function(x) x * 1e+06)

    Data.use <- data_libsize

    gene.id <- rownames(Data.use)
    gene.ref <- rownames(cell.line)
    common <- intersect(gene.id, gene.ref)
    logxx <- apply(Data.use[common, ], 2, function(x) {
        log(x + 0.1)
    })
    selected.cell.line <- apply(cell.line[common, ], 2, function(x) {
        x - mean(x)
    })
    n <- ncol(logxx)
    m <- ncol(cell.line)
    cor.mat <- matrix(0, nrow = n, ncol = m)
    for (j in 1:m) {
        for (i in 1:n) {
            cor.mat[i,j] <- cor(logxx[,i], selected.cell.line[,j], method=CorMethod)
        }
    }
    rownames(cor.mat) <- colnames(Data.use)
    colnames(cor.mat) <- colnames(cell.line)
	return(cor.mat)
}

#Compute the correlation matrix use a seurat object and cell line matrix as input
hpca <- read.table("/project2/gilad/ghousman/skeletal-human-chimp/human-chimp-skeletal-scRNA/data/HumanPrimaryCellAtlasData-SingleR")
hpca.names <- read.csv("/project2/gilad/ghousman/skeletal-human-chimp/human-chimp-skeletal-scRNA/data/HumanPrimaryCellAtlasLabels-SingleR.csv", row.names=1)
hpca.names$label <- paste0(hpca.names$geo.currated,"_",hpca.names$label.currated)
i=1
while(i <= length(colnames(hpca))) {
  colnames(hpca)[i] <- hpca.names$label[rownames(hpca.names)==colnames(hpca)[i]]
  i=i+1
}

#Define directory
dir <- '/project2/gilad/ghousman/skeletal-human-chimp/'

#Compute correlation matrix for each data set
indv <- readRDS(paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.indv-cell.log.10k.int.ipsc.rds"))  #integrate across individuals - ipsc
cor.mat.hpca=ComputeCorMat(indv,hpca)
write.table(cor.mat.hpca, file="/project2/gilad/ghousman/skeletal-human-chimp/human-chimp-skeletal-scRNA/data/corrmatrix.ips")
rm(indv,cor.mat.hpca)

indv <- readRDS(paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.indv-cell.log.10k.int.msc.rds"))   #integrate across individuals - msc
cor.mat.hpca=ComputeCorMat(indv,hpca)
write.table(cor.mat.hpca, file="/project2/gilad/ghousman/skeletal-human-chimp/human-chimp-skeletal-scRNA/data/corrmatrix.msc")
rm(indv,cor.mat.hpca)

indv <- readRDS(paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.indv-cell.log.10k.int.osteo.rds")) #integrate across individuals - osteo
cor.mat.hpca=ComputeCorMat(indv,hpca)
write.table(cor.mat.hpca, file="/project2/gilad/ghousman/skeletal-human-chimp/human-chimp-skeletal-scRNA/data/corrmatrix.ost")
rm(indv,cor.mat.hpca)

indv <- readRDS(paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.indv.log.10k.int.rds"))            #integrate across individuals - total
cor.mat.hpca=ComputeCorMat(indv,hpca)
write.table(cor.mat.hpca, file="/project2/gilad/ghousman/skeletal-human-chimp/human-chimp-skeletal-scRNA/data/corrmatrix.tot")
rm(indv,cor.mat.hpca)
