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
indv <- readRDS(paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.indv.log.10k.int.rds"))            #integrate across individuals - total

#Fit GoM Model
Combined_GoM <- CountClust::FitGoM(t(as.matrix(indv@assays$RNA@data)), K=c(3:7), tol=1, control=list(tmax=100))

#Save data
save(Combined_GoM, file=paste0(dir,"human-chimp-skeletal-scRNA/data/cellranger-data-full/gom3to7.indv.log.10k.int.rds"))

