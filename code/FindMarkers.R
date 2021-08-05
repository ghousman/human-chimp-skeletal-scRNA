#Read in command line arguments as list of character vectors
args=(commandArgs(TRUE))

#Check if arguments are passed and cycle through to evaluate each element
if(length(args)==0){
  print("No arguments supplied.")
  #supply default values
  res="0.1"
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
print(paste0("res: ",res))

#Load libraries
library(Seurat)
library(cvms)
library(broom)
library(tibble)
library(dplyr)
library(stringi)
library(stringr)

#Load data
#data integrated across individuals - total (conservative cell filter + non-zero genes + regress out UMI/mito)
scrna <- "../data/cellranger-data-full/data.filterC.log.indv-cell.intNo0.reg-tot.assign.rds"
data <- readRDS(scrna)
data@meta.data$Cell <- rownames(data@meta.data)

#Get total genes included in topics modeling
counts <- as.matrix(data@assays$RNA@counts) #do not use integrated data
counts <- counts[rowSums(counts>0)>0,]          #remove genes with 0 counts across cells
geneTot <- rownames(counts)

#Find top 100 markers for every cluster compared to all remaining clusters, report only the positive ones
if(res=="0.1") {Idents(data) <- data@meta.data$integrated_snn_res.0.1}
if(res=="0.25") {Idents(data) <- data@meta.data$integrated_snn_res.0.25}
if(res=="0.5") {Idents(data) <- data@meta.data$integrated_snn_res.0.5}
if(res=="0.75") {Idents(data) <- data@meta.data$integrated_snn_res.0.75}
if(res=="1")   {Idents(data) <- data@meta.data$integrated_snn_res.1}

markers <- FindAllMarkers(data, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25, features=geneTot)
markers %>% group_by(cluster) %>% top_n(n=100, wt=avg_logFC)
top5 <- markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_logFC)
saveRDS(markers, file=paste0("../output/markers",res,".data.filterC.log.indv-cell.intNo0.reg-tot.rds"))
