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
library(edgeR)

#Define main directory
dir <- '/project2/gilad/ghousman/skeletal-human-chimp/'
dir2 <- paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full")

#Load data
data <- readRDS(paste0(dir2,"/data.filterMT.merge.l-TEMP.rds"))
data@meta.data$Collection <- as.factor(data@meta.data$Collection)

#Find differentially expressed genes between humans and chimpanzees
Idents(data) <- factor(data@meta.data$species, levels=c("Chimp","Human"))

#Isolate ribosomal genes
genes <- rownames(data@assays$RNA@counts)
genes.ribo <- grep('^RP',genes,value=T)
genes.no.ribo <- genes[which(!(genes %in% genes.ribo))]

#Define EDGER Function
library(scran)
library(SingleCellExperiment)
runEDGER <- function(dataSub, genes.no.ribo) {
  
  library(edgeR)

  #count matrix
  counts <- as.matrix(GetAssayData(dataSub, assay="RNA", slot="counts"))
  counts <- counts[which(rownames(counts) %in% genes.no.ribo),] #remove ribosomal genes
  
  #metadata
  metadata <- dataSub@meta.data[,c("species","Collection","nCount_RNA","Phase")]
  metadata <- metadata[colnames(counts),]
  
  #make edgeR object
  dge <- DGEList(counts)
  meta_dge <- dge$samples[,c("lib.size","norm.factors")]
  meta_dge <- cbind(meta_dge, metadata)
  dge$samples <- meta_dge
  
  #filter genes
  keep <- filterByExpr(dge, group=dge$samples$species, min.count=1, min.total.count=15)
  table(keep)
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  
  #normalize data
  dge <- calcNormFactors(dge, method="TMM")
  summary(dge$samples$norm.factors)
  
  #calculate cellular detection rate
  cdr <- scale(colMeans(dge$count > 0))
  
  #design matrix
  dge$samples$Collection <- as.factor(as.numeric(dge$samples$Collection))
  dge$samples$Phase <- as.factor(as.numeric(dge$samples$Phase))
  dge$samples$species <- as.factor(as.character(dge$samples$species))
  design <- model.matrix(~Collection+cdr+nCount_RNA+Phase+species, data=dge$samples)
  dge <- estimateDisp(dge, design=design)
  print(paste0("common dispersion: ",dge$common.dispersion))

  #QLF model fit (quasi-likelihood ratio test) - reflects uncertainty in estimating dispersion for each gene
  fit <- glmQLFit(dge, design=design, robust=TRUE)
  head(fit$coefficients)

  #differential expression
  qlf <- glmQLFTest(fit, coef=dim(design)[2]) #detect DE genes between species
  tr <- glmTreat(fit, coef=dim(design)[2], lfc=log2(1.2)) #narrow list of DE genes to those with at least a fold-change of 1.2
  
  #assess output
  tt <- topTags(qlf, n=Inf, adjust.method="BH", p.value=1)
  print("All DE Genes")
  print(dim(tt$table))
  print(table(tt$table$FDR<0.05))
  print(summary(decideTests(qlf, adjust.method="BH", p.value=0.05)))

  tt <- topTags(tr, n=Inf, adjust.method="BH", sort.by="PValue", p.value=1)
  print("DE Genes with log2(1.2) Fold Change")
  print(dim(tt$table))
  print(table(tt$table$FDR<0.05))
  print(summary(decideTests(tr, adjust.method="BH", p.value=0.05)))

  return(tt)

}

#Find differentially expressed genes using differentiation stage at collection
edgr <- runEDGER(subset(data,subset=Cell.Type=="iPSC"), genes.no.ribo)
saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.ips.rds"))
edgr <- runEDGER(subset(data,subset=Cell.Type=="MSC"), genes.no.ribo)
saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.msc.rds"))
edgr <- runEDGER(subset(data,subset=Cell.Type=="Osteoblast"), genes.no.ribo)
saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.ost.rds"))

#Find differentially expressed genes using cluster definitions
edgr <- runEDGER(subset(data,subset=(seurat_clusters==1|seurat_clusters==4|seurat_clusters==6)), genes.no.ribo)
saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.i.rds"))
edgr <- runEDGER(subset(data,subset=seurat_clusters==0), genes.no.ribo)
saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.m.rds"))
edgr <- runEDGER(subset(data,subset=(seurat_clusters==2|seurat_clusters==3|seurat_clusters==5)), genes.no.ribo)
saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.o.rds"))
edgr <- runEDGER(subset(data,subset=seurat_clusters==1), genes.no.ribo)
saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.i1.rds"))
edgr <- runEDGER(subset(data,subset=seurat_clusters==4), genes.no.ribo)
saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.i2.rds"))
edgr <- runEDGER(subset(data,subset=seurat_clusters==6), genes.no.ribo)
saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.i3.rds"))
edgr <- runEDGER(subset(data,subset=seurat_clusters==2), genes.no.ribo)
saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.o1.rds"))
edgr <- runEDGER(subset(data,subset=seurat_clusters==3), genes.no.ribo)
saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.o2.rds"))
edgr <- runEDGER(subset(data,subset=seurat_clusters==5), genes.no.ribo)
saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.o3.rds"))

#Find differentially expressed genes using osteogenic ad hoc definitions
edgr <- runEDGER(subset(data,subset=OstAdHoc.Assign=="preosteoblast"), genes.no.ribo)
saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.pre.rds"))
edgr <- runEDGER(subset(data,subset=OstAdHoc.Assign=="osteoblast"), genes.no.ribo)
saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.obl.rds"))
edgr <- runEDGER(subset(data,subset=OstAdHoc.Assign=="embedding osteoblast"), genes.no.ribo)
saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.emb.rds"))
#edgr <- runEDGER(subset(data,subset=OstAdHoc.Assign=="mineralizing osteocyte"), genes.no.ribo)
#saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.min.rds"))
#edgr <- runEDGER(subset(data,subset=OstAdHoc.Assign=="mature osteocyte"), genes.no.ribo)
#saveRDS(edgr.mat, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.mat.rds"))

#Find differentially expressed genes using grade of membership definitions
edgr <- runEDGER(subset(data,subset=GoM.Assign=="iPSC"), genes.no.ribo)
saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.gomi.rds"))
edgr <- runEDGER(subset(data,subset=GoM.Assign=="MSC"), genes.no.ribo)
saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.gomm.rds"))
edgr <- runEDGER(subset(data,subset=GoM.Assign=="Osteoblast"), genes.no.ribo)
saveRDS(edgr, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.edger-compare.gomo.rds"))


