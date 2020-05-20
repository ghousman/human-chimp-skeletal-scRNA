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

#Define MAST Function
library(MAST)
runMAST <- function(dataSub, genes.no.ribo) {
  
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
  dge <- edgeR::calcNormFactors(dge, method="TMM")
  summary(dge$samples$norm.factors)
  cpms <- edgeR::cpm(dge)
  
  #calculate cellular detection rate
  cdr <- scale(colMeans(dge$count > 0))
  
  #design matrix
  dge$samples$Collection <- as.factor(as.numeric(dge$samples$Collection))
  dge$samples$Phase <- as.factor(as.numeric(dge$samples$Phase))
  dge$samples$species <- as.factor(as.character(dge$samples$species))
  sca <- FromMatrix(exprsArray=log2(cpms+1),
                    cData=data.frame(wellKey=rownames(dge$samples),
                                     species=dge$samples$species,
                                     Collection=dge$samples$Collection,
                                     nCount_RNA=dge$samples$nCount_RNA,
                                     Phase=dge$samples$Phase,
                                     cdr=cdr))
  
  #hurdle model fit
  zlmdata <- zlm(~Collection+cdr+nCount_RNA+Phase+species, sca)
  show(zlmdata)

  #differential expression
  mast <- lrTest(zlmdata, "species")
  
  #assess output
  print(length(mast[, "hurdle", "Pr(>Chisq)"]))
  print(table(mast[, "hurdle", "Pr(>Chisq)"]<0.05))

  return(mast)

}


#Find differentially expressed genes using differentiation stage at collection
mast <- runMAST(subset(data,subset=Cell.Type=="iPSC"), genes.no.ribo)
saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.ips.rds"))
mast <- runMAST(subset(data,subset=Cell.Type=="MSC"), genes.no.ribo)
saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.msc.rds"))
mast <- runMAST(subset(data,subset=Cell.Type=="Osteoblast"), genes.no.ribo)
saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.ost.rds"))

#Find differentially expressed genes using cluster definitions
mast <- runMAST(subset(data,subset=(seurat_clusters==1|seurat_clusters==4|seurat_clusters==6)), genes.no.ribo)
saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.i.rds"))
mast <- runMAST(subset(data,subset=seurat_clusters==0), genes.no.ribo)
saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.m.rds"))
mast <- runMAST(subset(data,subset=(seurat_clusters==2|seurat_clusters==3|seurat_clusters==5)), genes.no.ribo)
saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.o.rds"))
mast <- runMAST(subset(data,subset=seurat_clusters==1), genes.no.ribo)
saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.i1.rds"))
mast <- runMAST(subset(data,subset=seurat_clusters==4), genes.no.ribo)
saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.i2.rds"))
mast <- runMAST(subset(data,subset=seurat_clusters==6), genes.no.ribo)
saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.i3.rds"))
mast <- runMAST(subset(data,subset=seurat_clusters==2), genes.no.ribo)
saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.o1.rds"))
mast <- runMAST(subset(data,subset=seurat_clusters==3), genes.no.ribo)
saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.o2.rds"))
mast <- runMAST(subset(data,subset=seurat_clusters==5), genes.no.ribo)
saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.o3.rds"))

#Find differentially expressed genes using osteogenic ad hoc definitions
mast <- runMAST(subset(data,subset=OstAdHoc.Assign=="preosteoblast"), genes.no.ribo)
saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.pre.rds"))
mast <- runMAST(subset(data,subset=OstAdHoc.Assign=="osteoblast"), genes.no.ribo)
saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.obl.rds"))
mast <- runMAST(subset(data,subset=OstAdHoc.Assign=="embedding osteoblast"), genes.no.ribo)
saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.emb.rds"))
#mast <- runMAST(subset(data,subset=OstAdHoc.Assign=="mineralizing osteocyte"), genes.no.ribo)
#saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.min.rds"))
#mast <- runMAST(subset(data,subset=OstAdHoc.Assign=="mature osteocyte"), genes.no.ribo)
#saveRDS(mast.mat, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.mat.rds"))

#Find differentially expressed genes using grade of membership definitions
mast <- runMAST(subset(data,subset=GoM.Assign=="iPSC"), genes.no.ribo)
saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.gomi.rds"))
mast <- runMAST(subset(data,subset=GoM.Assign=="MSC"), genes.no.ribo)
saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.gomm.rds"))
mast <- runMAST(subset(data,subset=GoM.Assign=="Osteoblast"), genes.no.ribo)
saveRDS(mast, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.MAST-compare.gomo.rds"))


