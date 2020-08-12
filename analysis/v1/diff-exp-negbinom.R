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

#Define SEURATnegbinom Function
runSEURATnegbinom <- function(dataSub, genes.no.ribo) {
  deGenes <- FindMarkers(dataSub,
                         assay="RNA",
                         slot="counts",
                         ident.1="Human",
                         ident.2="Chimp",
                         features=genes.no.ribo,
                         logfc.threshold=0.25,
                         min.pct=0.25,
                         min.cells.feature=3,
                         test.use="negbinom",
                         latent.vars=c("Collection","nCount_RNA","Phase"))
  return(deGenes)
}

##Find differentially expressed genes using differentiation stage at collection
#negb <- runSEURATnegbinom(subset(data,subset=Cell.Type=="iPSC"), genes.no.ribo)
#saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.ips.rds"))
#negb <- runSEURATnegbinom(subset(data,subset=Cell.Type=="MSC"), genes.no.ribo)
#saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.msc.rds"))
#negb <- runSEURATnegbinom(subset(data,subset=Cell.Type=="Osteoblast"), genes.no.ribo)
#saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.ost.rds"))

##Find differentially expressed genes using cluster definitions
#negb <- runSEURATnegbinom(subset(data,subset=(seurat_clusters==1|seurat_clusters==4|seurat_clusters==6)), genes.no.ribo)
#saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.i.rds"))
#negb <- runSEURATnegbinom(subset(data,subset=seurat_clusters==0), genes.no.ribo)
#saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.m.rds"))
#negb <- runSEURATnegbinom(subset(data,subset=(seurat_clusters==2|seurat_clusters==3|seurat_clusters==5)), genes.no.ribo)
#saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.o.rds"))
#negb <- runSEURATnegbinom(subset(data,subset=seurat_clusters==1), genes.no.ribo)
#saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.i1.rds"))
#negb <- runSEURATnegbinom(subset(data,subset=seurat_clusters==4), genes.no.ribo)
#saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.i2.rds"))
#negb <- runSEURATnegbinom(subset(data,subset=seurat_clusters==6), genes.no.ribo)
#saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.i3.rds"))
#negb <- runSEURATnegbinom(subset(data,subset=seurat_clusters==2), genes.no.ribo)
#saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.o1.rds"))
#negb <- runSEURATnegbinom(subset(data,subset=seurat_clusters==3), genes.no.ribo)
#saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.o2.rds"))
#negb <- runSEURATnegbinom(subset(data,subset=seurat_clusters==5), genes.no.ribo)
#saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.o3.rds"))

#Find differentially expressed genes using ad hoc definitions
negb <- runSEURATnegbinom(subset(data,subset=AdHoc.Assign=="iPSC"), genes.no.ribo)
saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.hoci.rds"))
negb <- runSEURATnegbinom(subset(data,subset=AdHoc.Assign=="MSC"), genes.no.ribo)
saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.hocm.rds"))
negb <- runSEURATnegbinom(subset(data,subset=AdHoc.Assign=="Osteoblast"), genes.no.ribo)
saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.hoco.rds"))

##Find differentially expressed genes using osteogenic ad hoc definitions
#negb <- runSEURATnegbinom(subset(data,subset=OstAdHoc.Assign=="preosteoblast"), genes.no.ribo)
#saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.pre.rds"))
#negb <- runSEURATnegbinom(subset(data,subset=OstAdHoc.Assign=="osteoblast"), genes.no.ribo)
#saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.obl.rds"))
#negb <- runSEURATnegbinom(subset(data,subset=OstAdHoc.Assign=="embedding osteoblast"), genes.no.ribo)
#saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.emb.rds"))
##negb <- runSEURATnegbinom(subset(data,subset=OstAdHoc.Assign=="mineralizing osteocyte"), genes.no.ribo)
##saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.min.rds"))
##negb <- runSEURATnegbinom(subset(data,subset=OstAdHoc.Assign=="mature osteocyte"), genes.no.ribo)
##saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.mat.rds"))

##Find differentially expressed genes using grade of membership definitions
#negb <- runSEURATnegbinom(subset(data,subset=GoM.Assign=="iPSC"), genes.no.ribo)
#saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.gomi.rds"))
#negb <- runSEURATnegbinom(subset(data,subset=GoM.Assign=="MSC"), genes.no.ribo)
#saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.gomm.rds"))
#negb <- runSEURATnegbinom(subset(data,subset=GoM.Assign=="Osteoblast"), genes.no.ribo)
#saveRDS(negb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.n-compare.gomo.rds"))

