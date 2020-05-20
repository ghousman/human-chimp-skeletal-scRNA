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
dir2 <- paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full")

#Load data
#scrna <- "/data.filterMT.regMT.indv.log.10k.int.rds" #integrate across individuals - total
scrna <- "/data.filterMT.merge.l-TEMP.rds"          #merge all data without integration - total
indv <- readRDS(paste0(dir2,scrna))

#Find differentially expressed genes between humans and chimpanzees
Idents(indv) <- factor(indv@meta.data$species, levels=c("Chimp","Human"))

#Find differentially expressed genes using differentiation stage at collection
compare.ips  <- FindMarkers(subset(indv,subset=Cell.Type=="iPSC"),
                            assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.25,
                            test.use="wilcox")
saveRDS(compare.ips, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.ips.rds"))
compare.msc  <- FindMarkers(subset(indv,subset=Cell.Type=="MSC"),
                            assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.25,
                            test.use="wilcox")
saveRDS(compare.msc, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.msc.rds"))
compare.ost  <- FindMarkers(subset(indv,subset=Cell.Type=="Osteoblast"),
                            assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.25,
                            test.use="wilcox")
saveRDS(compare.ost, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.ost.rds"))
rm(compare.ips,compare.msc,compare.ost)

#Find differentially expressed genes using cluster definitions
compare.i  <- FindMarkers(subset(indv,subset=(seurat_clusters==1|seurat_clusters==4|seurat_clusters==6)),
                          assay="RNA", slot="data",
                          ident.1="Human", ident.2="Chimp",
                          logfc.threshold=0.25, min.pct=0.25,
                          test.use="wilcox")
saveRDS(compare.i, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.i.rds"))
compare.m  <- FindMarkers(subset(indv,subset=seurat_clusters==0),
                          assay="RNA", slot="data",
                          ident.1="Human", ident.2="Chimp",
                          logfc.threshold=0.25, min.pct=0.25,
                          test.use="wilcox")
saveRDS(compare.m, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.m.rds"))
compare.o  <- FindMarkers(subset(indv,subset=(seurat_clusters==2|seurat_clusters==3|seurat_clusters==5)),
                          assay="RNA", slot="data",
                          ident.1="Human", ident.2="Chimp",
                          logfc.threshold=0.25, min.pct=0.25,
                          test.use="wilcox")
saveRDS(compare.o, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.o.rds"))
rm(compare.i,compare.m,compare.o)
compare.i1  <- FindMarkers(subset(indv,subset=seurat_clusters==1),
                           assay="RNA", slot="data",
                           ident.1="Human", ident.2="Chimp",
                           logfc.threshold=0.25, min.pct=0.25,
                           test.use="wilcox")
saveRDS(compare.i1, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.i1.rds"))
compare.i2  <- FindMarkers(subset(indv,subset=seurat_clusters==4),
                           assay="RNA", slot="data",
                           ident.1="Human", ident.2="Chimp",
                           logfc.threshold=0.25, min.pct=0.25,
                           test.use="wilcox")
saveRDS(compare.i2, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.i2.rds"))
compare.i3  <- FindMarkers(subset(indv,subset=seurat_clusters==6),
                           assay="RNA", slot="data",
                           ident.1="Human", ident.2="Chimp",
                           logfc.threshold=0.25, min.pct=0.25,
                           test.use="wilcox")
saveRDS(compare.i3, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.i3.rds"))
rm(compare.i1,compare.i2,compare.i3)
compare.o1  <- FindMarkers(subset(indv,subset=seurat_clusters==2),
                           assay="RNA", slot="data",
                           ident.1="Human", ident.2="Chimp",
                           logfc.threshold=0.25, min.pct=0.25,
                           test.use="wilcox")
saveRDS(compare.o1, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.o1.rds"))
compare.o2  <- FindMarkers(subset(indv,subset=seurat_clusters==3),
                           assay="RNA", slot="data",
                           ident.1="Human", ident.2="Chimp",
                           logfc.threshold=0.25, min.pct=0.25,
                           test.use="wilcox")
saveRDS(compare.o2, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.o2.rds"))
compare.o3  <- FindMarkers(subset(indv,subset=seurat_clusters==5),
                           assay="RNA", slot="data",
                           ident.1="Human", ident.2="Chimp",
                           logfc.threshold=0.25, min.pct=0.25,
                           test.use="wilcox")
saveRDS(compare.o3, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.o3.rds"))
rm(compare.o1,compare.o2,compare.o3)

#Find differentially expressed genes using ad hoc definitions
compare.hoci  <- FindMarkers(subset(indv,subset=AdHoc.Assign=="iPSC"),
                             assay="RNA", slot="data",
                             ident.1="Human", ident.2="Chimp",
                             logfc.threshold=0.25, min.pct=0.25,
                             test.use="wilcox")
saveRDS(compare.hoci, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.hoci.rds"))
compare.hocm  <- FindMarkers(subset(indv,subset=AdHoc.Assign=="MSC"),
                             assay="RNA", slot="data",
                             ident.1="Human", ident.2="Chimp",
                             logfc.threshold=0.25, min.pct=0.25,
                             test.use="wilcox")
saveRDS(compare.hocm, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.hocm.rds"))
compare.hoco  <- FindMarkers(subset(indv,subset=AdHoc.Assign=="Osteoblast"),
                             assay="RNA", slot="data",
                             ident.1="Human", ident.2="Chimp",
                             logfc.threshold=0.25, min.pct=0.25,
                             test.use="wilcox")
saveRDS(compare.hoco, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.hoco.rds"))

#Find differentially expressed genes using osteogenic ad hoc definitions
compare.pre  <- FindMarkers(subset(indv,subset=OstAdHoc.Assign=="preosteoblast"),
                            assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.25,
                            test.use="wilcox")
saveRDS(compare.pre, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.pre.rds"))
compare.obl  <- FindMarkers(subset(indv,subset=OstAdHoc.Assign=="osteoblast"),
                            assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.25,
                            test.use="wilcox")
saveRDS(compare.obl, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.obl.rds"))
compare.emb  <- FindMarkers(subset(indv,subset=OstAdHoc.Assign=="embedding osteoblast"),
                            assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.25,
                            test.use="wilcox")
saveRDS(compare.emb, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.emb.rds"))
#compare.min  <- FindMarkers(subset(indv,subset=OstAdHoc.Assign=="mineralizing osteocyte"),
#                            assay="RNA", slot="data",
#                            ident.1="Human", ident.2="Chimp",
#                            logfc.threshold=0.25, min.pct=0.25,
#                            test.use="wilcox")
#saveRDS(compare.min, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.min.rds"))
#compare.mat  <- FindMarkers(subset(indv,subset=OstAdHoc.Assign=="mature osteocyte"),
#                            assay="RNA", slot="data",
#                            ident.1="Human", ident.2="Chimp",
#                            logfc.threshold=0.25, min.pct=0.25,
#                            test.use="wilcox")
#saveRDS(compare.mat, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.mat.rds"))
rm(compare.pre,compare.obl,compare.emb)
#rm(compare.pre,compare.obl,compare.emb,compare.min,compare.mat)

#Find differentially expressed genes using grade of membership definitions
compare.gomi  <- FindMarkers(subset(indv,subset=GoM.Assign=="iPSC"),
                             assay="RNA", slot="data",
                             ident.1="Human", ident.2="Chimp",
                             logfc.threshold=0.25, min.pct=0.25,
                             test.use="wilcox")
saveRDS(compare.gomi, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.gomi.rds"))
compare.gomm  <- FindMarkers(subset(indv,subset=GoM.Assign=="MSC"),
                             assay="RNA", slot="data",
                             ident.1="Human", ident.2="Chimp",
                             logfc.threshold=0.25, min.pct=0.25,
                             test.use="wilcox")
saveRDS(compare.gomm, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.gomm.rds"))
compare.gomo  <- FindMarkers(subset(indv,subset=GoM.Assign=="Osteoblast"),
                             assay="RNA", slot="data",
                             ident.1="Human", ident.2="Chimp",
                             logfc.threshold=0.25, min.pct=0.25,
                             test.use="wilcox")
saveRDS(compare.gomo, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de.w-compare.gomo.rds"))
rm(compare.gomi,compare.gomm,compare.gomo)
