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
scrna <- "/data.filterMT.regMT.indv.log.10k.int.rds"            #integrate across individuals - total
indv <- readRDS(paste0(dir2,scrna))

#Find differentially expressed genes between humans and chimpanzees
Idents(indv) <- factor(indv@meta.data$species, levels=c("Chimp","Human"))

#Find differentially expressed genes using differentiation stage at collection
compare.ips  <- FindMarkers(subset(indv,subset=Cell.Type=="iPSC"), assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.1,
                            test.use="negbinom")
saveRDS(compare.ips, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de-ips0.rds"))
compare.msc  <- FindMarkers(subset(indv,subset=Cell.Type=="MSC"), assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.1,
                            test.use="negbinom")
saveRDS(compare.msc, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de-msc0.rds"))
compare.ost  <- FindMarkers(subset(indv,subset=Cell.Type=="Osteoblast"), assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.1,
                            test.use="negbinom")
saveRDS(compare.ost, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de-ost0.rds"))

compare.ips  <- FindMarkers(subset(indv,subset=Cell.Type=="iPSC"), assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.1,
                            test.use="negbinom", latent.vars="Sample")
saveRDS(compare.ips, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de-ips1.rds"))
compare.msc  <- FindMarkers(subset(indv,subset=Cell.Type=="MSC"), assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.1,
                            test.use="negbinom", latent.vars="Sample")
saveRDS(compare.msc, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de-msc1.rds"))
compare.ost  <- FindMarkers(subset(indv,subset=Cell.Type=="Osteoblast"), assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.1,
                            test.use="negbinom", latent.vars="Sample")
saveRDS(compare.ost, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de-ost1.rds"))

compare.ips  <- FindMarkers(subset(indv,subset=Cell.Type=="iPSC"), assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.1,
                            test.use="negbinom", latent.vars="nCount_RNA")
saveRDS(compare.ips, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de-ips2.rds"))
compare.msc  <- FindMarkers(subset(indv,subset=Cell.Type=="MSC"), assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.1,
                            test.use="negbinom", latent.vars="nCount_RNA")
saveRDS(compare.msc, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de-msc2.rds"))
compare.ost  <- FindMarkers(subset(indv,subset=Cell.Type=="Osteoblast"), assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.1,
                            test.use="negbinom", latent.vars="nCount_RNA")
saveRDS(compare.ost, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de-ost2.rds"))

compare.ips  <- FindMarkers(subset(indv,subset=Cell.Type=="iPSC"), assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.1,
                            test.use="negbinom", latent.vars=c("Sample","nCount_RNA"))
saveRDS(compare.ips, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de-ips3.rds"))
compare.msc  <- FindMarkers(subset(indv,subset=Cell.Type=="MSC"), assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.1,
                            test.use="negbinom", latent.vars=c("Sample","nCount_RNA"))
saveRDS(compare.msc, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de-msc3.rds"))
compare.ost  <- FindMarkers(subset(indv,subset=Cell.Type=="Osteoblast"), assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.1,
                            test.use="negbinom", latent.vars=c("Sample","nCount_RNA"))
saveRDS(compare.ost, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de-ost3.rds"))

compare.ips  <- FindMarkers(subset(indv,subset=Cell.Type=="iPSC"), assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.1,
                            test.use="negbinom", latent.vars=c("Sample","nCount_RNA"))
saveRDS(compare.ips, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de-ips3.rds"))
compare.msc  <- FindMarkers(subset(indv,subset=Cell.Type=="MSC"), assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.1,
                            test.use="negbinom", latent.vars=c("Sample","nCount_RNA"))
saveRDS(compare.msc, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de-msc3.rds"))
compare.ost  <- FindMarkers(subset(indv,subset=Cell.Type=="Osteoblast"), assay="RNA", slot="data",
                            ident.1="Human", ident.2="Chimp",
                            logfc.threshold=0.25, min.pct=0.1,
                            test.use="negbinom", latent.vars=c("Sample","nCount_RNA"))
saveRDS(compare.ost, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/de-ost3.rds"))

