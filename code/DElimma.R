#Read in command line arguments as list of character vectors
args=(commandArgs(TRUE))

#Check if arguments are passed and cycle through to evaluate each element
if(length(args)==0){
  print("No arguments supplied.")
  #supply default values
  name="noarg"
  filter.arg=FALSE
  min.count=0.05
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
print(paste0("file tag: ",name))
print(paste0("cell assignment: ",assign))
print(paste0("filter genes: ",filter.arg))
print(paste0("gene min.count: ",min.count))

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
library(scran)
library(SingleCellExperiment)

#Load data (integrate across individuals - total (15k genes))
scrna <- "./../data/cellranger-data-full/data.filter.log.indv-cell.int19k-tot.assign.rds"
data <- readRDS(scrna)
data

#Isolate mitochondrial and ribosomal genes
genes <- rownames(data@assays$RNA@counts)
genes.mito <- c("MT-ATP6","MT-ATP8","MT-CO1","MT-CO2","MT-CO3","MT-CYB","MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5")
genes.ribo <- grep('^RP',genes,value=T)
genes.no.mito.ribo <- genes[which(!(genes %in% c(genes.mito,genes.ribo)))]
rm(genes,genes.mito,genes.ribo)

#Define EDGER Function
runEDGER <- function(dataSub, cell.assign, cell.subset, genes.no.mito.ribo, filter.arg, min.count) {

  #count matrix
  counts <- as.matrix(GetAssayData(dataSub, assay="RNA", slot="counts"))

  #remove mitochodrial and ribosomal genes
  counts <- counts[which(rownames(counts) %in% genes.no.mito.ribo),]

  #metadata
  metadata <- dataSub@meta.data[,c("Species","Collection","nCount_RNA","percent.mt","Phase")]
  metadata <- metadata[colnames(counts),]

  #make edgeR object
  dge <- DGEList(counts)
  meta_dge <- dge$samples[,c("lib.size","norm.factors")]
  meta_dge <- cbind(meta_dge, metadata)
  dge$samples <- meta_dge
  rm(dataSub,counts,metadata,meta_dge)

  #filter genes
  if (filter.arg==TRUE){
    keep <- filterByExpr(dge, group=dge$samples$Species, min.count=min.count, min.total.count=15)
    table(keep)
    dge <- dge[keep, , keep.lib.sizes=FALSE]
    rm(keep)
  }

  #normalize data
  dge <- calcNormFactors(dge, method="TMM")
  summary(dge$samples$norm.factors)

  #calculate cellular detection rate
  cdr <- scale(colMeans(dge$count > 0))

  #design matrix
  dge$samples$Collection <- as.factor(as.numeric(dge$samples$Collection))
  dge$samples$Phase <- as.factor(dge$samples$Phase)
  dge$samples$Species <- as.factor(as.character(dge$samples$Species))
  model <- "~Collection+cdr+nCount_RNA+percent.mt+Phase+Species"
  design <- model.matrix(~Collection+cdr+nCount_RNA+percent.mt+Phase+Species, data=dge$samples)

  dge <- estimateDisp(dge, design=design)
  print(paste0("common dispersion: ",dge$common.dispersion))

  #QLF model fit (quasi-likelihood ratio test) - reflects uncertainty in estimating dispersion for each gene
  fit <- glmQLFit(dge, design=design, robust=TRUE)
  head(fit$coefficients)

  #differential expression between species
  qlf <- glmQLFTest(fit, coef=dim(design)[2], poisson.bound=TRUE)

  #assess output
  tt <- topTags(qlf, n=Inf, adjust.method="BH", p.value=1)
  print("All DE Genes")
  print(dim(tt$table))
  print(table(tt$table$FDR<0.01))
  print(summary(decideTests(qlf, adjust.method="BH", p.value=0.01)))

  #add details to output
  tt$table$data <- rep(scrna, dim(tt$table)[1])
  tt$table$cell.assign <- rep(cell.assign, dim(tt$table)[1])
  tt$table$cell.subset <- rep(cell.subset, dim(tt$table)[1])
  tt$table$gene.filter <- rep(filter.arg, dim(tt$table)[1])
  tt$table$gene.filter.args <- rep(paste0("min.count=",min.count), dim(tt$table)[1])
  tt$table$model <- rep(model, dim(tt$table)[1])
  tt$table$comparison <- rep(tt$comparison, dim(tt$table)[1])
  tt$table$test <- rep(tt$test, dim(tt$table)[1])
  tt$table$adjust.method <- rep(tt$adjust.method, dim(tt$table)[1])

  return(tt$table)

}

#Find differentially expressed genes for different cell assignment subtypes
if (assign=="all") {
  subset.list <- list(stage=list("stage-Time 0","stage-Time 1","stage-Time 2"),
                      clust=list("cluster-iPSC.c1","cluster-iPSC.c2","cluster-iPSC.c3","cluster-MSC.c1","cluster-Osteoblast.c1","cluster-Osteoblast.c2"),
                      combI=c("cluster-iPSC.c1_iPSC.c2_iPSC.c3"),
                      combO=c("cluster-Osteoblast.c1_Osteoblast.c2"),
                      adhoc=list("adhoc-iPSC","adhoc-MSC","adhoc-Osteoblast"),
                      ostadhoc=list("ostadhoc-preosteoblast","ostadhoc-embedding osteoblast"))
}
if (assign=="stage") {
  subset.list <- list(stage=list("stage-Time 0","stage-Time 1","stage-Time 2"))
}
if (assign=="cluster") {
  subset.list <- list(clust=list("cluster-iPSC.c1","cluster-iPSC.c2","cluster-iPSC.c3","cluster-MSC.c1","cluster-Osteoblast.c1","cluster-Osteoblast.c2"),
                      combI=c("cluster-iPSC.c1_iPSC.c2_iPSC.c3"),
                      combO=c("cluster-Osteoblast.c1_Osteoblast.c2"))
}
if (assign=="adhoc") {
  subset.list <- list(adhoc=list("adhoc-iPSC","adhoc-MSC","adhoc-Osteoblast"))
}
if (assign=="ostadhoc") {
  subset.list <- list(ostadhoc=list("ostadhoc-preosteoblast","ostadhoc-embedding osteoblast"))
}

edgr.total <- data.frame()

for (i in subset.list) {

  for (j in i) {

    print(j)

    #define cell assignment and subtype of interest
    cell.assign <- str_split(j,"-")[[1]][1]
    cell.subset <- str_split(j,"-")[[1]][2]

    #define data subset
    if (str_count(cell.subset,"_")>=1){

      if (length(str_split(cell.subset,"_")[[1]])==2) {
        dataSub <- subset(data,
                          subset=Cluster==str_split(cell.subset,"_")[[1]][1]|
                            Cluster==str_split(cell.subset,"_")[[1]][2])
      }
      else {
        dataSub <- subset(data,
                          subset=Cluster==str_split(cell.subset,"_")[[1]][1]|
                            Cluster==str_split(cell.subset,"_")[[1]][2]|
                            Cluster==str_split(cell.subset,"_")[[1]][3])
      }
    }

    else {
      if(cell.assign=="stage") { dataSub <- subset(data,subset=Stage==cell.subset) }
      if(cell.assign=="cluster") { dataSub <- subset(data,subset=Cluster==cell.subset) }
      if(cell.assign=="adhoc") { dataSub <- subset(data,subset=AdHoc.Assign==cell.subset) }
      if(cell.assign=="ostadhoc") { dataSub <- subset(data,subset=OstAdHoc.Assign==cell.subset) }
    }

    #de analysis
    edgr <- runEDGER(dataSub=dataSub,
                     cell.assign=cell.assign,
                     cell.subset=cell.subset,
                     genes.no.mito.ribo=genes.no.mito.ribo,
                     filter.arg=filter.arg,
                     min.count=min.count)

    edgr.total <- rbind(edgr.total,edgr)

    print(dim(edgr.total))

    saveRDS(edgr.total, file=paste0("./../data/de-data/DEedgeR.",name,".",assign,".rds"))

  }

}


