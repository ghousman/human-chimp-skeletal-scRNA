#Read in command line arguments as list of character vectors
args=(commandArgs(TRUE))

#Check if arguments are passed and cycle through to evaluate each element
if(length(args)==0){
  print("No arguments supplied.")
  #supply default values
  name="noarg"
  data.source="tot"
  assign="all"
  data.type="sc"
  filter.arg=TRUE
  filter.type="pch"
  filter.param=0.2
  fml="formula1"
  lfc=FALSE
  fdr="BH"
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
print(paste0("file tag: ",name))
print(paste0("data source: ",data.source))
print(paste0("data type: ",data.type))
print(paste0("cell assignment: ",assign))
print(paste0("filter genes: ",filter.arg))
print(paste0("filter type: ",filter.type))
print(paste0("gene filter.param: ",filter.param))
print(paste0("adjuste logFC in test: ",lfc))
print(paste0("fdr method: ",fdr))

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
library(limma)
library(variancePartition)

#detectCores()
param = SnowParam(4, "SOCK", progressbar=TRUE)
register(param)

#Load data (integrate across individuals - total (15k genes))
if (data.source=="tot") {
  #data.dir <- "./../data/cellranger-data-full/data.filter.log.indv-cell.int19k.reg-tot.assign.rds"
  data.dir <- "./../data/cellranger-data-full/data.filter.log.indv-cell.int19k-tot.assign.rds"
  data <- readRDS(data.dir)
}
if (data.source=="t0") {
  #data.dir <- "./../data/cellranger-data-full/data.filter.log.indv-cell.int19k.reg-t0.assign.rds"
  data.dir <- "./../data/cellranger-data-full/data.filter.log.indv-cell.int19k-t0.assign.rds"
  data <- readRDS(data.dir)
}
if (data.source=="t1") {
  #data.dir <- "./../data/cellranger-data-full/data.filter.log.indv-cell.int19k.reg-t1.assign.rds"
  data.dir <- "./../data/cellranger-data-full/data.filter.log.indv-cell.int19k-t1.assign.rds"
  data <- readRDS(data.dir)
}
if (data.source=="t2") {
  #data.dir <- "./../data/cellranger-data-full/data.filter.log.indv-cell.int19k.reg-t2.assign.rds"
  data.dir <- "./../data/cellranger-data-full/data.filter.log.indv-cell.int19k-t2.assign.rds"
  data <- readRDS(data.dir)
}
print(data.dir)
data

#Isolate mitochondrial and ribosomal genes
genes <- rownames(data@assays$RNA@counts)
genes.mito <- c("MT-ATP6","MT-ATP8","MT-CO1","MT-CO2","MT-CO3","MT-CYB","MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5")
genes.ribo <- grep('^RP',genes,value=T)
genes.no.mito.ribo <- genes[which(!(genes %in% c(genes.mito,genes.ribo)))]
rm(genes,genes.mito,genes.ribo)

#Define Limma Voom (Dream) Function
runDREAM <- function(dataSub, data.dir, cell.assign, cell.subset, genes.no.mito.ribo, data.type, filter.arg, filter.type, filter.param, fml) {

  #make count matrix and metadata for edgeR object
  if (data.type=="pseudo") {
    counts <- c()
    metadata <- c()
    labels <- c()
    for (i in unique(dataSub@meta.data$Sample)) {
      x.lab <- i
      w <- which(dataSub@meta.data$Sample==i)
      if (length(w)>0) {
        x.spp <- dataSub@meta.data$Species[w][1]
        x.col <- dataSub@meta.data$Collection[w][1]
        x.ind <- dataSub@meta.data$Individual[w][1]
        x.rep <- dataSub@meta.data$Replicate[w][1]
        if (length(w)==1) {
          x.cnt <- dataSub@assays$RNA@counts[,w]
        } else {
          x.cnt <- Matrix::rowSums(dataSub@assays$RNA@counts[,w])
        }
        counts <- cbind(counts, x.cnt)
        metadata <- rbind(metadata, c(x.lab, x.spp, x.col, x.ind, x.rep))
        labels <- c(labels, x.lab)
      }
    }
    colnames(counts) <- labels
    rownames(metadata) <- labels
    colnames(metadata) <- c("Sample","Species","Collection","Individual","Replicate")
    metadata <- as.data.frame(metadata)
    rm(labels)
  }
  if (data.type=="sc") {
    counts <- as.matrix(GetAssayData(dataSub, assay="RNA", slot="counts"))
    metadata <- dataSub@meta.data[,c("Species","Individual","Replicate","Collection","percent.mt","Phase","Sample")]
    metadata <- metadata[colnames(counts),]
  }

  #remove mitochodrial and ribosomal genes
  counts <- counts[which(rownames(counts) %in% genes.no.mito.ribo),]

  #make edgeR object
  dge <- DGEList(counts)
  meta_dge <- dge$samples[,c("lib.size","norm.factors")]
  meta_dge <- cbind(meta_dge, metadata)
  dge$samples <- meta_dge
  rm(dataSub,counts,metadata,meta_dge)

  #filter genes
  if (filter.arg==TRUE){
    if (filter.type=="min.count"){
      keep <- filterByExpr(dge, group=dge$samples$Species, min.count=filter.param, min.total.count=15)
      table(keep)
      dge <- dge[keep, , keep.lib.sizes=FALSE]
      rm(keep)
    }
    if (filter.type=="pch"){
      dge$counts <- dge$counts[rowSums(dge$counts!=0)>=(filter.param*dim(dge$counts)[2]),]
    }
  }

  #normalize data
  dge <- calcNormFactors(dge, method="TMM")
  summary(dge$samples$norm.factors)

  #prep design matrix variables
  dge$samples$Individual <- as.factor(dge$samples$Individual)
  dge$samples$Replicate <- as.factor(as.character(dge$samples$Replicate))
  dge$samples$Species <- as.factor(dge$samples$Species)
  dge$samples$Collection <- as.factor(as.numeric(dge$samples$Collection))
  if (data.type=="sc") {
    dge$samples$Phase <- as.factor(dge$samples$Phase)
    dge$samples$Sample <- as.factor(dge$samples$Sample)
  }

  #design matrix
  if (fml=="formula1") {
    model <- "~Species+(1|Individual)+(1|Replicate)"
    formula <- ~Species+(1|Individual)+(1|Replicate)
  }
  if (fml=="formula2") {
    model <- "~Species+Collection+percent.mt+Phase+(1|Individual)"
    formula <- ~Species+Collection+percent.mt+Phase+(1|Individual)
  }
  if (fml=="formula3") {
    model <- "~Species+percent.mt+Phase+(1|Sample)"
    formula <- ~Species+percent.mt+Phase+(1|Sample)
  }
  if (fml=="formula4") {
    model <- "~Species+percent.mt+Phase+(1|Individual)+(1|Replicate)"
    formula <- ~Species+percent.mt+Phase+(1|Individual)+(1|Replicate)
  }
  if (fml=="formula5") {
    model <- "~Species+percent.mt+Phase+(1|Individual)"
    formula <- ~Species+percent.mt+Phase+(1|Individual)
  }

  #estimate weights using linear mixed model of dream (voom)
  vobjDream = voomWithDreamWeights(counts=dge, formula=formula, data=dge$samples)

  #define contrasts
  #L = getContrast(vobjDream, formula, dge$samples, coefficient="SpeciesHuman")
  #plotContrasts(L)

  #fit dream model on each gene
  fitmm = dream(vobjDream, formula, dge$samples)

  #assess differential expression output for particular contrast
  tt <- topTable(fitmm, n=Inf, coef='SpeciesHuman', adjust.method="BH", p.value=1)
  print("All DE Genes")
  print(dim(tt))
  print(table(tt$adj.P.Val<0.05))
  print(summary(decideTests(fitmm, adjust.method=fdr, p.value=0.05)))

  #add details to output
  tt$data <- rep(data.dir, dim(tt)[1])
  tt$data.type <- rep(data.type, dim(tt)[1])
  tt$cell.assign <- rep(cell.assign, dim(tt)[1])
  tt$cell.subset <- rep(cell.subset, dim(tt)[1])
  tt$gene.filter <- rep(filter.arg, dim(tt)[1])
  tt$gene.filter.args <- rep(paste0("filter.param=",filter.param), dim(tt)[1])
  tt$model <- rep(model, dim(tt)[1])
  tt$comparison <- rep(tt$comparison, dim(tt)[1])
  tt$test <- rep(tt$test, dim(tt)[1])
  tt$adjust.method <- rep(tt$adjust.method, dim(tt)[1])
  tt$gene <- rownames(tt)

  return(tt)

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

DEgene.total <- data.frame()

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
    DEgene <- runDREAM(dataSub=dataSub,
                       data.dir=data.dir,
                       cell.assign=cell.assign,
                       cell.subset=cell.subset,
                       genes.no.mito.ribo=genes.no.mito.ribo,
                       data.type=data.type,
                       filter.arg=filter.arg,
                       filter.type=filter.type,
                       filter.param=filter.param,
                       fml=fml)

    DEgene.total <- rbind(DEgene.total,DEgene)

    print(dim(DEgene.total))

    saveRDS(DEgene.total, file=paste0("./../data/de-data/DEdream.",name,".",assign,".rds"))

  }

}


