#Read in command line arguments as list of character vectors
args=(commandArgs(TRUE))

#Check if arguments are passed and cycle through to evaluate each element
if(length(args)==0){
  print("No arguments supplied.")
  #supply default values
  name="noarg"
  data.source="tot"
  data.type="pseudo"
  subset.rep="NULL"
  subset.num=FALSE
  filter.arg=TRUE
  filter.type="logcpm"
  filter.param=0
  fml=formula1
  lfc=FALSE
  fdr="BH"
  assign="all"
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
print(paste0("file tag: ",name))
print(paste0("data source: ",data.source))
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
library(BiocParallel)

#detectCores()
register(SnowParam(4, "SOCK", progressbar=TRUE))

#Load data (conservative cell filter + integrated across individuals + non-zero genes + regress out UMI/mito)
if (data.source=="tot") {
  data.dir <- "./../data/cellranger-data-full/data.filterC.log.indv-cell.intNo0.reg-tot.assign.rds"
  data <- readRDS(data.dir)
}
if (data.source=="t0") {
  data.dir <- "./../data/cellranger-data-full/data.filterC.log.indv-cell.intNo0.reg-t0.assign.rds"
  data <- readRDS(data.dir)
}
if (data.source=="t1") {
  data.dir <- "./../data/cellranger-data-full/data.filterC.log.indv-cell.intNo0.reg-t1.assign.rds"
  data <- readRDS(data.dir)
}
if (data.source=="t2") {
  data.dir <- "./../data/cellranger-data-full/data.filterC.log.indv-cell.intNo0.reg-t2.assign.rds"
  data <- readRDS(data.dir)
}
print(data.dir)
data

#Define individual-replicate sets
data@meta.data$Individual.Replicate <- paste0(data@meta.data$Individual,".",data@meta.data$Replicate)

#Isolate mitochondrial and ribosomal genes
genes <- rownames(data@assays$RNA@counts)
genes.mito <- c("MT-ATP6","MT-ATP8","MT-CO1","MT-CO2","MT-CO3","MT-CYB","MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L","MT-ND5")
genes.ribo <- grep('^RP',genes,value=T)
genes.no.mito.ribo <- genes[which(!(genes %in% c(genes.mito,genes.ribo)))]
rm(genes,genes.mito,genes.ribo)

#Define Limma Voom (Dream) Function
runDREAM <- function(dataSub, data.dir, subset.rep, cell.assign, cell.subset, genes.no.mito.ribo, filter.arg, filter.type, filter.param, fml) {

  #make count matrix and metadata for edgeR object
  if (data.type=="pseudo") {
    counts <- c()
    metadata <- c()
    labels <- c()
    for (i in unique(dataSub@meta.data$Individual.Replicate)) {
      x.lab <- i
      w <- which(dataSub@meta.data$Individual.Replicate==i)
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
    colnames(metadata) <- c("Individual.Replicate","Species","Collection","Individual","Replicate")
    metadata <- as.data.frame(metadata)
    rm(labels)
  }
  if (data.type=="sc") {
    counts <- as.matrix(GetAssayData(dataSub, assay="RNA", slot="counts"))
    metadata <- dataSub@meta.data[,c("Species","Individual.Replicate","Individual","Replicate","Collection","percent.mt","Phase","Sample")]
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
    if (filter.type=="logcpm"){
      keep <- rowMeans(edgeR::cpm(dge,log=TRUE,prior.count=0.25))>filter.param
      dge$counts <- dge$counts[keep,]
    }
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
  if (fml=="formula1" & subset.rep %in% c("keepRep1","keepRep2")) {
    model <- "~Species+(1|Individual)"
    formula <- ~Species+(1|Individual)
  }
  if (fml=="formula2") {
    model <- "~Species+Collection+percent.mt+Phase+(1|Individual)"
    formula <- ~Species+Collection+percent.mt+Phase+(1|Individual)
  }
  if (fml=="formula3") {
    model <- "~Species+percent.mt+Phase+(1|Sample)"
    formula <- ~Species+percent.mt+Phase+(1|Sample)
  }
  if (fml=="formula4" & data.type=="pseudo") {
    model <- "~Species+(1|Individual)+(1|Replicate)"
    formula <- ~Species+(1|Individual)+(1|Replicate)
  }
  if (fml=="formula4" & data.type=="pseudo" & (subset.rep=="keepRep1" | subset.rep=="keepRep2")) {
    model <- "~Species+(1|Individual)"
    formula <- ~Species+(1|Individual)
  }
  if (fml=="formula4" & data.type=="sc") {
    model <- "~Species+percent.mt+Phase+(1|Individual)+(1|Replicate)"
    formula <- ~Species+percent.mt+Phase+(1|Individual)+(1|Replicate)
  }
  if (fml=="formula4" & data.type=="sc" & (subset.rep=="keepRep1" | subset.rep=="keepRep2")) {
    model <- "~Species+percent.mt+Phase+(1|Individual)"
    formula <- ~Species+percent.mt+Phase+(1|Individual)
  }
  if (fml=="formula5" & data.type=="pseudo") {
    model <- "~Species+(1|Individual)"
    formula <- ~Species+(1|Individual)
  }
  if (fml=="formula5" & data.type=="sc") {
    model <- "~Species+percent.mt+Phase+(1|Individual)"
    formula <- ~Species+percent.mt+Phase+(1|Individual)
  }
  if (fml=="formula6" & data.type=="pseudo") {
    model <- "~Species+(1|Individual)"
    formula <- ~Species+(1|Individual)
  }
  if (fml=="formula6" & data.type=="pseudo" & (subset.rep=="keepRep1" | subset.rep=="keepRep2")) {
    model <- "~Species"
    formula <- ~Species
  }
  if (fml=="formula6" & data.type=="pseudo" & (cell.assign=="ostcluster0.50" | cell.assign=="ostcluster0.50.T2") & cell.subset=="Osteogenic.c2") {
    model <- "~Species"
    formula <- ~Species
  }
  if (fml=="formula6" & data.type=="pseudo" & (cell.assign=="ostcluster0.75" | cell.assign=="ostcluster0.75.T2") & cell.subset=="Osteogenic.c6") {
    model <- "~Species"
    formula <- ~Species
  }
  if (fml=="formula6" & data.type=="pseudo" & (cell.assign=="ostadhocX" | cell.assign=="ostadhocX.T2") & cell.subset=="mineralizing osteoblast") {
    model <- "~Species"
    formula <- ~Species
  }
  if (fml=="formula6" & data.type=="sc") {
    model <- "~Species+percent.mt+Phase+(1|Individual.Replicate)"
    formula <- ~Species+percent.mt+Phase+(1|Individual.Replicate)
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
  print(table(tt$adj.P.Val<0.01))
  print(summary(decideTests(fitmm, adjust.method=fdr, p.value=0.01)))

  #rename final column of toptable to enable rbind()
  #column name is "B" for mineralizing osteoblasts in ostadhocX and ostadhocX.T2
  colnames(tt)[6] <- "z.std"

  #add details to output
  tt$data <- rep(data.dir, dim(tt)[1])
  tt$data.type <- rep(data.type, dim(tt)[1])
  tt$subset.rep <- rep(subset.rep, dim(tt)[1])
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
if (assign=="stage") {
  subset.list <- list(stage=list("stage-Time 0",
                                 "stage-Time 1",
                                 "stage-Time 2"))
}
if (assign=="cluster"  & data.source=="tot") {
  subset.list <- list(clust=list("cluster-iPSC.c1",
                                 "cluster-iPSC.c2",
                                 "cluster-iPSC.c3",
                                 "cluster-MSC.c1",
                                 "cluster-Osteogenic.c1",
                                 "cluster-Osteogenic.c2"),
                      combI=c("cluster-iPSC.c1_iPSC.c2_iPSC.c3"),
                      combO=c("cluster-Osteogenic.c1_Osteogenic.c2"))
}
if (assign=="adhoc") {
  subset.list <- list(adhoc=list("adhoc-iPSC",
                                 "adhoc-MSC",
                                 "adhoc-Osteogenic"))
}
if (assign=="ostcluster0.50"  & data.source=="tot") {
  subset.list <- list(clust=list("ostcluster0.50-Osteogenic.c1",
                                 "ostcluster0.50-Osteogenic.c2",
                                 "ostcluster0.50-Osteogenic.c3",
                                 "ostcluster0.50-Osteogenic.c4"))
}
if (assign=="ostcluster0.75"  & data.source=="tot") {
  subset.list <- list(clust=list("ostcluster0.75-Osteogenic.c1",
                                 "ostcluster0.75-Osteogenic.c2",
                                 "ostcluster0.75-Osteogenic.c3",
                                 "ostcluster0.75-Osteogenic.c4",
                                 "ostcluster0.75-Osteogenic.c5",
                                 "ostcluster0.75-Osteogenic.c6",
                                 "ostcluster0.75-Osteogenic.c7",
                                 "ostcluster0.75-Osteogenic.c8"))
}
if (assign=="ostcluster0.50.T2"  & data.source=="tot") {
  subset.list <- list(clust=list("ostcluster0.50.T2-Osteogenic.c1",
                                 "ostcluster0.50.T2-Osteogenic.c2",
                                 "ostcluster0.50.T2-Osteogenic.c3",
                                 "ostcluster0.50.T2-Osteogenic.c4"))
}
if (assign=="ostcluster0.75.T2"  & data.source=="tot") {
  subset.list <- list(clust=list("ostcluster0.75.T2-Osteogenic.c1",
                                 "ostcluster0.75.T2-Osteogenic.c2",
                                 "ostcluster0.75.T2-Osteogenic.c3",
                                 "ostcluster0.75.T2-Osteogenic.c4",
                                 "ostcluster0.75.T2-Osteogenic.c5",
                                 "ostcluster0.75.T2-Osteogenic.c6",
                                 "ostcluster0.75.T2-Osteogenic.c7",
                                 "ostcluster0.75.T2-Osteogenic.c8"))
}
if (assign=="ostadhoc") {
  subset.list <- list(ostadhoc=list("ostadhoc-preosteoblast",
                                    "ostadhoc-osteoblast",
                                    "ostadhoc-embedding osteoblast",
                                    "ostadhoc-mineralizing osteoblast",
                                    "ostadhoc-maturing osteocyte"))
}
if (assign=="ostadhocX") {
  subset.list <- list(ostadhoc=list("ostadhocX-preosteoblast",
                                    "ostadhocX-osteoblast",
                                    "ostadhocX-embedding osteoblast",
                                    "ostadhocX-mineralizing osteoblast",
                                    "ostadhocX-maturing osteocyte"))
}
if (assign=="ostadhoc.T2") {
  subset.list <- list(ostadhoc=list("ostadhoc.T2-preosteoblast",
                                    "ostadhoc.T2-osteoblast",
                                    "ostadhoc.T2-embedding osteoblast",
                                    "ostadhoc.T2-mineralizing osteoblast",
                                    "ostadhoc.T2-maturing osteocyte"))
}
if (assign=="ostadhocX.T2") {
  subset.list <- list(ostadhoc=list("ostadhocX.T2-preosteoblast",
                                    "ostadhocX.T2-osteoblast",
                                    "ostadhocX.T2-embedding osteoblast",
                                    "ostadhocX.T2-mineralizing osteoblast",
                                    "ostadhocX.T2-maturing osteocyte"))
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
                          subset=Cluster0.05==str_split(cell.subset,"_")[[1]][1]|
                            Cluster0.05==str_split(cell.subset,"_")[[1]][2])
      }
      else {
        dataSub <- subset(data,
                          subset=Cluster0.05==str_split(cell.subset,"_")[[1]][1]|
                            Cluster0.05==str_split(cell.subset,"_")[[1]][2]|
                            Cluster0.05==str_split(cell.subset,"_")[[1]][3])
      }
    }

    else {
      if(cell.assign=="stage")             { dataSub <- subset(data,subset=Stage==cell.subset) }
      if(cell.assign=="cluster")           { dataSub <- subset(data,subset=Cluster0.05==cell.subset) }
      if(cell.assign=="adhoc")             { dataSub <- subset(data,subset=AdHoc.Assign.thresh0==cell.subset) }
      if(cell.assign=="ostcluster0.50")    { dataSub <- subset(data,subset=Cluster0.50==cell.subset) }
      if(cell.assign=="ostcluster0.75")    { dataSub <- subset(data,subset=Cluster0.75==cell.subset) }
      if(cell.assign=="ostcluster0.50.T2") { dataSub <- subset(data,subset=Cluster0.50==cell.subset & Stage=="Time 2") }
      if(cell.assign=="ostcluster0.75.T2") { dataSub <- subset(data,subset=Cluster0.75==cell.subset & Stage=="Time 2") }
      if(cell.assign=="ostadhoc")          { dataSub <- subset(data,subset=OstAdHoc.Assign.4gene.thresh0==cell.subset) }
      if(cell.assign=="ostadhocX")         { dataSub <- subset(data,subset=OstAdHoc.Assign.x.4gene.thresh0==cell.subset) }
      if(cell.assign=="ostadhoc.T2")       { dataSub <- subset(data,subset=OstAdHoc.Assign.4gene.thresh0==cell.subset & Stage=="Time 2") }
      if(cell.assign=="ostadhocX.T2")      { dataSub <- subset(data,subset=OstAdHoc.Assign.x.4gene.thresh0==cell.subset & Stage=="Time 2") }
    }

    #further define data subsampling based on sample replicates
    if(subset.rep=="keepRep1") {
      print(subset.rep)
      dataSub <- subset(dataSub,subset=Individual.Replicate!="H1.1")
      dataSub <- subset(dataSub,subset=Individual.Replicate!="C1.1")
    }
    if(subset.rep=="keepRep2") {
      print(subset.rep)
      dataSub <- subset(dataSub,subset=Individual.Replicate!="H1.2")
      dataSub <- subset(dataSub,subset=Individual.Replicate!="C1.2")
    }

    #further define data subsampling based on number of cells to randomly sample for DE analysis
    if(subset.num==TRUE) {
      human.cells <- rownames(dataSub@meta.data)[which(dataSub@meta.data$Species=="Human")]
      chimp.cells <- rownames(dataSub@meta.data)[which(dataSub@meta.data$Species=="Chimp")]
      min.cell.count <- min(c(length(human.cells),length(chimp.cells)))
      print(paste0("Minimum number of cells in comparative groups: ",min.cell.count))
      h.cell.sub <- human.cells[sample(1:length(human.cells), min.cell.count, replace=FALSE)]
      c.cell.sub <- chimp.cells[sample(1:length(chimp.cells), min.cell.count, replace=FALSE)]
      dataSub <- subset(dataSub,cells=c(h.cell.sub,c.cell.sub))
    }

    #de analysis
    DEgene <- runDREAM(dataSub=dataSub,
                       data.dir=data.dir,
                       subset.rep=subset.rep,
                       cell.assign=cell.assign,
                       cell.subset=cell.subset,
                       genes.no.mito.ribo=genes.no.mito.ribo,
                       filter.arg=filter.arg,
                       filter.type=filter.type,
                       filter.param=filter.param,
                       fml=fml)

    DEgene.total <- rbind(DEgene.total,DEgene)

    print(dim(DEgene.total))

    saveRDS(DEgene.total, file=paste0("./../data/de-data/DEdream.",name,".",assign,".rds"))

  }

}


