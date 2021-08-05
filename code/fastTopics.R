#Read in command line arguments as list of character vectors
args=(commandArgs(TRUE))

#Check if arguments are passed and cycle through to evaluate each element
if(length(args)==0){
  print("No arguments supplied.")
  #supply default values
  clustnum=3
  data.source="tot"
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
print(paste0("k: ",clustnum))
print(paste0("data source: ",data.source))

#Load libraries
library(Seurat)
library(CountClust)
library(classtpx)
library(fastTopics)

#Load data
scrna <- paste0("../data/cellranger-data-full/data.filterC.log.indv-cell.intNo0.reg-",data.source,".rds")
data <- readRDS(scrna)

#Prepare data (unintegrated count data)
counts <- as.matrix(data@assays$RNA@counts) #do not use integrated data
counts <- counts[rowSums(counts>0)>0,]      #remove genes with 0 counts across cells
counts <- t(as.matrix(counts))              #fix matrix orientation (barcodes x features)
rm(data)

#Fit fastTopics model and save data
fit <- fit_poisson_nmf(counts, k=clustnum, numiter=100)
saveRDS(fit, file=paste0("../data/gom-data/topics",clustnum,".data.filterC.log.indv-cell.intNo0.reg-",data.source,".rds"))
