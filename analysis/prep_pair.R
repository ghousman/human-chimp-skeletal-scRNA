#Load libraries
library(Seurat)
library(dplyr)
library(stringi)
library(stringr)
library(ggplot2)
library(colorspace)
library(RColorBrewer)

#Define main directory
dir <- '/project2/gilad/ghousman/skeletal-human-chimp/'

#Load batch info
batch <- read.csv(file=paste0(dir,'human-chimp-skeletal-scRNA/data/scrna-batch.csv'), header=TRUE, sep=",")

#Read in files
data <- readRDS(paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.rds"))

#Merge all cell types in each human-chimp pair
h1c1 <- merge(data[[1]], y=c(data[[2]],data[[3]]), add.cell.ids=c("H1C1.I","H1C1.M","H1C1.O"))
h1c1r2 <- merge(data[[4]], y=c(data[[5]],data[[6]]), add.cell.ids=c("H1C1r2.I","H1C1r2.M","H1C1r2.O"))
h2c2 <- merge(data[[7]], y=c(data[[8]],data[[9]]), add.cell.ids=c("H2C2.I","H2C2.M","H2C2.O"))
h3c3 <- merge(data[[10]], y=c(data[[11]],data[[12]]), add.cell.ids=c("H3C3.I","H3C3.M","H3C3.O"))
h4c4 <- merge(data[[13]], y=c(data[[14]],data[[15]]), add.cell.ids=c("H4C4.I","H4C4.M","H4C4.O"))
h5c5 <- merge(data[[16]], y=c(data[[17]],data[[18]]), add.cell.ids=c("H5C5.I","H5C5.M","H5C5.O"))
h6c6 <- merge(data[[19]], y=c(data[[20]],data[[21]]), add.cell.ids=c("H6C6.I","H6C6.M","H6C6.O"))

objects <- list(h1c1,h1c1r2,h2c2,h3c3,h4c4,h5c5,h6c6)

#Perform SCTransform normalization on each collection
for (i in 1:length(objects)) {
  print(i)
  objects[[i]] <- SCTransform(objects[[i]], conserve.memory=TRUE)
}

#Save data
saveRDS(objects, file=paste0(dir,"/human-chimp-skeletal-scRNA/data/cellranger-data-full/data.pair.sct.rds"))

