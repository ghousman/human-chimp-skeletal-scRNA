###compute correlation matrix between sample and cell lines in cellAtlas
###object: Seurat object
###cell.line: cell Atlas expression matrix
###CorMethod: correlation method to use
ComputeCorMat=function(object,cell.line,CorMethod = "pearson"){

	data=GetAssayData(object)
    data_norm <- apply(data, 2, function(x) x/sum(x))
    data_libsize <- apply(data_norm, 2, function(x) x * 1e+06)

    Data.use <- data_libsize

    gene.id <- rownames(Data.use)
    gene.ref <- rownames(cell.line)
    common <- intersect(gene.id, gene.ref)
    logxx <- apply(Data.use[common, ], 2, function(x) {
        log(x + 0.1)
    })
    selected.cell.line <- apply(cell.line[common, ], 2, function(x) {
        x - mean(x)
    })
    n <- ncol(logxx)
    m <- ncol(cell.line)
    cor.mat <- matrix(0, nrow = n, ncol = m)
    for (j in 1:m) {
        for (i in 1:n) {
            cor.mat[i, j] <- cor(logxx[, i], selected.cell.line[, 
                j], method = CorMethod)
        }
    }
    rownames(cor.mat) <- colnames(Data.use)
    colnames(cor.mat) <- colnames(cell.line)
	return(cor.mat)
}

###Assign clustered cell atlas as cell type to all cells in the sample
###object: seurat object
###cor.mat: correlation matrix
###K: number of clusters generated from cell atlas
###topCorrelation: number of highest correlated cell lines in cell atlas to be assigned, set to 1 to extract only the top correlated cell line
AssignClusteredLabel=function (object,cor.mat,K=10,topCorrelation=5,	dist.method="euclidean",hclust.method="ward.D"){
	
	dist_mat <- dist(t(cor.mat), method = dist.method)
	hclust_avg_cellline <- hclust(dist_mat, method = hclust.method)
	#plot(hclust_avg_cellline)

	cut_cellline <- cutree(hclust_avg_cellline, k = K)
	cellline_abstract=sapply(strsplit(as.character(colnames(cor.mat)), "\\_"), "[[", 2)
	clusteredCellline=rbind(cut_cellline,cellline_abstract)
	nametable=rep('n',K)
	for(i in 1:K){
		indexy=as.character(i)
		nametable[i]=paste0(names( sort(table(clusteredCellline[2,clusteredCellline[1,]==indexy]),decreasing=TRUE)),collapse='_' )
	}
	assignm=apply(cor.mat,1,function(x){
		indexy=names(which.max(table(clusteredCellline[1,names(sort(x,decreasing=T)[1:topCorrelation])])))
		y=paste0(names( sort(table(clusteredCellline[2,clusteredCellline[1,]==indexy]),decreasing=TRUE)),collapse='_' )
		return(y)
	})
	object$clusteredCellAtlas=assignm
	return(object)
}

#load cell atlas data matrix
load("cell_atlas_ref_panel")

#Compute the correlation matrix use a seurat object and cell line matrix as input
cor.mat=ComputeCorMat(object,cell.line)

#Add a metadata column "clusteredCellAtlas" to the object
object=AssignClusteredLabel (object,cor.mat)
