#' @title Reduce the feature dimension of the input data
#'
#' @description
#' \code{Feature.reduce}Reduce the feature dimension of the input data
#'
#' @param X genotype in numeric format, pure 0, 1, 2 matrix; geno is n individual rows by m snps columns
#' @param Y a phenotype vector of type numeric
#' @param method the method for reducing the feature; options are GRM, PCA
#' @param dimension confirm the dimensions after the reduction
#'
#' @return X data after dimensionality reduction
#' @export
#' @importFrom stats "prcomp"
#'
#' @examples
#' data(example_dataset)
#' X <- genotype
#' X.FR <- Feature.reduce(X = X, method="GRM")
#' 


Feature.reduce <- function(X, Y=NULL, method="GRM", dimension="Ind.num"){
	#--------------------------------------------------------------------------------------------------------
	# Object: Reduce the feature dimension of the input data
	#
	# Input:
	# X: genotype in numeric format, pure 0, 1, 2 matrix; geno is n individual rows by m snps columns
	# Y: a phenotype vector of type numeric
	# method: the method for reducing the feature; options are GRM, PCA, GWAS(MLM)
	# dimension: confirm the dimensions after the reduction
	#
	# Output:
	# X data after dimensionality reduction
	#--------------------------------------------------------------------------------------------------------
	
	#Feature_index <- function(Y=Y,X=X,method="MLM",K=G,dimension=dimension){
		#Y <- cbind(c(1:length(Y)),Y)
		#if(method=="MLM"){
			#nind <- dim(X)[1]
			#nmarker <- dim(X)[2]
			#pop.map <- simer::generate.map(pop.marker = nmarker)
			#geno <- t(X)#转置，适配rMVP
			#genotype <- bigmemory::as.big.matrix(geno)#转成bigmatrix
			#Kinship <- bigmemory::as.big.matrix(K)#转成bigmatrix
			#imMVP <- rMVP::MVP(phe=Y,geno=genotype,map=pop.map,K=Kinship,nPC.MLM=NULL,priority="speed",file.output=FALSE,method=c("MLM"))
			#p_matrix <- imMVP$mlm.results[,3]
			#names(p_matrix) <- c(1:length(p_matrix))
			#if(dimension == "Ind.num"){index <- names(sort(p_matrix))[1:nind]}else{index <- names(sort(p_matrix))[1:dimension]}
		#}
		#return(as.numeric(index))
	#}

	if(!method %in% c("GRM","PCA")){stop("method shuold be one of 'GRM', 'PCA'.")}
	
	if(method == "GRM"){
		result <- VanRaden.K(geno=X)
	}
	if(method == "PCA"){
		result <- prcomp(X)$x
	}
	#if(method == "GWAS"){
		#trainnames <- names(Y)
		#if(is.null(trainnames)){stop("names(Y) is NULL!")}
		#if(is.null(rownames(X))){stop("rownames(X) is NULL!")}
		#if(!all(is.element(names(Y),rownames(X)))){stop("rownames(X) do not contain names(Y) all!")}
		#Xtr <- X[trainnames,]
		#G <- VanRaden.K(geno=Xtr)
		#index <- Feature_index(Y=Y,X=Xtr,method="MLM",K=G,dimension=dimension)
		#result <- X[,index]
	#}	
	
	return(result)
}

