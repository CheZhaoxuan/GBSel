#' @title construct a kinship matrix by VanRaden method
#'
#' @description
#' \code{VanRaden.K}use VanRaden method to construct the kinship matrix
#'
#' @param geno a genotype matrix of type numeric
#' @param maf a maf is allele frequency of second allele vector
#'
#' @return a kinship matrix
#' @export
#'
#' @examples
#' data(example_dataset)
#' X <- genotype
#' G <- VanRaden.K(geno = X,maf=NULL)
#'


VanRaden.K <- function(geno,maf=NULL){
	#--------------------------------------------------------------------------------------------------------
	# Object: To calculate the VanRaden1 matrix
	# Method source:
	# VanRaden PM. Efficient methods to compute genomic predictions, J Dairy Sci 2008;91:4414-4423
	#
	# Input:
	# geno: genotype in numeric format, pure 0, 1, 2 matrix; geno is n individual rows by m snps columns
	# maf: a maf is allele frequency of second allele vector
	#
	# Output:
	# VanRaden1 matrix
	#--------------------------------------------------------------------------------------------------------
    if (!is.matrix(geno)) stop("geno should be a matrix")
    if (!is.numeric(geno)) stop ("geno contains non numeric values")
    if (any(is.na(geno))) warning("geno contains NA values")
    if (is.null(maf)) maf<-colMeans(geno,na.rm=T)/2
    names(maf)<-colnames(geno)
    if (ncol(geno)!=length(maf)) stop("length of maf vector should equal number of columns in geno")
    if (sum(names(maf)!=colnames(geno))!=0) stop ("SNP names on vector maf and geno don't match")
    if (is.na(sum(maf))) stop ("allelic maf vector can't have NA values") 
	
	#Calculate Z matrix
	M <- geno
	P <- maf
	Z <- sweep(M,2,2*P)
	Z <- as.matrix(Z)

	#G=tcrossprod((geno), (geno))
	G <- tcrossprod(Z,Z)

	#Adjust
	Adjust <- function(x){x*(1-x)}
	Adj <- 2*sum(Adjust(P))
	VanRaden1 <- G/Adj

	return(VanRaden1)}
