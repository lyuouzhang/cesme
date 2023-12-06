#' @name buettner
#' @aliases buettner
#'
#' @title Real data example for CESME
#' @description The data contains gene expression data of 182 objects on 8989 genes.
#' @usage data(buettner)
#' @docType data
#' @format
#' \describe{
#'  \item{Rows}{Objects}
#'  \item{Columns}{Genes}
#' }
#'
#'
#' @details
#' This is the gene expression data from Buettner
#' @importFrom stats cov kmeans model.matrix qnorm sd
#' @import Matrix
#' @useDynLib cesme, .registration=TRUE
#' @export
#' 
#' @keywords
#' dataset
#' @examples
#' data(buettner)
#' require(kernlab)
#' require(bigmemory)
#' buettner_label = data_buettner[,1]
#' buettner_covariate = data_buettner[,-1]
#' n = dim(buettner_covariate)[1]
#' p = dim(buettner_covariate)[2]
#' K = max(buettner_label)
#' z.initial = as.numeric(specc(buettner_covariate, centers = K))
#' npn.def = npn.clust.bic.large(buettner_covariate, K=K)
#' npn.clust = npn.def$z[npn.def$iter+1,,1]



data_example = function(){
  data(buettner)
  require(kernlab)
  require(bigmemory)
  buettner_label = data_buettner[,1]
  buettner_covariate = data_buettner[,-1]
  n = dim(buettner_covariate)[1]
  p = dim(buettner_covariate)[2]
  K = max(buettner_label)
  z.initial = as.numeric(specc(buettner_covariate, centers = K))
  npn.def = npn.clust.bic.large(buettner_covariate, K=K)
  npn.clust = npn.def$z[npn.def$iter+1,,1]
}
