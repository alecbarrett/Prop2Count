
#' Converts proportions to a counts-like space
#'
#' @param proportions_matrix a genes x replicates matrix of proportions from 0,1.
#' @param nCell_vector a named vector with integers showing the number of cells in each replicate.
#' @param round a boolean value indicating whether or not to round the output down to an integer value
#' @export
#' @return a matrix of values from 0,infinity to be used as counts for downstream applications
prop2count <- function(proportions_matrix, nCell_vector, round = FALSE){

  if(ncol(proportions_matrix) != length(nCell_vector)){
    stop('The number of entries in the nCell vector, and the number of samples in the proportions matrix, do not match up')
  }
  if(sum(proportions_matrix > 1) > 0 | sum(proportions_matrix < 0) > 0){
    stop('proportions input needs to be in the range of 0 and 1')
  }
  if(sum(proportions_matrix==1) > 0 ){
    proportions_matrix <- prevent_infinite(proportions_matrix, nCell_vector)
  }

  prop2counts_mat <- (proportions_matrix/(1-proportions_matrix)) %*% diag(nCell_vector)

  colnames(prop2counts_mat) <- colnames(proportions_matrix)

  if(round){
    prop2counts_mat <- round(prop2counts_mat)
  }

  return(prop2counts_mat)

}

