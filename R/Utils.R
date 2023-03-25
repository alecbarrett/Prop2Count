

#' Checks for instances where proportions == 1, and moves the value down to prevent infinite counts
#'
#' @param proportions_matrix a genes x replicates matrix of proportions from 0,1.
#' @param nCell_vector a named vector with integers showing the number of cells in each replicate.
#'
#' @return A matrix of proportions with all proportions less than 1, modified proportions will be returned as nCells-1/nCells
#'
#' @export
#'
prevent_infinite <- function(proportions_matrix, nCell_vector){

  if(0 %in% nCell_vector){
    stop('Cannot use an input where the number of cells in a replicate is 0')
  }
  if(NA %in% nCell_vector){
    stop('NAs found in nCell vector')
  }

  new_prop <- sapply(colnames(proportions_matrix), function(x){

    prop <- proportions_matrix[,x, drop = F]

    nCell <- nCell_vector[x]
    nCell_minus_1 <- nCell -1

    max_prop <- nCell_minus_1/nCell

    prop[prop==1] <- max_prop

    return(prop)

  })

  rownames(new_prop) <- rownames(proportions_matrix)

  return(new_prop)

}

#' Collect aggregated data from Seurat Object, an internal function in Prop2Count
#'
#' @param Seurat_object A Seurat object.
#' @param idents The colname(s) in the metadata to use as identifiers for aggregation & calculating proportions, what I call replicates here
#' @param min_cells A number, the minimum number of barcodes in a unique replicate (see idents) that will be aggregated, fewer than this number will be excluded from the data
#' @param sep A character to separate the ident names, if there's more than one metadata column being used. Pick a pretty unique thing that's easy to separate on later
#'
#' @return A list containing the following items:
#' \itemize{
#'    \item `counts`: the total counts per replicate in each gene
#'    \item `proportions`: the proportion of cells in the replicate with >0 counts in each gene
#'    \item `nCells`: the total cells in the replicate
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#'    L4_CeNGEN_data <- Data(L4_CeNGEN_data)
#'    example_list <- get_replicate_counts_proportions_Seurat(L4_CeNGEN_data, c('Cell.type', 'Experiment'))
#' }


get_replicate_counts_proportions_Seurat <- function(Seurat_object, idents, min_cells = 10, sep = '__'){
  if(length(idents > 1)){
    ident <- paste(idents, sep = sep)
    Seurat_object@meta.data <- tidyr::unite(Seurat_object@meta.data, ident, idents, sep = sep)
  }else{ident <- idents}

  nCells <- Seurat_object$ident |> table()
  nCells <- nCells[nCells >=10]
  replicates <- names(nCells) |> sort()
  nCells <- nCells[replicates]


  counts <- sapply(replicates, function(x){
    Matrix::rowSums(Seurat_object@assays$RNA@counts[,Seurat_object$ident == x])
  })


  proportions <- sapply(replicates, function(x){
    Matrix::rowSums(Seurat_object@assays$RNA@counts[,Seurat_object$ident == x] > 0)/nCells[x]
  })

  return(list('counts' = counts, 'proportions' = proportions, 'nCells' = nCells))

}



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




#' plot the relationship between proportions and log-counts in the original dataset
#'
#' @param Prop2Count_object a genes x replicates matrix of proportions from 0,1.
#' @param replicate Either 'highest' or 'random' which will select a single replicate to plot. 'highest' returns the replicate with the most cells, and 'random' returns a random replicate
#' @param model a character string indicating which transformation model to use. Current options are 'logit' and 'poisson'
#' @export
#' @return a matrix of values from 0,infinity to be used as counts for downstream applications
Plot_prop2count <- function(Prop2Count_object, replicate = 'highest', model = 'logit'){

  counts <- Prop2Count_object@aggregated_counts
  proportions <- Prop2Count_object@proportions
  nCells <- Prop2Count_object@nCells


  if(replicate == 'highest'){
    replicate = nCells[nCells == nCells |> max()] |> names()
  }
  if(replicate == 'random'){
    replicate = sample(names(nCells), size = 1)
  }


  line_prop <- seq(0.001,.999, 0.001)
  if(model == 'logit'){
    line_model <- log(line_prop/(1-line_prop))
  }
  if(model == 'poisson'){
    line_model <-log(log(1/(1-line_prop)))
  }

  nCell_replicate <- nCells[replicate]
  counts_replicate <- counts[,replicate]
  proportions_replicate <- proportions[, replicate]

  counts_replicate <- counts_replicate[proportions_replicate > 0 & proportions_replicate < 1]
  proportions_replicate <- proportions_replicate[proportions_replicate > 0 & proportions_replicate < 1]

  plot(proportions_replicate, log(counts_replicate), main = paste(replicate, model, 'model'),
       xlab = 'proportion of non-zero cells', ylab = 'log(counts)')
  lines( x = line_prop,
         y = (line_model) + (log(nCell_replicate)) , col = 'red', lwd = 2)


}

