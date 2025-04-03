

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
  if(is.null(names(nCell_vector))){
    stop('nCell_vector must be named')
  }
  if(!all(names(nCell_vector) %in% colnames(proportions_matrix))){
    stop('names in the ncell vector and proportions matrix columns must match')
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





#' plot the relationship between proportions and log-counts in the original dataset
#'
#' @param Prop2Count_object a genes x replicates matrix of proportions from 0,1.
#' @param replicate Either 'highest' or 'random' which will select a single replicate to plot. 'highest' returns the replicate with the most cells, and 'random' returns a random replicate
#' @param model a character string indicating which transformation model to use. Current options are 'logit' and 'poisson'
#'
#' @importFrom graphics plot lines
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

