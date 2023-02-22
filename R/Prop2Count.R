#'# Class definitions

#' The Prop2Count Class
#'
#' @slot aggregated_counts genes x replicates matrix of integer counts
#' @slot proportions genes x replicates matrix of the proportion of cells with > 0 counts
#' @slot nCells vector of the number of cells per replicate
#' @slot prop2count genes x replicates matrix of transformed proportions as doubles, or integers (if round == T)
#'
#'
#' @exportClass Prop2Count
#' @importFrom methods setClass
#' @importFrom methods is
#' @importClassesFrom Matrix dgCMatrix
#'



Prop2Count = methods::setClass(
  Class = "Prop2Count",
  slots = c(
    aggregated_counts = "Matrix",
    proportions = "Matrix",
    nCells = "table",
    prop2count = "Matrix"
  )
)

setMethod(f = "show", signature = "Prop2Count", definition = function(object) {
  cat("An object of class", class(object))
  invisible(x = NULL)
})


#' Collect aggregated data from Seurat Object, an internal function in Prop2Count
#'
#' @param object A Seurat object.
#' @param idents The colname(s) in the metadata to use as identifiers for aggregation & calculating proportions, what I call replicates here
#' @param min_cells A number, the minimum number of barcodes in a unique replicate (see idents) that will be aggregated, fewer than this number will be excluded from the data
#' @param sep A character to separate the ident names, if there's more than one metadata column being used. Pick a pretty unique thing that's easy to separate on later
#' @param round A boolean value, indicating whether or not to round the transformed proportions output down to an integer. Default = FALSE
#' @param keep_counts A boolean value, indicating whether to keep the oriignal aggregated counts matrix. Default = TRUE
#'
#' @return a Prop2Count object with transformed counts

create_Prop2Count <- function(object, idents, sep = '__', min_cells = 10, round = FALSE, keep_counts = TRUE){

  if(is(tmp_sc, 'Seurat')){
    data_object <- get_replicate_counts_proportions_Seurat(Seurat_object = object,
                                                         idents = idents,
                                                         min_cells = min_cells,
                                                         sep = '__')

    prop2count_ <- prop2count(proportions_matrix = data_object$proportions,
                              nCell_vector = data_object$nCells,
                              round = round)
    if(keep_counts){
      data_object = methods::new(
        Class = "Prop2Count",
        aggregated_counts = data_object$counts,
        proportions = data_object$proportions,
        nCells = data_object$nCells,
        prop2count = prop2count_)
    }else(
      data_object = methods::new(
        Class = "Prop2Count",
        proportions = data_object$proportions,
        nCells = data_object$nCells,
        prop2count = prop2count_)
    )




  }

  return(data_object)

}







