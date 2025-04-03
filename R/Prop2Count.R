
#' The Prop2Count Class
#'
#' @slot aggregated_counts genes x replicates matrix of integer counts
#' @slot proportions genes x replicates matrix of the proportion of cells with > 0 counts
#' @slot nCells vector of the number of cells per replicate
#' @slot prop2count genes x replicates matrix of transformed proportions as doubles, or integers (if round == T)
#' @slot replicate.data replicates x features data.frame with replicate level descriptions (original identifiers, total counts, etc...)
#' @slot meta.data replicates x features data.frame with replicate level descriptions (original identifiers, total counts, etc...)
#'
#' @exportClass Prop2Count
#' @importFrom methods setClass new
#' @importFrom methods is
#' @importClassesFrom Matrix dgCMatrix
#'
Prop2Count <- methods::setClass(
  Class = "Prop2Count",
  slots = c(
    aggregated_counts = "matrix",
    proportions = "matrix",
    nCells = "table",
    prop2count = "matrix",
    replicate.data = "data.frame",
    meta.data = "list"
  )
)

#' Show method for Prop2Count objects
#'
#' @param object A Prop2Count object
#'
#' @importFrom methods show
#' @exportMethod show


setMethod(
  f = "show",
  signature = "Prop2Count",
  definition = function(object) {
    cat("Prop2Count object with data for", nrow(object@proportions),
        "genes across", ncol(object@proportions),
        "replicates\n\n")
    cat("total genes: ", nrow(object@proportions), "\n")
    cat("total replicates: ", ncol(object@proportions), "\n")
    invisible(x = NULL)
  }
)

#' Create a Prop2Count object
#'
#' @param object A Seurat or Single Cell Experiment object
#' @param idents The colname(s) in the metadata to use as identifiers for aggregation & calculating proportions
#' @param sep A character to separate the ident names, if there's more than one metadata column being used
#' @param min_cells A number, the minimum number of barcodes in a unique replicate (see idents) that will be aggregated
#' @param round A boolean value, indicating whether or not to round the transformed proportions output down to an integer
#'
#' @return A Prop2Count object with transformed counts
#' @export
create_Prop2Count <- function(object,
                              idents,
                              sep = '__',
                              min_cells = 10,
                              round = FALSE) {
  if(!is(object, 'Seurat') & !is(object, 'SingleCellExperiment')) {
    stop("Input object must be a Seurat or Single Cell Experiment object")
  }

  if(is(object, 'Seurat')){

    data_object <- get_replicate_counts_proportions_Seurat(
      object = object,
      idents = idents,
      min_cells = min_cells,
      sep = sep
    )
  }

  if(is(object, 'SingleCellExperiment')){

    data_object <- get_replicate_counts_proportions_SCE(
      object = object,
      idents = idents,
      min_cells = min_cells,
      sep = sep
    )
  }


  prop2count_ <- prop2count(
    proportions_matrix = data_object$proportions,
    nCell_vector = data_object$nCells,
    round = round
  )

  result <- methods::new(
    Class = "Prop2Count",
    aggregated_counts = data_object$counts,
    proportions = data_object$proportions,
    nCells = data_object$nCells,
    prop2count = prop2count_,
    replicate.data = data_object$replicate.data,
    meta.data = data_object$meta.data
  )

  return(result)
}






