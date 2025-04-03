

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
#' @importFrom stats median
#' @export
#'
#' @examples
#' \dontrun{
#'    L4_CeNGEN_data <- Data(L4_CeNGEN_data)
#'    example_list <- get_replicate_counts_proportions_Seurat(L4_CeNGEN_data, c('Cell.type',
#'                                                            'Experiment'))
#' }


get_replicate_counts_proportions_Seurat <- function(Seurat_object,
                                                    idents,
                                                    min_cells = 10,
                                                    sep = '__'){

  metadata <- Seurat_object@meta.data
  if(length(idents) > 1){
    ident <- apply(metadata[,idents], 1, paste, collapse = sep)
    metadata_unique <- metadata[,idents] |> unique()
    rownames(metadata_unique) <- ident |> unique()
  }else{
    ident <- metadata[[ident]]
    metadata_unique <- metadata[, ident, drop = F]|> unique()
    rownames(metadata_unique) <- ident |> unique()
  }

  nCells <- ident |> table()
  nCells <- nCells[nCells >= min_cells]
  replicates <- names(nCells) |> sort()
  nCells <- nCells[replicates]
  metadata_unique <- metadata_unique[replicates,]

  counts <- sapply(replicates, function(x){
    Matrix::rowSums(Seurat_object@assays$RNA@counts[,ident == x])
  })


  proportions <- sapply(replicates, function(x){
    Matrix::rowSums(Seurat_object@assays$RNA@counts[,ident == x] > 0)/nCells[x]
  })

  metadata_unique$nCells <- nCells

  ret_list <- list('counts' = counts,
                   'proportions' = proportions,
                   'nCells' = nCells,
                   'replicate.data' =  metadata_unique,
                   'meta.data' = list(total_cells = sum(nCells),
                                      median_cells_per_replicate = median(nCells),
                                      total_replicates = ncol(proportions),
                                      min_cells_threshold = min_cells,
                                      original_name = deparse(substitute(Seurat_object)),
                                      orig_object_type = 'Seurat',
                                      object_package_version = Seurat_object@version
                   )
  )

  return(ret_list)
}

