

#' Collect aggregated data from Seurat Object, an internal function in Prop2Count
#'
#' @param object A Seurat object.
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


get_replicate_counts_proportions_Seurat <- function(object,
                                                    idents,
                                                    min_cells = 10,
                                                    sep = '__'){


  if(!is(object, 'Seurat')){
    stop('must input Seurat object')
  }


  major_version <- strsplit(as.character(object@version), '\\.')[[1]][1]


  if(major_version == "4"){
    stop('Support for Seurat V4 objects will resume shortly')
  }

  if(major_version == "5"){
    ret <- get_replicate_counts_proportions_Seurat_V5(
      object,
      idents,
      min_cells = 10,
      sep = '__'
    )

  }

  }


#' Collect aggregated data from Seurat Object V4, an internal function in Prop2Count
#'
#' @param object A Seurat object.
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

get_replicate_counts_proportions_Seurat_V5 <- function(object,
                                                    idents,
                                                    min_cells = 10,
                                                    sep = '__'){

  metadata <- object@meta.data
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
    Matrix::rowSums(object@assays$RNA@counts[,ident == x])
  })


  proportions <- sapply(replicates, function(x){
    Matrix::rowSums(object@assays$RNA@counts[,ident == x] > 0)/nCells[x]
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
                                      original_name = deparse(substitute(Seurat_object))                   )
  )

  return(ret_list)
}




#' Collect aggregated data from SingleCellExperiment object, an internal function in Prop2Count
#'
#' @param object A SingleCellExperiment object.
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
#' @importFrom SummarizedExperiment colData

get_replicate_counts_proportions_SCE <- function(object,
                                                       idents,
                                                       min_cells = 10,
                                                       sep = '__'){

  metadata <- colData(object) |> as.data.frame()
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
    Matrix::rowSums(counts(object)[,ident == x])
  })


  proportions <- sapply(replicates, function(x){
    Matrix::rowSums(counts(object)[,ident == x] > 0)/nCells[x]
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
                                      original_name = deparse(substitute(object))

                   )
  )

  return(ret_list)
}


