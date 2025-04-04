


## load Seurat example

example_seurat <- readRDS(system.file("extdata", "Vignette_Seurat_object_V5.rds", package = "Prop2CountR"))

example_sce <- Seurat::as.SingleCellExperiment(example_seurat)

test_that('seurat object processed', {

  expect_no_error(get_replicate_counts_proportions_Seurat(example_seurat, c('Cell.type', 'Experiment')))

})


test_that('Seurat idents that do not exist throw an error', {
  expect_error(get_replicate_counts_proportions_Seurat(example_seurat, 'NonExistentColumn'))
})

metadata <- example_seurat@meta.data
cell_counts <- table(metadata$Cell.type)
max_cells <- max(cell_counts) + 10

test_that('too high of cell minimum throws an error', {

  metadata <- example_seurat@meta.data
  cell_counts <- table(metadata$Cell.type)
  max_cells <- max(cell_counts) + 100

  expect_error(get_replicate_counts_proportions_Seurat(example_seurat, 'Cell.type', max_cells))


})

test_that('get_replicate_counts_proportions_SCE processes data correctly', {
  result_sce <- get_replicate_counts_proportions_SCE(example_sce, c('Cell.type', 'Experiment'))

  expect_type(result_sce, "list")
  #expect_true(all(c('counts', 'proportions', 'nCells', 'replicate.data', 'meta.data') %in% names(result_sce)))

  result_single <- get_replicate_counts_proportions_SCE(example_sce, c('Cell.type'))
  #expect_type(result_single, "list")

  result_multiple <- get_replicate_counts_proportions_SCE(example_sce, c('Cell.type', 'Experiment'))
  #expect_true(grepl('__', colnames(result_multiple$counts)[1])) # Check separator in column names

})
library


