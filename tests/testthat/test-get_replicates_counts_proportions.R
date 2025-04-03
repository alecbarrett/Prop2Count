
## load Seurat example

example_seurat <- readRDS(system.file("extdata", "Vignette_Seurat_object_V5.rds", package = "Prop2CountR"))



### generate dummy data

bad_proportions_matrix <- matrix(rnorm(500), ncol = 10)
colnames(bad_proportions_matrix) <- paste(rep('celltype', 10), 1:10, sep = '_')
rownames(bad_proportions_matrix) <- paste(rep('gene', 50), 1:50, sep = '_')

valid_proportions_matrix <- matrix(runif(500), ncol = 10)
colnames(valid_proportions_matrix) <- paste(rep('celltype', 10), 1:10, sep = '_')
rownames(valid_proportions_matrix) <- paste(rep('gene', 50), 1:50, sep = '_')


runif(5)

valid_ncell <- rpois(10, 10)
names(valid_ncell) <- colnames(valid_proportions_matrix)

short_ncell <- rpois(5,10)
names(short_ncell) <- names(valid_ncell)[1:5]

test_that('prop2count function requires 0 to 1 values', {
  expect_error(Prop2CountR::prop2count(bad_proportions_matrix, valid_ncell),
               'proportions input needs to be in the range of 0 and 1')
})


test_that('prop2count function requires 0 to 1 values', {
  expect_error(Prop2CountR::prop2count(valid_proportions_matrix, short_ncell),
               "The number of entries in the nCell vector, and the number of samples in the proportions matrix, do not match up")
})

test_that('valid outputs', {

  res <- Prop2CountR::prop2count(valid_proportions_matrix, valid_ncell)

  expect_true(all(is.finite(res)))

})


test_that('infinite values are fixed', {
  valid_proportions_matrix_1 <- valid_proportions_matrix
  valid_proportions_matrix_1[valid_proportions_matrix_1 > 0.9] <- 1

  res <- prop2count(valid_proportions_matrix_1, valid_ncell)

  expect_true(all(is.finite(res)))

})

library(Prop2CountR)

