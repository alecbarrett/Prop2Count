

### generate dummy data ----

bad_proportions_matrix <- matrix(rnorm(500), ncol = 10)
colnames(bad_proportions_matrix) <- paste(rep('celltype', 10), 1:10, sep = '_')
rownames(bad_proportions_matrix) <- paste(rep('gene', 50), 1:50, sep = '_')

valid_proportions_matrix <- matrix(runif(500), ncol = 10)
colnames(valid_proportions_matrix) <- paste(rep('celltype', 10), 1:10, sep = '_')
rownames(valid_proportions_matrix) <- paste(rep('gene', 50), 1:50, sep = '_')


valid_ncell <- rpois(10, 10)
names(valid_ncell) <- colnames(valid_proportions_matrix)

valid_ncell_allzero <- rep(0,10)
names(valid_ncell_allzero) <- colnames(valid_proportions_matrix)

short_ncell <- rpois(5,10)
names(short_ncell) <- names(valid_ncell)[1:5]



### run tests ----

test_that('prop2count function requires 0 to 1 values', {
  expect_error(prop2count(bad_proportions_matrix, valid_ncell),
               'proportions input needs to be in the range of 0 and 1')
})


test_that('prop2count function requires 0 to 1 values', {
  expect_error(prop2count(valid_proportions_matrix, short_ncell),
               "The number of entries in the nCell vector, and the number of samples in the proportions matrix, do not match up")
})

test_that('valid outputs', {

  res <- prop2count(valid_proportions_matrix, valid_ncell)

  expect_true(all(is.finite(res)))

})


test_that('infinite values are fixed', {
  valid_proportions_matrix_1 <- valid_proportions_matrix
  valid_proportions_matrix_1[valid_proportions_matrix_1 > 0.9] <- 1

  res <- prop2count(valid_proportions_matrix_1, valid_ncell)

  expect_true(all(is.finite(res)))

})


test_that("prop2count rounds correctly", {

  result_no_round <- prop2count(valid_proportions_matrix, valid_ncell, round = FALSE)

  result_round <- prop2count(valid_proportions_matrix, valid_ncell, round = TRUE)

  expect_equal(result_round, round(result_no_round))
  expect_true(all(result_round == floor(result_round)))
})


test_that("prop2count catches 0 cells values", {

  expect_error(prop2count(valid_proportions_matrix, valid_ncell_allzero))

})

test_that('NAs in proportions matrix throws an error', {

  valid_proportions_matrix_NA <- valid_proportions_matrix
  valid_proportions_matrix_NA[sample(1:length(valid_proportions_matrix_NA), 10)] <- NA

  expect_error(prop2count(valid_proportions_matrix_NA, valid_ncell))

})

test_that('NAs in nCell vector throws an error', {

  valid_ncell_NA <- valid_ncell
  valid_ncell_NA[sample(1:length(valid_ncell_NA), 2)] <- NA

  expect_error(prop2count(valid_proportions_matrix, valid_ncell_NA))

})


test_that('results matrix matches proportions input', {

  results <- prop2count(valid_proportions_matrix, valid_ncell)

  expect_true(all(dim(results) == dim(valid_proportions_matrix)))

})

