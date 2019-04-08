test_that("filter_genes() works faultlessly.", {
  data("input_data", package = "BANDITS")
  
  input_data = filter_genes(input_data, min_counts_per_gene = 20)
  
  expect_is(input_data, "BANDITS_data")
  expect_is(input_data@genes[[1]], "character")
  expect_is(input_data@transcripts[[1]], "character")
#  expect_is(input_data@effLen[[1]], "vector")
  expect_is(input_data@effLen[[1]], "numeric")
  expect_is(input_data@classes[[1]], "matrix")
  expect_is(input_data@counts[[1]], "data.frame")
  expect_is(input_data@uniqueId[[1]], "logical")
  expect_is(input_data@all_genes, "character")
  expect_true( length(input_data@all_genes) <= length(unlist(input_data@genes)) )
})
