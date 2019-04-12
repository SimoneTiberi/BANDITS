test_that("filter_genes() works faultlessly.", {
  data("input_data", package = "BANDITS")
  
  input_data = filter_genes(input_data, min_counts_per_gene = 20)
  
  expect_is(input_data, "BANDITS_data")
  expect_is(genes(input_data)[[1]], "character")
  expect_is(transcripts(input_data)[[1]], "character")
  expect_is(effLen(input_data)[[1]], "numeric")
  expect_is(classes(input_data)[[1]], "matrix")
  expect_is(counts(input_data)[[1]], "data.frame")
  expect_is(uniqueId(input_data)[[1]], "logical")
  expect_is(all_genes(input_data), "character")
  expect_true( all( all_genes(input_data) %in% unlist(genes(input_data)) ) )
})
