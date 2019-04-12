test_that("test_DTU() works faultlessly.", {
  data("gene_tr_id", package = "BANDITS")
  data("precision", package = "BANDITS")
  data("input_data", package = "BANDITS")
  
  input_data@genes = input_data@genes[1:2]
  input_data@transcripts = input_data@transcripts[1:2]
  input_data@effLen = input_data@effLen[1:2]
  input_data@classes = input_data@classes[1:2]
  input_data@counts = input_data@counts[1:2]
  input_data@uniqueId = input_data@uniqueId[1:2]
  input_data@all_genes = unlist(input_data@genes[1:2])

  sample_names = paste0("sample", seq_len(4))
  samples_design = data.frame(sample_id = sample_names,
                              group = c("A", "A", "B", "B"))

  set.seed(61217)
  results = test_DTU(BANDITS_data = input_data,
                     precision = precision$prior,
                     samples_design = samples_design,
                     group_col_name = "group",
                     R = 10^4, burn_in = 2*10^3, n_cores = 2,
                     gene_to_transcript = gene_tr_id)
  
  expect_is(results, "BANDITS_test")
  expect_is(Gene_results(results), "data.frame")
  expect_is(Transcript_results(results), "data.frame")
  expect_is(convergence(results), "data.frame")
})
