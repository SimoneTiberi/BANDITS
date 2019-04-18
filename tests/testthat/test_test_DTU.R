test_that("test_DTU() works faultlessly.", {
  data("gene_tr_id", package = "BANDITS")
  data("precision", package = "BANDITS")
  data("input_data", package = "BANDITS")
  
  SEL = 1:2
  
  input_data = new("BANDITS_data",
                     genes       = genes(input_data)[SEL], 
                     transcripts = transcripts(input_data)[SEL],
                     effLen      = effLen(input_data)[SEL],
                     classes     = classes(input_data)[SEL],
                     counts      = counts(input_data)[SEL], 
                     uniqueId    = uniqueId(input_data)[SEL],
                     all_genes   = all_genes(input_data)[ all_genes(input_data) %in% unlist(genes(input_data)[SEL]) ] )
  
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
