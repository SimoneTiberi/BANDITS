test_that("create_data() works faultlessly.", {
  data("gene_tr_id", package = "BANDITS")

  data_dir = system.file("extdata", package = "BANDITS")
  sample_name = paste0("sample", 1)
  equiv_classes_files = file.path(data_dir, "STAR-salmon", sample_name, 
                                  "aux_info", "eq_classes.txt")
  
  quant_files = file.path(data_dir, "STAR-salmon", sample_name, "quant.sf")
  
  txi = tximport::tximport(files = quant_files, type = "salmon", txOut = TRUE)
  
  eff_len = eff_len_compute(x_eff_len = txi$length)

  input_data = create_data(salmon_or_kallisto = "salmon",
                           gene_to_transcript = gene_tr_id,
                           salmon_path_to_eq_classes = equiv_classes_files,
                           eff_len = eff_len, 
                           n_cores = 1,
                           transcripts_to_keep = NULL)
  
  
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
