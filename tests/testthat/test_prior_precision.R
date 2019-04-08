test_that("prior_precision() works faultlessly.", {
  gene = rep(LETTERS[1:2], each = 3)
  transcript = LETTERS[1:6]
  gene_tr_id = data.frame(gene = gene, transcript = transcript)
  
  # 2 samples, 6 transcripts
  counts <- 10 * matrix(1:6, nrow=6, ncol=2, byrow = FALSE)
  rownames(counts) = transcript
  
  precision = prior_precision(gene_to_transcript = gene_tr_id,
                              transcript_counts = counts,
                              min_transcript_proportion = 0.01,
                              min_transcript_counts = 10,
                              min_gene_counts = 20, n_cores = 2)

  expect_is(precision, "list")
#  expect_is(precision$prior, "vector")
  expect_is(precision$prior, "numeric")
#  expect_is(precision$genewise_log_precision, "vector")
  expect_is(precision$genewise_log_precision, "numeric")
  expect_true(length(precision$prior) == 2)
})
