test_that("filter_transcripts() works faultlessly.", {
  gene = rep(LETTERS[1:2], each = 3)
  transcript = LETTERS[1:6]
  gene_tr_id = data.frame(gene = gene, transcript = transcript)
  
  # 2 samples, 6 transcripts
  counts <- matrix(1:6, nrow=6, ncol=2, byrow = FALSE)
  rownames(counts) = transcript
  
  transcripts_to_keep = filter_transcripts(gene_to_transcript = gene_tr_id,
                                           transcript_counts = counts, 
                                           min_transcript_proportion = 0.01,
                                           min_transcript_counts = 10, 
                                           min_gene_counts = 20)
  # Gene A has 12 counts, transcritps A, B and C will be removed (min_gene_counts = 20)
  # transcript D has 8 counts, it will be removed (min_transcript_counts = 10)
  # only transcripts E and F respect the requirements

#  expect_is(transcripts_to_keep, "vector")
  expect_is(transcripts_to_keep, "character")
  expect_true(all(transcripts_to_keep %in% rownames(counts)))
  expect_true(all(transcripts_to_keep %in% gene_tr_id[,2]))
  expect_true(length(transcripts_to_keep) <= nrow(counts) )
  expect_true(length(transcripts_to_keep) <= nrow(gene_tr_id) )
})
