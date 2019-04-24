test_that("eff_len_compute() works faultlessly.", {
  # 2 samples, 6 transcripts
  set.seed(61217)
  len <- jitter( matrix(1:6, nrow=6, ncol=2, byrow = FALSE) )
  rownames(len) = LETTERS[1:6]
  
  eff_len = eff_len_compute(x_eff_len = len)

  expect_is(eff_len, "numeric")
  expect_true(all(names(eff_len) %in% rownames(len)))
  expect_true(length(eff_len) == nrow(len) )
})
