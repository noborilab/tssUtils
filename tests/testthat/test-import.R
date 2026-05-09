test_that("readSignal round-trips a single bedGraph", {
  gr <- GRanges("1", IRanges(c(10, 50), width = 5), score = c(2.5, -1.5))
  f <- tempfile(fileext = ".bedGraph")
  rtracklayer::export.bedGraph(gr, f)
  out <- readSignal(f)
  expect_s4_class(out, "GRanges")
  expect_equal(sort(out$score), sort(gr$score))
})

test_that("readSignal forces negative-strand sign with nScoreIsNegative", {
  pos <- GRanges("1", IRanges(10, width = 5), score = 3)
  neg <- GRanges("1", IRanges(50, width = 5), score = 4)
  fp <- tempfile(fileext = ".bedGraph")
  fn <- tempfile(fileext = ".bedGraph")
  rtracklayer::export.bedGraph(pos, fp)
  rtracklayer::export.bedGraph(neg, fn)
  outNeg <- readSignal(fp, fn, nScoreIsNegative = TRUE)
  outPos <- readSignal(fp, fn, nScoreIsNegative = FALSE)
  expect_true(all(outNeg$score[as.character(strand(outNeg)) == "-"] <= 0))
  expect_true(all(outPos$score[as.character(strand(outPos)) == "-"] >= 0))
})

test_that("readWindowsStranded source no longer references the neg_eg typo", {
  # Regression: R/import.R:203 referenced an undefined `neg_eg` variable,
  # which would error at runtime when keepTopPctile > 0.
  src <- paste(deparse(readWindowsStranded), collapse = "\n")
  expect_false(grepl("neg_eg", src, fixed = TRUE))
  expect_true(grepl("neg_neg", src, fixed = TRUE))
})
