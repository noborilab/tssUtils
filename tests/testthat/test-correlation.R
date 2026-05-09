test_that("correlateTSSEnhancers filters by minPCC and attaches distance", {
  TSS <- mkTSS(starts = c(1000, 2000), strands = c("+", "+"),
    names = c("t1", "t2"))
  enh <- GRanges("1", IRanges(c(5000, 8000), width = 100),
    name = c("e1", "e2"))
  names(enh) <- enh$name

  quantTSS <- matrix(c(1, 2, 3,
                       3, 2, 1), nrow = 2, byrow = TRUE,
    dimnames = list(c("t1", "t2"), c("s1", "s2", "s3")))
  quantEnh <- matrix(c(1, 2, 3,
                       3, 2, 1), nrow = 2, byrow = TRUE,
    dimnames = list(c("e1", "e2"), c("s1", "s2", "s3")))

  out <- correlateTSSEnhancers(quantTSS, quantEnh, TSS, enh, minPCC = 0.9)
  expect_true(all(out$PCC >= 0.9))
  expect_true(all(out$Dist > 0))
  expect_setequal(out$TSS, c("t1", "t2"))
})

test_that("correlateTSSEnhancers honours codingOnly via anno", {
  TSS <- mkTSS(starts = c(1000, 2000), strands = c("+", "+"),
    names = c("t1", "t2"))
  enh <- GRanges("1", IRanges(5000, width = 100), name = "e1")
  names(enh) <- enh$name
  quantTSS <- matrix(c(1, 2, 3, 3, 2, 1), nrow = 2, byrow = TRUE,
    dimnames = list(c("t1", "t2"), c("s1", "s2", "s3")))
  quantEnh <- matrix(c(1, 2, 3), nrow = 1,
    dimnames = list("e1", c("s1", "s2", "s3")))
  anno <- data.frame(row.names = c("t1", "t2"),
    TSS = c("t1", "t2"),
    TSSTypeCoding = c("pcTSS", "ncTSS"),
    stringsAsFactors = FALSE)
  out <- correlateTSSEnhancers(quantTSS, quantEnh, TSS, enh,
    minPCC = -1, anno = anno, codingOnly = TRUE)
  expect_equal(out$TSS, "t1")
})
