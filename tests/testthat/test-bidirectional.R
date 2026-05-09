test_that("findDivergent flags pcTSS/ncTSS pairs within maxDist", {
  pcTSS <- GRanges("1", IRanges(1000, 1010), strand = "+", name = "pc1")
  ncTSS <- GRanges("1", IRanges(900, 910), strand = "-", name = "nc1")
  div <- findDivergent(pcTSS, ncTSS, maxDist = 200, returnMerged = FALSE)
  expect_equal(length(div$pcTSS), 1)
  expect_equal(length(div$ncTSS), 1)
  expect_equal(div$pcTSS$DivTSS, "nc1")
  expect_equal(div$ncTSS$DivTSS, "pc1")
})

test_that("findDivergent drops pairs further than maxDist", {
  pcTSS <- GRanges("1", IRanges(1000, 1010), strand = "+", name = "pc1")
  ncTSS <- GRanges("1", IRanges(100, 110), strand = "-", name = "nc1")
  div <- findDivergent(pcTSS, ncTSS, maxDist = 200, returnMerged = FALSE)
  expect_equal(length(div$pcTSS), 0)
  expect_equal(length(div$ncTSS), 0)
})

test_that("findBidirectional collapses opposite-strand TSS pairs", {
  TSS <- GRanges("1", IRanges(c(1000, 900), c(1010, 910)),
    strand = c("+", "-"), name = c("a", "b"))
  out <- findBidirectional(TSS, maxDist = 200)
  expect_s4_class(out, "GRanges")
  expect_true(length(out) >= 1)
  expect_true(all(as.character(strand(out)) == "*"))
})
