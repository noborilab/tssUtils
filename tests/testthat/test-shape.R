test_that("shannonEntropy returns NA when all scores fall below the cutoff", {
  expect_true(is.na(shannonEntropy(c(0, 0, 0), minScore = 0.5)))
})

test_that("shannonEntropy is maximal for a uniform distribution", {
  uniformH <- shannonEntropy(rep(1, 4))
  spikeH <- shannonEntropy(c(10, 0.0001, 0.0001, 0.0001))
  expect_gt(uniformH, spikeH)
  expect_equal(uniformH, log(4))
})

test_that("shannonEntropy is 0 for a single-position score", {
  expect_equal(shannonEntropy(c(5)), 0)
})

test_that("simpsonDiversity is 0 for a single-position score and >0 for spread", {
  expect_equal(simpsonDiversity(c(5)), 0)
  expect_gt(simpsonDiversity(rep(1, 4)), 0)
})

test_that("tssShape returns matrices of (length(TSS) x length(signalList)) shape", {
  TSS <- mkTSS(starts = 100, strands = "+", names = "t1")
  end(TSS) <- 200
  sig <- mkSignalList(samples = c("s1", "s2"))
  out <- tssShape(TSS, sig, percentiles = c(0.1, 0.5, 0.9), minScore = 0)
  expect_named(out, c("maxPos", "maxScore", "shannon", "simpson",
    "pct1", "pct2", "pct3"), ignore.order = TRUE)
  for (nm in names(out)) {
    expect_equal(dim(out[[nm]]), c(1, 2))
  }
  expect_equal(out$maxPos[1, 1], 150)
})

test_that("aggregateTSSShape picks max-score sample for method=max", {
  shape <- list(
    maxPos = matrix(c(100, 150), nrow = 1, dimnames = list("t1", c("s1", "s2"))),
    maxScore = matrix(c(2, 5), nrow = 1, dimnames = list("t1", c("s1", "s2")))
  )
  TSS <- mkTSS(starts = 100, strands = "+", names = "t1")
  end(TSS) <- 200
  out <- aggregateTSSShape(shape, TSS, method = "max")
  expect_equal(out$thickStart, 150L)
})

test_that("aggregateTSSShape returns the median peak position for method=median", {
  shape <- list(
    maxPos = matrix(c(100, 150, 200), nrow = 1,
      dimnames = list("t1", c("s1", "s2", "s3"))),
    maxScore = matrix(c(1, 1, 1), nrow = 1,
      dimnames = list("t1", c("s1", "s2", "s3")))
  )
  TSS <- mkTSS(starts = 100, strands = "+", names = "t1")
  end(TSS) <- 200
  out <- aggregateTSSShape(shape, TSS, method = "median")
  expect_equal(out$thickStart, 150L)
})

test_that("tssThickBED sets coordinates and thick mcols on the result", {
  TSS <- mkTSS(starts = 100, strands = "+", names = "t1")
  end(TSS) <- 200
  out <- tssThickBED(TSS, start = 110, end = 190,
    thickStart = 130, thickEnd = 170)
  expect_equal(start(out), 110L)
  expect_equal(end(out), 190L)
  expect_equal(out$thickStart, 130L)
  expect_equal(out$thickEnd, 170L)
})
