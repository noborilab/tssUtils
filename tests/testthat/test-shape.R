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

test_that("tssShape antisense=FALSE output is unchanged (no new matrices)", {
  TSS <- mkTSS(starts = 100, strands = "+", names = "t1")
  end(TSS) <- 200
  sig <- mkSignalList(samples = c("s1", "s2"))
  out <- tssShape(TSS, sig, percentiles = c(0.1, 0.5, 0.9), minScore = 0)
  expect_false(any(c("sumScore", "sumScoreAnti", "maxScoreAnti",
    "shannonAnti", "simpsonAnti") %in% names(out)))
})

test_that("tssShape antisense=TRUE adds opposite-strand matrices", {
  TSS <- mkTSS(starts = 100, strands = "+", names = "t1")
  end(TSS) <- 200
  # Sense reads on + strand, antisense reads on - strand within the window.
  sense <- GRanges("1", IRanges(120:180, 120:180), strand = "+")
  sense$score <- 100 - abs(150 - (120:180))
  anti <- GRanges("1", IRanges(140:160, 140:160), strand = "-")
  anti$score <- rep(2, length(140:160))
  sig <- list(s1 = sort(c(sense, anti)))
  out <- tssShape(TSS, sig, percentiles = c(0.5), minScore = 0, antisense = TRUE)
  expect_true(all(c("sumScore", "sumScoreAnti", "maxScoreAnti",
    "shannonAnti", "simpsonAnti") %in% names(out)))
  expect_equal(out$sumScore[1, 1], sum(sense$score))
  expect_equal(out$sumScoreAnti[1, 1], sum(anti$score))
  expect_equal(dim(out$sumScoreAnti), c(1, 1))
})

test_that("tssShape antisense sums are zero-filled when one strand is empty", {
  TSS <- mkTSS(starts = 100, strands = "+", names = "t1")
  end(TSS) <- 200
  # Antisense reads only (on the - strand); no sense reads in the window.
  anti <- GRanges("1", IRanges(140:160, 140:160), strand = "-")
  anti$score <- rep(3, length(140:160))
  out <- tssShape(TSS, list(s1 = anti), percentiles = 0.5, minScore = 0,
    antisense = TRUE)
  expect_equal(out$sumScore[1, 1], 0)                 # sense zero, not NA
  expect_equal(out$sumScoreAnti[1, 1], sum(anti$score))
  # Directionality is then defined and fully antisense.
  d <- (out$sumScore - out$sumScoreAnti) / (out$sumScore + out$sumScoreAnti)
  expect_equal(d[1, 1], -1)
})
