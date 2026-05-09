test_that("swapStrand flips +/- and is involutive", {
  gr <- GRanges("1", IRanges(c(10, 50, 100), width = 5),
    strand = c("+", "-", "*"))
  out <- swapStrand(gr)
  expect_equal(as.character(strand(sort(out))),
    as.character(strand(sort(GRanges("1",
      IRanges(c(10, 50, 100), width = 5),
      strand = c("-", "+", "*"))))))
  expect_equal(sort(swapStrand(swapStrand(gr))), sort(gr))
})

test_that("calcPctiles returns coordinates inside range", {
  out <- calcPctiles(rangeStart = 100, rangeEnd = 200,
    targetPos = c(110, 150, 180),
    targetScores = c(1, 2, 1),
    pctiles = c(0.25, 0.5, 0.9))
  expect_true(all(out >= 100 & out <= 200))
  expect_true(all(diff(out) >= 0))
})
