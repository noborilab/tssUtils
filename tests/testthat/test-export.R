test_that("writeBEDPE no longer errors on the trailing-comma data.frame call", {
  # Regression: R/export.R:20 had a trailing comma inside data.frame(), which
  # would cause an "argument is empty" runtime error before any output.
  x <- GRanges("1", IRanges(100, 200))
  y <- GRanges("1", IRanges(500, 600))
  f <- tempfile(fileext = ".bedpe")
  expect_error(writeBEDPE(x, y, f), NA)
  expect_true(file.exists(f))
})

test_that("writeBEDPE generates one name per row (regression for length(z) bug)", {
  # Regression: R/export.R:24 used length(z) on a data.frame, which returns
  # ncol (= 6), so the names vector was always 6 long regardless of nrow.
  x <- GRanges("1", IRanges(c(100, 300, 700), width = 50))
  y <- GRanges("1", IRanges(c(500, 200, 900), width = 50))
  f <- tempfile(fileext = ".bedpe")
  z <- writeBEDPE(x, y, f)
  expect_equal(nrow(z), 3)
  expect_equal(length(z$name), 3)
  expect_true(all(grepl("^INT_\\d+$", z$name)))
})

test_that("writeBEDPE writes the leftmost interval as a whole into chrom1/start1/end1", {
  x <- GRanges("1", IRanges(100, 200))
  y <- GRanges("1", IRanges(500, 600))
  f <- tempfile(fileext = ".bedpe")
  z <- writeBEDPE(x, y, f)
  expect_equal(z$start1, 100)
  expect_equal(z$end1, 200)
  expect_equal(z$start2, 500)
  expect_equal(z$end2, 600)

  # Reversed input pair: y is leftmost — it becomes interval 1.
  z2 <- writeBEDPE(y, x, tempfile(fileext = ".bedpe"))
  expect_equal(z2$start1, 100)
  expect_equal(z2$end1, 200)
  expect_equal(z2$start2, 500)
  expect_equal(z2$end2, 600)
})

test_that("writeBEDPE honours the names argument", {
  x <- GRanges("1", IRanges(100, 200))
  y <- GRanges("1", IRanges(500, 600))
  z <- writeBEDPE(x, y, tempfile(fileext = ".bedpe"), names = "custom_1")
  expect_equal(z$name, "custom_1")
})
