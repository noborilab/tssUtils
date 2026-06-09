skip_if_not_installed("Biostrings")
skip_if_not_installed("universalmotif")
skip_if_not_installed("Rsamtools")

mkTATA <- function() universalmotif::create_motif("TATAAA", name = "TATA")

test_that("scanPromoterMotifs finds a sense motif at the right offset (+ TSS)", {
  fa <- mkFastaFixture(len = 400, plantAt = 50)  # TATAAA at 50-55
  TSS <- GRanges("1", IRanges(60, width = 1), strand = "+", name = "p")
  GenomeInfoDb::seqinfo(TSS) <- GenomeInfoDb::Seqinfo("1", 400)
  res <- scanPromoterMotifs(TSS, fa, mkTATA(), window = c(-50L, 10L),
    threshold = 0.9, thresholdType = "logodds")
  expect_s3_class(res, "promoterMotifScan")
  occ <- res$occurrences
  expect_equal(length(occ), 1L)
  expect_equal(occ$relStrand, "sense")
  expect_equal(occ$offset, -10L)            # motif start 50, anchor 60
  expect_equal(start(occ), 50L)
  expect_equal(end(occ), 55L)
  expect_equal(as.character(strand(occ)), "+")
})

test_that("scanPromoterMotifs labels an antisense hit and maps coordinates", {
  fa <- mkFastaFixture(len = 400, plantAt = 250)  # forward TATAAA at 250-255
  TSS <- GRanges("1", IRanges(240, width = 1), strand = "-", name = "p")
  GenomeInfoDb::seqinfo(TSS) <- GenomeInfoDb::Seqinfo("1", 400)
  res <- scanPromoterMotifs(TSS, fa, mkTATA(), window = c(-50L, 10L),
    threshold = 0.9, thresholdType = "logodds")
  occ <- res$occurrences
  expect_equal(length(occ), 1L)
  # Forward motif is antisense to a minus-strand TSS, on the + genomic strand.
  expect_equal(occ$relStrand, "antisense")
  expect_equal(as.character(strand(occ)), "+")
  expect_equal(start(occ), 250L)
  expect_equal(end(occ), 255L)
  expect_equal(occ$offset, -15L)            # 290 - 255 + 1 = 36; -50 + 35
})

test_that("scanPromoterMotifs matrices are consistent with occurrences", {
  fa <- mkFastaFixture(len = 400, plantAt = c(50, 250))
  TSS <- GRanges("1", IRanges(c(60, 240), width = 1), strand = c("+", "-"),
    name = c("a", "b"))
  GenomeInfoDb::seqinfo(TSS) <- GenomeInfoDb::Seqinfo("1", 400)
  res <- scanPromoterMotifs(TSS, fa, mkTATA(), window = c(-50L, 10L),
    threshold = 0.9, thresholdType = "logodds")
  mats <- res$matrices
  expect_equal(rownames(mats$count), c("a", "b"))
  expect_equal(colnames(mats$count), "TATA")
  expect_equal(unname(mats$count[, 1]),
    as.integer(table(factor(res$occurrences$promoter, levels = c("a", "b")))))
  expect_equal(mats$presence, (mats$count > 0L) * 1L, ignore_attr = TRUE)
  expect_true(all(mats$bestScore[mats$count > 0] > 0))
})

test_that("scanPromoterMotifs returns empty structure with no hits", {
  fa <- mkFastaFixture(len = 400, plantAt = 50)
  TSS <- GRanges("1", IRanges(60, width = 1), strand = "+", name = "p")
  GenomeInfoDb::seqinfo(TSS) <- GenomeInfoDb::Seqinfo("1", 400)
  # A GC-rich motif that won't match the planted TATA region.
  m <- universalmotif::create_motif("GGGCCC", name = "GC")
  res <- scanPromoterMotifs(TSS, fa, m, window = c(-50L, 10L),
    threshold = 1e-6, thresholdType = "pvalue")
  expect_equal(length(res$occurrences), 0L)
  expect_equal(dim(res$matrices$count), c(1L, 1L))
  expect_equal(res$matrices$count[1, 1], 0L)
  expect_true(is.na(res$matrices$bestScore[1, 1]))
})

test_that("scanPromoterMotifs reports motifs longer than the window", {
  fa <- mkFastaFixture(len = 400, plantAt = 50)
  TSS <- GRanges("1", IRanges(60, width = 1), strand = "+", name = "p")
  GenomeInfoDb::seqinfo(TSS) <- GenomeInfoDb::Seqinfo("1", 400)
  longMotif <- universalmotif::create_motif(paste(rep("A", 40), collapse = ""),
    name = "longA")
  expect_message(
    res <- scanPromoterMotifs(TSS, fa, longMotif, window = c(-10L, 5L),
      threshold = 0.9, thresholdType = "logodds"),
    "longer than")
  expect_equal(res$matrices$count[1, 1], 0L)
})

test_that("scanPromoterMotifs validates the window offsets", {
  fa <- mkFastaFixture(len = 400, plantAt = 50)
  TSS <- GRanges("1", IRanges(60, width = 1), strand = "+", name = "p")
  GenomeInfoDb::seqinfo(TSS) <- GenomeInfoDb::Seqinfo("1", 400)
  expect_error(scanPromoterMotifs(TSS, fa, mkTATA(), window = c(10L, -10L)),
    "upstream <= downstream")
  expect_error(scanPromoterMotifs(TSS, fa, mkTATA(), window = c(5L, 20L)),
    "window\\[1\\] must be <= 0")
})

test_that("scanPromoterMotifs drops NA dominantPos without crashing", {
  fa <- mkFastaFixture(len = 400, plantAt = 50)
  TSS <- GRanges("1", IRanges(c(60, 160), width = 1), strand = "+",
    name = c("a", "b"))
  GenomeInfoDb::seqinfo(TSS) <- GenomeInfoDb::Seqinfo("1", 400)
  expect_message(
    res <- scanPromoterMotifs(TSS, fa, mkTATA(), window = c(-50L, 10L),
      dominantPos = c(60L, NA_integer_),
      threshold = 0.9, thresholdType = "logodds"),
    "NA anchor")
  # Promoter "a" still scored, both promoters present in the matrices.
  expect_equal(rownames(res$matrices$count), c("a", "b"))
  expect_equal(res$matrices$count["a", "TATA"], 1L)
  expect_equal(res$matrices$count["b", "TATA"], 0L)
})

test_that("scanPromoterMotifs errors helpfully on a seqlevels-style mismatch", {
  fa <- mkFastaFixture(len = 400, plantAt = 50)
  TSS <- GRanges("chr1", IRanges(60, width = 1), strand = "+", name = "p")
  expect_error(
    scanPromoterMotifs(TSS, fa, mkTATA(), window = c(-50L, 10L)),
    "seqlevels")
})
