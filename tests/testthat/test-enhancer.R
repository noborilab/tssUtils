test_that("defineEnhancers returns the three-element list shape", {
  TSS <- mkTSS(starts = c(50000, 50100, 60000),
    strands = c("+", "-", "+"),
    names = c("a", "b", "c"))
  anno <- data.frame(row.names = TSS$name,
    TSS = TSS$name,
    Chr = "1",
    Start = start(TSS),
    End = end(TSS),
    Strand = as.character(strand(TSS)),
    TSSTypeCoding = c("ncTSS", "ncTSS", "ncTSS"),
    TSSTypeLocation = c("Intergenic", "Intergenic", "Intergenic"),
    FeatureID = NA_character_,
    FeatureSymbol = NA_character_,
    FeatureType = NA_character_,
    TSSFeatOrientation = NA_character_,
    DistanceToFeature = NA_integer_,
    IsDivergent = c(FALSE, FALSE, FALSE),
    DivergentDistance = NA_integer_,
    DivergentTSS = NA_character_,
    stringsAsFactors = FALSE)
  exprFilt <- c(TRUE, TRUE, TRUE)
  out <- defineEnhancers(TSS, anno, exprFilt)
  expect_named(out, c("bdGR", "intTSS", "enh"))
  expect_s4_class(out$bdGR, "GRanges")
  expect_s4_class(out$intTSS, "GRanges")
  expect_s4_class(out$enh, "GRanges")
})

test_that("quantifyEnhancers sums TSS quantifications into enhancer totals", {
  enh <- GRanges("1", IRanges(start = c(50000, 70000),
    end = c(50500, 70500)), name = c("E1", "E2"))
  TSS <- mkTSS(starts = c(50100, 50200, 70300),
    strands = c("+", "-", "+"),
    names = c("a", "b", "c"))
  quant <- mkQuant(c("a", "b", "c"), nSamples = 2, fill = 1)
  out <- quantifyEnhancers(enh, TSS, quant)
  expect_equal(dim(out), c(2, 2))
  expect_equal(unname(out["E1", ]), c(2, 2))
  expect_equal(unname(out["E2", ]), c(1, 1))
})

test_that("enhancerStrandStats flags bidirectional regions", {
  enh <- GRanges("1", IRanges(50000, 50500), name = "E1")
  TSS <- mkTSS(starts = c(50100, 50200),
    strands = c("+", "-"),
    names = c("p1", "n1"))
  quant <- mkQuant(c("p1", "n1"), nSamples = 3, fill = 5)
  exprFilt <- c(TRUE, TRUE)
  out <- enhancerStrandStats(enh, TSS, quant, exprFilt)
  expect_true(out$IsBidirectional[1])
  expect_equal(out$Enhancer, "E1")
})
