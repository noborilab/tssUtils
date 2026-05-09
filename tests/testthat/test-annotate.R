test_that("prepFeatureData returns the documented slot shape", {
  fd <- mkFeatureData()
  expect_named(fd, c("feat", "featTSS", "featRC", "featRCTSS",
    "featAnn", "nameCol", "typeCol", "defaultType"),
    ignore.order = TRUE)
  expect_s4_class(fd$feat, "GRanges")
  expect_s4_class(fd$featTSS, "GRanges")
  expect_equal(width(fd$featTSS)[1], 1)
})

test_that("filterByExpression matches the documented threshold logic", {
  q <- matrix(c(0.6, 0.7, 0.8,
                0.4, 0.4, 0.4,
                0.6, 0.6, 0.4), nrow = 3, byrow = TRUE,
    dimnames = list(c("a", "b", "c"), c("s1", "s2", "s3")))
  expect_equal(unname(filterByExpression(q)),
    c(TRUE, FALSE, FALSE))
  expect_equal(unname(filterByExpression(q, minSamples = 2)),
    c(TRUE, FALSE, TRUE))
  expect_equal(unname(filterByExpression(q, minCpm = 0.4)),
    c(TRUE, TRUE, TRUE))
})

test_that("annotateTSS returns the full column set on a tiny synthetic input", {
  txData <- mkTxData()
  feData <- mkFeatureData()
  TSS <- mkTSS(starts = c(1000, 8000, 15000, 1500, 1500, 15300, 50000),
    strands = c("+", "+", "+", "-", "+", "-", "+"))
  anno <- annotateTSS(TSS, txData,
    featureSets = list(TE = feData),
    featureCoding = c(TE = "teTSS"))
  expected <- c("TSS", "Chr", "Start", "End", "Strand", "TSSTypeCoding",
    "TSSTypeLocation", "FeatureID", "FeatureSymbol", "FeatureType",
    "TSSFeatOrientation", "DistanceToFeature", "IsDivergent",
    "DivergentDistance", "DivergentTSS")
  expect_true(all(expected %in% colnames(anno)))
})

test_that("annotateTSS classifies hand-placed TSSs into the expected buckets", {
  txData <- mkTxData()
  feData <- mkFeatureData()
  TSS <- mkTSS(
    starts =   c(1000,  8000,  15000, 1500,  1500,  15300, 50000),
    strands =  c("+",   "+",   "+",   "-",   "+",   "-",   "+"),
    names = c("pcProm", "ncProm", "tePromSense", "antiTx", "intraTx", "antiTE", "ig"))
  anno <- annotateTSS(TSS, txData,
    featureSets = list(TE = feData),
    featureCoding = c(TE = "teTSS"))

  expect_equal(anno["pcProm", "TSSTypeCoding"], "pcTSS")
  expect_equal(anno["pcProm", "TSSTypeLocation"], "Promoter")
  expect_equal(anno["pcProm", "TSSFeatOrientation"], "Sense")

  expect_equal(anno["ncProm", "TSSTypeCoding"], "ncTSS")
  expect_equal(anno["ncProm", "TSSTypeLocation"], "Promoter")

  expect_equal(anno["tePromSense", "TSSTypeCoding"], "teTSS")
  expect_equal(anno["tePromSense", "TSSTypeLocation"], "Promoter")

  expect_equal(anno["antiTx", "TSSTypeLocation"], "Intragenic")
  expect_equal(anno["antiTx", "TSSFeatOrientation"], "Antisense")

  expect_equal(anno["intraTx", "TSSTypeLocation"], "Intragenic")
  expect_equal(anno["intraTx", "TSSFeatOrientation"], "Sense")

  expect_equal(anno["antiTE", "TSSTypeLocation"], "Intragenic")
  expect_equal(anno["antiTE", "TSSFeatOrientation"], "Antisense")
  expect_equal(anno["antiTE", "TSSTypeCoding"], "teTSS")

  expect_equal(anno["ig", "TSSTypeLocation"], "Intergenic")
})

test_that("annotateTSS dispatches over multiple named featureSets", {
  txData <- mkTxData()
  te <- mkFeatureData()

  snoFeat <- GRanges("1", IRanges(40000, 40500), strand = "+", name = "sno1")
  snoData <- prepFeatureData(snoFeat, defaultType = "snoRNA")

  TSS <- mkTSS(starts = c(40000, 15000),
    strands = c("+", "+"),
    names = c("snoSense", "teSense"))
  anno <- annotateTSS(TSS, txData,
    featureSets = list(TE = te, sno = snoData),
    featureCoding = c(TE = "teTSS", sno = "ncTSS"))

  expect_equal(anno["snoSense", "FeatureType"], "snoRNA")
  expect_equal(anno["snoSense", "TSSTypeCoding"], "ncTSS")
  expect_equal(anno["teSense", "FeatureType"], "transposable_element")
  expect_equal(anno["teSense", "TSSTypeCoding"], "teTSS")
})
