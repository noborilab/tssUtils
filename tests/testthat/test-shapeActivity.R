# --- tssShapeMatched ---------------------------------------------------------

mkRegionTSS <- function() {
  TSS <- GRanges("1", IRanges(c(100, 400), c(200, 500)), strand = c("+", "-"),
    name = c("A", "B"))
  GenomeInfoDb::seqinfo(TSS) <- GenomeInfoDb::Seqinfo("1", 5000)
  TSS
}

test_that("tssShapeMatched at full depth matches tssShape on retained promoters", {
  TSS <- mkRegionTSS()
  mk <- function() {
    pos <- c(100:200, 400:500)
    sc <- c(round(10 * pmax(50 - abs(150 - (100:200)), 0)) + 1L,
            round(10 * pmax(30 - abs(450 - (400:500)), 0)) + 1L)
    st <- c(rep("+", 101), rep("-", 101))
    gr <- GRanges("1", IRanges(pos, pos), strand = st)
    gr$score <- as.integer(sc)
    GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::Seqinfo("1", 5000)
    sort(gr)
  }
  csl <- list(s1 = mk(), s2 = mk())
  ref <- tssShape(TSS, csl, minScore = 0.5)
  set.seed(1)
  m <- suppressMessages(
    tssShapeMatched(TSS, csl, target = 1e7, draws = 1, minCount = 1))
  expect_equal(m$maxPos, ref$maxPos)
  expect_equal(m$shannon, ref$shannon)
  expect_equal(m$pct4, ref$pct4)
})

test_that("tssShapeMatched is reproducible under a seed and adds antisense matrices", {
  TSS <- GRanges("1", IRanges(100, 200), strand = "+", name = "A")
  GenomeInfoDb::seqinfo(TSS) <- GenomeInfoDb::Seqinfo("1", 5000)
  csl <- mkCountSignalList(depths = c(deep = 20, shallow = 5), antisense = TRUE)
  set.seed(7)
  m1 <- suppressMessages(
    tssShapeMatched(TSS, csl, draws = 5, minCount = 5, antisense = TRUE))
  set.seed(7)
  m2 <- suppressMessages(
    tssShapeMatched(TSS, csl, draws = 5, minCount = 5, antisense = TRUE))
  expect_equal(m1$shannon, m2$shannon)
  expect_true(all(c("sumScore", "sumScoreAnti") %in% names(m1)))
  # Sense depth is matched to the shallow sample's total in both samples.
  expect_equal(m1$sumScore[1, "deep"], m1$sumScore[1, "shallow"])
})

test_that("tssShapeMatched drops promoters below minCount and logs it", {
  TSS <- mkRegionTSS()
  csl <- list(s1 = local({
    gr <- GRanges("1", IRanges(c(150, 450), c(150, 450)), strand = c("+", "-"))
    gr$score <- c(3L, 100L)  # promoter A has only 3 counts
    GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::Seqinfo("1", 5000)
    gr
  }))
  expect_message(
    m <- tssShapeMatched(TSS, csl, minCount = 10, draws = 1),
    "retaining 1/2")
  expect_equal(rownames(m$shannon), "B")
})

test_that("tssShapeMatched warns on non-integer (RPM-like) scores", {
  TSS <- GRanges("1", IRanges(100, 200), strand = "+", name = "A")
  GenomeInfoDb::seqinfo(TSS) <- GenomeInfoDb::Seqinfo("1", 5000)
  gr <- GRanges("1", IRanges(140:160, 140:160), strand = "+")
  gr$score <- rep(2.5, length(140:160))
  GenomeInfoDb::seqinfo(gr) <- GenomeInfoDb::Seqinfo("1", 5000)
  expect_warning(
    suppressMessages(tssShapeMatched(TSS, list(s1 = gr), minCount = 1, draws = 1)),
    "not RPM")
})

test_that("tssShapeMatched is reproducible via set.seed()", {
  TSS <- GRanges("1", IRanges(100, 200), strand = "+", name = "A")
  GenomeInfoDb::seqinfo(TSS) <- GenomeInfoDb::Seqinfo("1", 5000)
  csl <- mkCountSignalList(depths = c(deep = 20, shallow = 5))
  set.seed(123)
  m1 <- suppressMessages(tssShapeMatched(TSS, csl, draws = 4, minCount = 5))
  set.seed(123)
  m2 <- suppressMessages(tssShapeMatched(TSS, csl, draws = 4, minCount = 5))
  expect_equal(m1$shannon, m2$shannon)
})

test_that(".subsample returns counts unchanged when target exceeds the total", {
  n <- c(3L, 5L, 2L)
  expect_equal(tssUtils:::.subsample(n, size = 100), n)
  set.seed(1)
  s <- tssUtils:::.subsample(n, size = 5)
  expect_equal(sum(s), 5L)
  expect_true(all(s <= n))
})

# --- deltaShapeVector --------------------------------------------------------

test_that("deltaShapeVector computes each metric on a toy", {
  rn <- c("a", "b"); cn <- c("s1", "s2")
  sA <- list(shannon = matrix(c(1, 2, 3, 4), 2, dimnames = list(rn, cn)),
    maxPos = matrix(c(10, 20, 30, 40), 2, dimnames = list(rn, cn)),
    pct1 = matrix(c(5, 5, 5, 5), 2, dimnames = list(rn, cn)),
    pct2 = matrix(c(15, 25, 15, 25), 2, dimnames = list(rn, cn)),
    sumScore = matrix(c(8, 8, 8, 8), 2, dimnames = list(rn, cn)),
    sumScoreAnti = matrix(c(2, 2, 2, 2), 2, dimnames = list(rn, cn)))
  sB <- sA
  sB$shannon <- sB$shannon + matrix(c(1, 0, 1, 0), 2, dimnames = list(rn, cn))
  d <- deltaShapeVector(sA, sB, metric = "entropy")
  expect_equal(d[["a"]], 1)
  expect_equal(d[["b"]], 0)
  w <- deltaShapeVector(sA, sA, metric = "width", pctLo = 1, pctHi = 2)
  expect_equal(unname(w), c(0, 0))
  dir <- deltaShapeVector(sA, sA, metric = "directionality")
  expect_equal(unname(dir), c(0, 0))
})

test_that("deltaShapeVector enforces required arguments", {
  sh <- list(shannon = matrix(1, 1, 1, dimnames = list("a", "s1")))
  expect_error(deltaShapeVector(sh, sh, metric = "width"), "pctLo")
  expect_error(deltaShapeVector(sh, sh, metric = "directionality"), "antisense")
})

# --- shapeMotifActivity ------------------------------------------------------

test_that("shapeMotifActivity recovers planted motif coefficients", {
  skip_if_not_installed("glmnet")
  set.seed(123)
  nP <- 200; prom <- paste0("p", seq_len(nP))
  M <- matrix(rbinom(nP * 6, 1, 0.4), nP, 6,
    dimnames = list(prom, paste0("m", 1:6)))
  y <- 2 * M[, "m1"] - 1.5 * M[, "m3"] + rnorm(nP, 0, 0.3)
  names(y) <- prom
  set.seed(9)
  res <- suppressMessages(
    shapeMotifActivity(y, M, alpha = 0.5, significance = "permutation",
      nperm = 100))
  expect_s3_class(res, "data.frame")
  expect_equal(res$motif[1:2], c("m1", "m3"))
  expect_gt(res$coef[res$motif == "m1"], 0)
  expect_lt(res$coef[res$motif == "m3"], 0)
  expect_lt(res$pPerm[res$motif == "m1"], 0.05)
  expect_equal(attr(res, "nUsed"), nP)
})

test_that("shapeMotifActivity supports unpenalized covariates and drops bad columns", {
  skip_if_not_installed("glmnet")
  set.seed(1)
  nP <- 150; prom <- paste0("p", seq_len(nP))
  M <- cbind(
    m1 = rbinom(nP, 1, 0.5),
    const = rep(1L, nP))            # zero-variance, must be dropped
  rownames(M) <- prom
  gc <- runif(nP)
  y <- 1.5 * M[, "m1"] + 0.5 * gc + rnorm(nP, 0, 0.2)
  names(y) <- prom
  cov <- data.frame(gc = gc, row.names = prom)
  set.seed(2)
  res <- suppressMessages(
    shapeMotifActivity(y, M, covariates = cov, significance = "none"))
  expect_false("const" %in% res$motif)
  expect_true("m1" %in% res$motif)
})

# --- motifAnchoredDelta ------------------------------------------------------

mkAnchorSignal <- function(bump = FALSE, seqlen = 5000) {
  si <- GenomeInfoDb::Seqinfo("1", seqlen)
  pos <- 800:2200
  grp <- GRanges("1", IRanges(pos, pos), strand = "+"); grp$score <- rep(1L, length(pos))
  grm <- GRanges("1", IRanges(pos, pos), strand = "-"); grm$score <- rep(1L, length(pos))
  out <- c(grp, grm)
  if (bump) {
    bp <- GRanges("1", IRanges(1032, width = 1), strand = "+"); bp$score <- 50L
    bm <- GRanges("1", IRanges(1972, width = 1), strand = "-"); bm$score <- 50L
    out <- c(out, bp, bm)
  }
  GenomeInfoDb::seqinfo(out) <- si
  sort(out)
}

test_that("motifAnchoredDelta recovers a bump at the right distance and mirrors strand", {
  A <- list(s1 = mkAnchorSignal(FALSE))
  B <- list(s1 = mkAnchorSignal(TRUE))
  motif <- GRanges("1", IRanges(c(1000, 2000), width = 6), strand = c("+", "-"),
    name = "M")
  GenomeInfoDb::seqinfo(motif) <- GenomeInfoDb::Seqinfo("1", 5000)
  set.seed(1)
  res <- suppressMessages(
    motifAnchoredDelta(A, B, motif, window = c(100, 100), ci = "normal",
      control = "none"))
  expect_equal(nrow(res), 201L)
  expect_equal(res$distance[which.max(res$deltaMean)], 30L)
  expect_equal(max(res$deltaMean), 50)
  expect_equal(res$n[1], 2L)
})

test_that("motifAnchoredDelta supports CI modes and controls", {
  A <- list(s1 = mkAnchorSignal(FALSE))
  B <- list(s1 = mkAnchorSignal(TRUE))
  motif <- GRanges("1", IRanges(c(1000, 2000), width = 6), strand = c("+", "-"),
    name = "M")
  GenomeInfoDb::seqinfo(motif) <- GenomeInfoDb::Seqinfo("1", 5000)
  bg <- GRanges("1", IRanges(c(3000, 3500), width = 6), strand = c("+", "-"))
  GenomeInfoDb::seqinfo(bg) <- GenomeInfoDb::Seqinfo("1", 5000)
  set.seed(1)
  res <- suppressMessages(
    motifAnchoredDelta(A, B, motif, window = c(50, 50), ci = "bootstrap",
      nboot = 100, control = "background", backgroundGR = bg))
  expect_true(all(c("deltaLo", "deltaHi", "control") %in% names(res)))
  expect_equal(res$control[res$distance == 30], 0)
})

test_that("motifAnchoredDelta warns on unequal sample counts with sum", {
  A <- list(s1 = mkAnchorSignal(FALSE), s2 = mkAnchorSignal(FALSE))
  B <- list(s1 = mkAnchorSignal(FALSE))   # different sample count, same biology
  motif <- GRanges("1", IRanges(1000, width = 6), strand = "+")
  GenomeInfoDb::seqinfo(motif) <- GenomeInfoDb::Seqinfo("1", 5000)
  expect_warning(
    motifAnchoredDelta(A, B, motif, window = c(20, 20), ci = "none",
      control = "none", combineSamples = "sum"),
    "library count")
  # With mean, identical biology gives a flat (zero) delta and no warning.
  prof <- motifAnchoredDelta(A, B, motif, window = c(20, 20), ci = "none",
    control = "none", combineSamples = "mean")
  expect_lt(max(abs(prof$deltaMean)), 1e-8)
})

test_that("motifAnchoredDelta validates inputs", {
  A <- list(s1 = mkAnchorSignal(FALSE))
  expect_error(motifAnchoredDelta(A, A, GRanges()), "empty")
  starGR <- GRanges("1", IRanges(100, 105), strand = "*")
  expect_error(motifAnchoredDelta(A, A, starGR), "stranded")
  okGR <- GRanges("1", IRanges(100, 105), strand = "+")
  expect_error(
    motifAnchoredDelta(A, A, okGR, control = "background"),
    "backgroundGR")
})
