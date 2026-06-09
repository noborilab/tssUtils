# Shape-and-motif TF-activity inference. Depth-matching (tssShapeMatched), the
# delta-shape response vector (deltaShapeVector), the MARA-style regression
# (shapeMotifActivity), and the position-resolved meta-profile
# (motifAnchoredDelta). See R/motif.R for the motif scan that feeds these.

#' Downsample a vector of per-position counts to a target total.
#'
#' @param n Integer vector of counts at each position.
#' @param size Target total number of reads to keep.
#' @param replace If `FALSE` (default), sample without replacement
#'   (multivariate hypergeometric, the correct library-downsampling model). If
#'   `TRUE`, use the multinomial approximation.
#'
#' @return An integer vector of downsampled counts, aligned with `n`. If `size`
#'   is at least `sum(n)` (or `n` is all zero), `n` is returned unchanged.
#'
#' @keywords internal
#' @noRd
.subsample <- function(n, size, replace = FALSE) {
  total <- sum(n)
  if (size >= total || total == 0 || size <= 0) {
    if (size <= 0) return(integer(length(n)))
    return(as.integer(n))
  }
  if (isTRUE(replace)) {
    as.integer(stats::rmultinom(1, size, prob = n / total))
  } else {
    urn <- rep.int(seq_along(n), n)
    as.integer(tabulate(sample(urn, size), nbins = length(n)))
  }
}

#' Extract per-promoter positions and integer counts from one signal track.
#'
#' Mirrors the strand-aware `findOverlaps()` block in [tssShape()], so the
#' reconstructed per-draw signal lines up with what [tssShape()] would see.
#'
#' @param TSS A (possibly strand-flipped) `GRanges` of TSS regions.
#' @param pos Integer vector of signal positions (`start(bw)`).
#' @param count Integer vector of signal counts aligned with `pos`.
#'
#' @return A list, one entry per `TSS`, each either `NULL` or a list with `pos`
#'   and `n`.
#'
#' @keywords internal
#' @noRd
.extractCountsPerPromoter <- function(TSS, bw, pos, count) {
  res <- vector("list", length(TSS))
  ovs <- findOverlaps(TSS, bw)
  if (length(ovs)) {
    qH <- queryHits(ovs)
    sH <- subjectHits(ovs)
    spl <- split(sH, qH)
    qi <- as.integer(names(spl))
    for (j in seq_along(spl)) {
      idx <- spl[[j]]
      res[[qi[j]]] <- list(pos = pos[idx], n = count[idx])
    }
  }
  res
}

#' NA-robust mean of a list of equally shaped numeric matrices.
#'
#' @param shapeList A list of [tssShape()] outputs (one per draw).
#'
#' @return A single [tssShape()]-shaped list, averaged across draws. Cells that
#'   are `NA` in every draw stay `NA`.
#'
#' @keywords internal
#' @noRd
.avgShapes <- function(shapeList) {
  nm <- names(shapeList[[1]])
  template <- shapeList[[1]]
  D <- length(shapeList)
  out <- vector("list", length(nm))
  names(out) <- nm
  for (k in nm) {
    nr <- nrow(template[[k]])
    nc <- ncol(template[[k]])
    arr <- array(NA_real_, dim = c(nr, nc, D))
    for (d in seq_len(D)) arr[, , d] <- shapeList[[d]][[k]]
    cnt <- rowSums(!is.na(arr), dims = 2L)
    arr[is.na(arr)] <- 0
    avg <- rowSums(arr, dims = 2L) / cnt
    avg[!is.finite(avg)] <- NA_real_
    dim(avg) <- c(nr, nc)
    dimnames(avg) <- dimnames(template[[k]])
    out[[k]] <- avg
  }
  out
}

#' Depth-matched per-TSS shape statistics.
#'
#' Shannon entropy and interquantile width are downward-biased at low tag
#' counts, so unequal coverage between samples or conditions can fake a shape
#' change. This recomputes [tssShape()] after downsampling every sample, per
#' promoter, to a common target depth, averaging over several independent draws
#' to cut the sampling variance. The result is a drop-in replacement for
#' [tssShape()] output, restricted to the promoters that survive the minimum
#' count filter.
#'
#' This must run on raw integer 5' counts, not RPM-normalized tracks; read the
#' un-normalized count bigWigs/bedGraphs with [readSignal()] and pass the
#' resulting list as `countSignalList`.
#'
#' @param TSS A `GRanges` of TSS regions (with `mcols$name`).
#' @param countSignalList Named list of stranded count `GRanges` (one per
#'   sample), with integer `$score`.
#' @param target Either `NULL` (default), to downsample each promoter to its own
#'   minimum total across samples, or a single integer cap. The effective target
#'   per promoter is `pmin(perPromoterMin, target)`.
#' @param draws Number of independent downsampling draws to average (default
#'   10).
#' @param minCount Promoters whose effective target is below this are dropped
#'   (default 10).
#' @param percentiles,minScore Passed to [tssShape()].
#' @param replace If `FALSE` (default), downsample without replacement; if
#'   `TRUE`, use the multinomial approximation.
#' @param antisense Passed to [tssShape()]. When `TRUE`, antisense reads are
#'   downsampled by the same per-sample factor as the sense reads, preserving
#'   the sense/antisense ratio while matching sense depth across samples.
#'
#' @return A list of matrices in the same shape as [tssShape()] (plus the
#'   antisense matrices when `antisense = TRUE`), restricted to retained
#'   promoters. Averaged `maxPos` and `pctN` entries are fractional expected
#'   coordinates.
#'
#' @details This function draws random subsamples; call [set.seed()] beforehand
#'   for reproducible results.
#'
#' @seealso [tssShape()], [deltaShapeVector()].
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
tssShapeMatched <- function(TSS, countSignalList, target = NULL, draws = 10L,
  minCount = 10L, percentiles = c(1e-99, 0.1, 0.2, 0.5, 0.8, 0.9, 1),
  minScore = 0.5, replace = FALSE, antisense = FALSE) {

  if (is.null(mcols(TSS)$name)) TSS$name <- paste0("TSS_", seq_along(TSS))
  names(TSS) <- TSS$name

  nT <- length(TSS)
  nS <- length(countSignalList)
  if (!nS) stop("countSignalList is empty")
  sampleNames <- names(countSignalList)
  if (is.null(sampleNames)) sampleNames <- paste0("sample", seq_len(nS))

  # Coerce each sample to integer counts (warn once per non-integer sample),
  # then extract per-promoter sense (and antisense) positions and counts.
  flip <- structure(c("-", "+", "*"), names = c("+", "-", "*"))
  TSSanti <- TSS
  strand(TSSanti) <- flip[as.character(strand(TSS))]

  senseList <- vector("list", nS)
  antiList <- vector("list", nS)
  seqByProm <- as.character(seqnames(TSS))
  strandByProm <- as.character(strand(TSS))
  antiStrandByProm <- flip[strandByProm]

  for (i in seq_len(nS)) {
    bw <- countSignalList[[i]]
    sc <- bw$score
    cnt <- abs(sc)
    if (any(cnt != round(cnt), na.rm = TRUE))
      warning("tssShapeMatched(): sample '", sampleNames[i],
        "' has non-integer scores; it expects raw 5' counts, not RPM. ",
        "Rounding to integers.")
    cnt <- as.integer(round(cnt))
    pos <- start(bw)
    senseList[[i]] <- .extractCountsPerPromoter(TSS, bw, pos, cnt)
    if (isTRUE(antisense))
      antiList[[i]] <- .extractCountsPerPromoter(TSSanti, bw, pos, cnt)
  }

  # Per-promoter sense totals across samples; decide retention.
  senseTotal <- matrix(0L, nT, nS, dimnames = list(TSS$name, sampleNames))
  for (i in seq_len(nS)) {
    li <- senseList[[i]]
    for (p in seq_len(nT))
      if (!is.null(li[[p]])) senseTotal[p, i] <- sum(li[[p]]$n)
  }
  perPromoterMin <- matrixStats::rowMins(senseTotal)
  effTarget <- if (is.null(target)) perPromoterMin
    else pmin(perPromoterMin, as.integer(target))
  retain <- perPromoterMin >= minCount & effTarget >= minCount

  nRetain <- sum(retain)
  message("tssShapeMatched(): retaining ", nRetain, "/", nT,
    " promoters (", nT - nRetain, " dropped below minCount = ", minCount, ")")
  if (!nRetain)
    stop("No promoters retained; lower minCount or check countSignalList depth")

  retIdx <- which(retain)
  TSSret <- TSS[retIdx]

  antiTotal <- NULL
  if (isTRUE(antisense)) {
    antiTotal <- matrix(0L, nT, nS, dimnames = list(TSS$name, sampleNames))
    for (i in seq_len(nS)) {
      li <- antiList[[i]]
      for (p in seq_len(nT))
        if (!is.null(li[[p]])) antiTotal[p, i] <- sum(li[[p]]$n)
    }
  }

  # Build one promoter's downsampled contribution (sense, plus antisense scaled
  # by the same per-sample factor) for sample i.
  promoterContribution <- function(i, p) {
    pos <- integer(0); sco <- integer(0); st <- character(0)
    ns <- senseList[[i]][[p]]
    if (!is.null(ns) && length(ns$n)) {
      sub <- .subsample(ns$n, effTarget[p], replace)
      k <- sub > 0L
      if (any(k)) { pos <- ns$pos[k]; sco <- sub[k]; st <- rep(strandByProm[p], sum(k)) }
    }
    if (isTRUE(antisense)) {
      na <- antiList[[i]][[p]]
      if (!is.null(na) && length(na$n) && senseTotal[p, i] > 0L) {
        aSize <- as.integer(round(antiTotal[p, i] * effTarget[p] / senseTotal[p, i]))
        suba <- .subsample(na$n, aSize, replace)
        ka <- suba > 0L
        if (any(ka)) {
          pos <- c(pos, na$pos[ka]); sco <- c(sco, suba[ka])
          st <- c(st, rep(antiStrandByProm[p], sum(ka)))
        }
      }
    }
    if (!length(pos)) return(NULL)
    list(pos = pos, score = sco, strand = st,
      seq = rep(seqByProm[p], length(pos)))
  }

  # For each draw, reconstruct a downsampled signal list and recompute shape.
  # Contributions are gathered per promoter and concatenated once per sample to
  # avoid quadratic vector growth.
  drawShapes <- vector("list", draws)
  for (d in seq_len(draws)) {
    subSignalList <- vector("list", nS)
    names(subSignalList) <- sampleNames
    for (i in seq_len(nS)) {
      contribs <- lapply(retIdx, function(p) promoterContribution(i, p))
      contribs <- contribs[!vapply(contribs, is.null, logical(1))]
      gr <- GRanges(
        unlist(lapply(contribs, `[[`, "seq")),
        IRanges(unlist(lapply(contribs, `[[`, "pos")),
          unlist(lapply(contribs, `[[`, "pos"))),
        strand = unlist(lapply(contribs, `[[`, "strand")))
      gr$score <- unlist(lapply(contribs, `[[`, "score"))
      subSignalList[[i]] <- gr
    }
    drawShapes[[d]] <- tssShape(TSSret, subSignalList, percentiles, minScore,
      antisense = antisense)
  }

  .avgShapes(drawShapes)
}

#' Per-promoter delta-shape response vector between two conditions.
#'
#' Reduces two (depth-matched) [tssShape()] outputs to a single number per
#' promoter, suitable as the response in [shapeMotifActivity()]. Promoters are
#' aligned by name, the chosen statistic is aggregated across samples within
#' each condition, and condition B minus condition A is returned.
#'
#' @param shapeA,shapeB [tssShape()] / [tssShapeMatched()] outputs for the two
#'   conditions (A is the baseline, B the contrast).
#' @param metric Which shape statistic to difference:
#'   \describe{
#'     \item{`"entropy"`}{Shannon entropy (`shannon`).}
#'     \item{`"width"`}{Interquantile width, `pctHi - pctLo` in bp; requires
#'       `pctLo` and `pctHi`.}
#'     \item{`"modeShift"`}{Dominant-peak position (`maxPos`), genome-oriented,
#'       so on a minus-strand promoter a positive delta is an upstream shift.}
#'     \item{`"directionality"`}{`(sumScore - sumScoreAnti) /
#'       (sumScore + sumScoreAnti)`; requires shapes built with
#'       `antisense = TRUE`.}
#'   }
#' @param pctLo,pctHi Integer percentile indices (the `N` in `pctN`) bounding
#'   the width window. Required for `metric = "width"`.
#' @param aggregate How to aggregate across samples within a condition, either
#'   `"mean"` (default) or `"median"`.
#'
#' @return A named numeric vector (names are promoters present in both
#'   conditions) of `aggB - aggA`. `NA` in either condition yields `NA`.
#'
#' @seealso [tssShapeMatched()], [shapeMotifActivity()].
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
deltaShapeVector <- function(shapeA, shapeB,
  metric = c("entropy", "width", "modeShift", "directionality"),
  pctLo = NULL, pctHi = NULL, aggregate = c("mean", "median")) {

  metric <- match.arg(metric)
  aggregate <- match.arg(aggregate)

  aggRows <- function(m) {
    v <- if (aggregate == "mean") rowMeans(m, na.rm = TRUE)
      else matrixStats::rowMedians(m, na.rm = TRUE)
    names(v) <- rownames(m)
    v[is.nan(v)] <- NA_real_
    v
  }

  condVal <- function(shape) {
    switch(metric,
      entropy = aggRows(shape$shannon),
      width = {
        if (is.null(pctLo) || is.null(pctHi))
          stop("metric = 'width' requires pctLo and pctHi")
        lo <- shape[[paste0("pct", pctLo)]]
        hi <- shape[[paste0("pct", pctHi)]]
        if (is.null(lo) || is.null(hi))
          stop("pctLo/pctHi index a percentile not present in the shape object")
        aggRows(hi - lo)
      },
      modeShift = aggRows(shape$maxPos),
      directionality = {
        if (is.null(shape$sumScore) || is.null(shape$sumScoreAnti))
          stop("metric = 'directionality' requires shapes built with antisense = TRUE")
        denom <- shape$sumScore + shape$sumScoreAnti
        d <- (shape$sumScore - shape$sumScoreAnti) / denom
        d[denom == 0] <- NA_real_
        aggRows(d)
      })
  }

  a <- condVal(shapeA)
  b <- condVal(shapeB)
  common <- intersect(names(a), names(b))
  if (!length(common))
    stop("shapeA and shapeB share no promoter names")
  out <- b[common] - a[common]
  names(out) <- common
  out
}

#' MARA-style regularized regression of delta-shape on promoter motif content.
#'
#' Fits an elastic-net (glmnet) model of a per-promoter delta-shape response on
#' the per-promoter motif matrix, optionally adjusting for covariates. The
#' signed coefficient per motif is its estimated shape-activity change between
#' conditions; the regularization deconvolves collinear motifs. Significance is
#' assessed by permutation or stability selection, since glmnet does not give
#' usable p-values directly.
#'
#' @param deltaShape Named numeric response vector, e.g. from
#'   [deltaShapeVector()].
#' @param motifMatrix A promoter-by-motif matrix (rownames = promoter), such as
#'   `scanPromoterMotifs(...)$matrices$presence`. Position-binned columns are
#'   treated as independent predictors.
#' @param covariates Optional `data.frame` of per-promoter covariates (aligned
#'   by rowname), turned into a model matrix and added as predictors.
#' @param alpha glmnet elastic-net mixing parameter (default 0.5).
#' @param penalizeCovariates If `FALSE` (default), covariate columns are left
#'   unpenalized (always retained); if `TRUE`, they are penalized like motifs.
#' @param standardize Passed to glmnet (default `TRUE`).
#' @param nfolds Number of cross-validation folds for lambda selection.
#' @param lambda Which cross-validated lambda to use, `"lambda.1se"` (default,
#'   conservative) or `"lambda.min"`.
#' @param significance `"permutation"` (default), `"stability"`, or `"none"`.
#' @param nperm Number of permutations for `significance = "permutation"`.
#' @param nstab Number of subsamples for `significance = "stability"`.
#' @param stabSubsample Fraction of promoters per stability subsample.
#' @param stabThreshold Selection-frequency threshold reported for stability.
#'
#' @return A `data.frame`, one row per motif column, sorted by descending
#'   absolute coefficient, with columns `motif`, `coef`, `sign`, `pPerm`,
#'   `selectionFreq`, and `lambda`. Attributes `nUsed`, `nDropped`, and `alpha`
#'   record the fit. Covariate coefficients are not returned.
#'
#' @details The cross-validation, permutations, and stability subsamples are
#'   random; call [set.seed()] beforehand for reproducible results.
#'
#' @seealso [deltaShapeVector()], [scanPromoterMotifs()].
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
shapeMotifActivity <- function(deltaShape, motifMatrix, covariates = NULL,
  alpha = 0.5, penalizeCovariates = FALSE, standardize = TRUE, nfolds = 10,
  lambda = c("lambda.1se", "lambda.min"),
  significance = c("permutation", "stability", "none"),
  nperm = 1000, nstab = 100, stabSubsample = 0.5, stabThreshold = 0.6) {

  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("Package 'glmnet' is required for shapeMotifActivity(); install it.",
      call. = FALSE)

  lambda <- match.arg(lambda)
  significance <- match.arg(significance)

  emptyOut <- function() {
    out <- data.frame(motif = character(), coef = numeric(), sign = numeric(),
      pPerm = numeric(), selectionFreq = numeric(), lambda = numeric(),
      stringsAsFactors = FALSE)
    attr(out, "nUsed") <- 0L
    attr(out, "nDropped") <- 0L
    attr(out, "alpha") <- alpha
    out
  }

  if (is.null(names(deltaShape)))
    stop("deltaShape must be a named vector")
  if (is.null(rownames(motifMatrix)))
    stop("motifMatrix must have rownames (promoters)")

  motifMatrix <- as.matrix(motifMatrix)
  motifIds <- colnames(motifMatrix)

  # Align response, motifs, and covariates by promoter name.
  common <- intersect(names(deltaShape), rownames(motifMatrix))
  if (!is.null(covariates)) {
    if (is.null(rownames(covariates)))
      stop("covariates must have rownames (promoters)")
    common <- intersect(common, rownames(covariates))
  }
  y <- deltaShape[common]
  X <- motifMatrix[common, , drop = FALSE]
  cov <- if (!is.null(covariates)) covariates[common, , drop = FALSE] else NULL

  ok <- !is.na(y) & apply(X, 1, function(r) all(!is.na(r)))
  if (!is.null(cov))
    ok <- ok & stats::complete.cases(cov)
  nDropped <- sum(!ok)
  if (nDropped)
    message("shapeMotifActivity(): dropping ", nDropped,
      " promoters with missing response/predictor values")
  y <- y[ok]; X <- X[ok, , drop = FALSE]
  if (!is.null(cov)) cov <- cov[ok, , drop = FALSE]

  # Drop zero-variance motif columns (glmnet cannot use them).
  motVar <- matrixStats::colVars(X, na.rm = TRUE)
  zero <- !is.finite(motVar) | motVar == 0
  if (any(zero))
    message("shapeMotifActivity(): dropping ", sum(zero),
      " zero-variance motif column(s)")
  X <- X[, !zero, drop = FALSE]
  motifIds <- motifIds[!zero]

  if (length(y) < 10L || !ncol(X))
    stop("Too few usable promoters or motifs after filtering")

  # Assemble the design matrix and penalty factors.
  covMM <- NULL
  if (!is.null(cov)) {
    covMM <- stats::model.matrix(~ ., data = as.data.frame(cov))
    covMM <- covMM[, colnames(covMM) != "(Intercept)", drop = FALSE]
  }
  xMat <- if (is.null(covMM)) X else cbind(X, covMM)
  penalty <- rep(1, ncol(xMat))
  if (!is.null(covMM) && !isTRUE(penalizeCovariates))
    penalty[(ncol(X) + 1L):ncol(xMat)] <- 0
  motCols <- seq_len(ncol(X))

  cv <- glmnet::cv.glmnet(xMat, y, alpha = alpha, nfolds = nfolds,
    penalty.factor = penalty, standardize = standardize)
  lambdaVal <- cv[[lambda]]
  # Fit the observed coefficients with the SAME single-lambda estimator used in
  # the permutation/stability refits, so observed and null are comparable.
  obsFit <- glmnet::glmnet(xMat, y, alpha = alpha,
    penalty.factor = penalty, standardize = standardize, lambda = lambdaVal)
  coefs <- as.numeric(stats::coef(obsFit))[-1]  # drop intercept
  motCoef <- coefs[motCols]
  names(motCoef) <- motifIds

  pPerm <- rep(NA_real_, length(motCoef))
  selectionFreq <- rep(NA_real_, length(motCoef))
  names(pPerm) <- names(selectionFreq) <- motifIds

  if (significance == "permutation") {
    ge <- integer(length(motCoef))
    obsAbs <- abs(motCoef)
    for (b in seq_len(nperm)) {
      yp <- y[sample.int(length(y))]
      fit <- glmnet::glmnet(xMat, yp, alpha = alpha,
        penalty.factor = penalty, standardize = standardize,
        lambda = lambdaVal)
      cp <- as.numeric(stats::coef(fit))[-1][motCols]
      ge <- ge + (abs(cp) >= obsAbs)
    }
    pPerm <- (ge + 1) / (nperm + 1)
  } else if (significance == "stability") {
    sel <- integer(length(motCoef))
    nSub <- max(2L, floor(length(y) * stabSubsample))
    for (b in seq_len(nstab)) {
      idx <- sample.int(length(y), nSub)
      fit <- glmnet::glmnet(xMat[idx, , drop = FALSE], y[idx], alpha = alpha,
        penalty.factor = penalty, standardize = standardize,
        lambda = lambdaVal)
      cp <- as.numeric(stats::coef(fit))[-1][motCols]
      sel <- sel + (cp != 0)
    }
    selectionFreq <- sel / nstab
  }

  out <- data.frame(row.names = NULL,
    motif = motifIds,
    coef = motCoef,
    sign = sign(motCoef),
    pPerm = pPerm,
    selectionFreq = selectionFreq,
    lambda = lambdaVal,
    stringsAsFactors = FALSE)
  out <- out[order(-abs(out$coef)), , drop = FALSE]
  rownames(out) <- NULL
  attr(out, "nUsed") <- length(y)
  attr(out, "nDropped") <- nDropped
  attr(out, "alpha") <- alpha
  attr(out, "stabThreshold") <- stabThreshold
  out
}

#' Zero-filled stranded score matrix over fixed-width windows.
#'
#' Wraps [genomation::ScoreMatrix()] but, unlike [readWindowsStranded()], keeps
#' every window, filling windows with no recovered tags with zeros (a missing 5'
#' tag is a genuine zero for count data).
#'
#' @param sig A single-strand signal `GRanges`.
#' @param win The (single-strand) window `GRanges` subset.
#' @param scoreCol Score column name.
#' @param W Window width.
#'
#' @return A `length(win)` x `W` numeric matrix, rows aligned to `win`.
#'
#' @keywords internal
#' @noRd
.scoreMat <- function(sig, win, scoreCol, W) {
  m <- matrix(0, nrow = length(win), ncol = W)
  if (!length(win) || !length(sig)) return(m)
  sm <- tryCatch(
    suppressWarnings(ScoreMatrix(sig, win, strand.aware = TRUE, weight.col = scoreCol)),
    error = function(e) NULL)
  if (is.null(sm)) return(m)
  d <- sm@.Data
  # genomation drops zero-signal windows and names survivors by their index
  # into `win`; a single window comes back degenerate (W x 1, no rownames).
  if (length(win) == 1L) {
    v <- as.numeric(d)
    if (length(v) == W) {
      v[is.na(v) | v < 0] <- 0
      m[1, ] <- v
    }
    return(m)
  }
  d[is.na(d) | d < 0] <- 0
  rn <- suppressWarnings(as.integer(rownames(d)))
  if (all(is.na(rn))) {
    if (nrow(d) == length(win)) m[] <- d
  } else {
    ok <- !is.na(rn) & rn >= 1L & rn <= length(win)
    if (any(ok)) m[rn[ok], ] <- d[ok, , drop = FALSE]
  }
  m
}

#' Per-condition sense/antisense window matrices, combined across samples.
#'
#' @param signalList A list of stranded signal `GRanges` (one per sample).
#' @param win The stranded window `GRanges` (with full-window orientation).
#' @param scoreCol Score column name.
#' @param W Window width.
#' @param combineSamples `"sum"` or `"mean"` across samples.
#'
#' @return A list with `sense` and `anti` matrices, each `length(win)` x `W`.
#'
#' @keywords internal
#' @noRd
.conditionMatrices <- function(signalList, win, scoreCol, W, combineSamples) {
  nWin <- length(win)
  posIdx <- which(as.character(strand(win)) == "+")
  negIdx <- which(as.character(strand(win)) == "-")
  senseSum <- matrix(0, nWin, W)
  antiSum <- matrix(0, nWin, W)
  for (sig in signalList) {
    posSig <- sig[as.character(strand(sig)) == "+"]
    negSig <- sig[as.character(strand(sig)) == "-"]
    sense <- matrix(0, nWin, W)
    anti <- matrix(0, nWin, W)
    if (length(posIdx)) {
      sense[posIdx, ] <- .scoreMat(posSig, win[posIdx], scoreCol, W)
      anti[posIdx, ] <- .scoreMat(negSig, win[posIdx], scoreCol, W)
    }
    if (length(negIdx)) {
      sense[negIdx, ] <- .scoreMat(negSig, win[negIdx], scoreCol, W)
      anti[negIdx, ] <- .scoreMat(posSig, win[negIdx], scoreCol, W)
    }
    senseSum <- senseSum + sense
    antiSum <- antiSum + anti
  }
  if (combineSamples == "mean" && length(signalList)) {
    senseSum <- senseSum / length(signalList)
    antiSum <- antiSum / length(signalList)
  }
  list(sense = senseSum, anti = antiSum)
}

#' Position-resolved meta-profile of delta 5' signal around a motif.
#'
#' For the occurrences of a single motif, extracts per-position 5' signal in two
#' conditions over a fixed, strand-oriented window centred on each motif, and
#' returns condition B minus condition A aggregated across occurrences, as a
#' meta-profile of delta-signal versus distance from the motif. A localized
#' peak or dip at a stereotyped offset is the signature of a direct local
#' effect; compare it against the control band.
#'
#' Pass a `motifGR` already restricted to a single motif (subset
#' `scanPromoterMotifs(...)$occurrences`); facet over motifs by looping.
#'
#' @param countSignalListA,countSignalListB Lists of stranded count `GRanges`
#'   (one per sample) for conditions A (baseline) and B (contrast).
#' @param motifGR A single-motif occurrence `GRanges`; strand is required (no
#'   `*`).
#' @param window Integer length-2 vector `c(upstream, downstream)` of the
#'   profile half-widths, in motif-5' orientation.
#' @param agg Across-occurrence aggregation, `"mean"` (default) or `"median"`.
#' @param ci Confidence band: `"bootstrap"` (default), `"normal"`, or `"none"`.
#' @param ciLevel Confidence level (default 0.95).
#' @param nboot Bootstrap resamples for `ci = "bootstrap"`.
#' @param combineSamples How to combine samples within a condition, `"sum"`
#'   (default, raw counts) or `"mean"`. With `"sum"` the two conditions should
#'   have the same number of samples, otherwise the delta is confounded by
#'   library count; use `"mean"` for unequal replicate counts (a warning is
#'   issued in that case).
#' @param control Null band: `"shuffle"` (default; shift windows by `shuffleBy`
#'   bp, strand-aware), `"background"` (use `backgroundGR`), or `"none"`.
#' @param backgroundGR For `control = "background"`, a stranded `GRanges` of
#'   motif-free anchor sites; subsampled to the occurrence count when larger.
#' @param nShuffle Number of shuffles to average for `control = "shuffle"`.
#' @param shuffleBy Shift size in bp for `control = "shuffle"`.
#' @param scoreCol Score column name in the signal `GRanges`.
#'
#' @return A `data.frame` with `window[1] + window[2] + 1` rows and columns
#'   `distance`, `deltaMean`, `deltaLo`, `deltaHi`, `signalA`, `signalB`,
#'   `deltaAntisense`, `control`, and `n`. Attributes `window`, `nOccurrences`,
#'   and `agg` record the call.
#'
#' @details The bootstrap confidence band and the shuffle/background control
#'   draw random samples; call [set.seed()] beforehand for reproducible results.
#'
#' @seealso [scanPromoterMotifs()].
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
motifAnchoredDelta <- function(countSignalListA, countSignalListB, motifGR,
  window = c(100, 100), agg = c("mean", "median"),
  ci = c("bootstrap", "normal", "none"), ciLevel = 0.95, nboot = 1000,
  combineSamples = c("sum", "mean"),
  control = c("shuffle", "background", "none"), backgroundGR = NULL,
  nShuffle = 1, shuffleBy = 200, scoreCol = "score") {

  agg <- match.arg(agg)
  ci <- match.arg(ci)
  combineSamples <- match.arg(combineSamples)
  control <- match.arg(control)

  if (!length(motifGR)) stop("motifGR is empty")
  if (any(as.character(strand(motifGR)) == "*"))
    stop("motifGR must be stranded (no '*')")
  if (control == "background") {
    if (is.null(backgroundGR)) stop("control = 'background' requires backgroundGR")
    if (any(as.character(strand(backgroundGR)) == "*"))
      stop("backgroundGR must be stranded (no '*')")
  }

  if (combineSamples == "sum" &&
      length(countSignalListA) != length(countSignalListB))
    warning("motifAnchoredDelta(): conditions have different sample counts (",
      length(countSignalListA), " vs ", length(countSignalListB),
      ") and combineSamples = 'sum'; the delta will be confounded by library ",
      "count. Use combineSamples = 'mean'.")

  window <- as.integer(window)
  W <- window[1] + window[2] + 1L
  distance <- seq.int(-window[1], window[2])

  aggCol <- function(m) {
    if (!nrow(m)) return(rep(NA_real_, ncol(m)))
    if (agg == "mean") colMeans(m, na.rm = TRUE)
    else matrixStats::colMedians(m, na.rm = TRUE)
  }

  makeWindows <- function(gr) {
    centers <- resize(gr, 1, fix = "center")
    win <- suppressWarnings(promoters(centers,
      upstream = window[1], downstream = window[2] + 1L))
    win$order <- seq_along(win)
    inBounds <- start(win) >= 1L
    sl <- GenomeInfoDb::seqlengths(win)
    if (any(!is.na(sl))) {
      lenForWin <- sl[as.character(seqnames(win))]
      inBounds <- inBounds & (is.na(lenForWin) | end(win) <= lenForWin)
    }
    list(win = win[inBounds], nDropped = sum(!inBounds))
  }

  mw <- makeWindows(motifGR)
  win <- mw$win
  if (mw$nDropped)
    message("motifAnchoredDelta(): dropped ", mw$nDropped,
      " window(s) off chromosome ends")
  if (!length(win)) stop("No full-width windows remain")

  matA <- .conditionMatrices(countSignalListA, win, scoreCol, W, combineSamples)
  matB <- .conditionMatrices(countSignalListB, win, scoreCol, W, combineSamples)

  Dsense <- matB$sense - matA$sense
  Danti <- matB$anti - matA$anti

  deltaMean <- aggCol(Dsense)
  signalA <- aggCol(matA$sense)
  signalB <- aggCol(matB$sense)
  deltaAntisense <- aggCol(Danti)
  nWin <- nrow(Dsense)

  deltaLo <- rep(NA_real_, W)
  deltaHi <- rep(NA_real_, W)
  if (ci == "bootstrap") {
    boot <- matrix(NA_real_, nboot, W)
    for (b in seq_len(nboot)) {
      idx <- sample.int(nWin, nWin, replace = TRUE)
      boot[b, ] <- aggCol(Dsense[idx, , drop = FALSE])
    }
    a <- (1 - ciLevel) / 2
    qs <- matrixStats::colQuantiles(boot, probs = c(a, 1 - a), na.rm = TRUE)
    deltaLo <- qs[, 1]
    deltaHi <- qs[, 2]
  } else if (ci == "normal") {
    se <- matrixStats::colSds(Dsense, na.rm = TRUE) / sqrt(nWin)
    z <- stats::qnorm(1 - (1 - ciLevel) / 2)
    deltaLo <- deltaMean - z * se
    deltaHi <- deltaMean + z * se
  }

  control_v <- rep(NA_real_, W)
  if (control == "shuffle") {
    acc <- matrix(0, nShuffle, W)
    for (s in seq_len(nShuffle)) {
      by <- ifelse(as.character(strand(win)) == "-", -shuffleBy * s, shuffleBy * s)
      shifted <- suppressWarnings(GenomicRanges::shift(win, by))
      cmA <- .conditionMatrices(countSignalListA, shifted, scoreCol, W, combineSamples)
      cmB <- .conditionMatrices(countSignalListB, shifted, scoreCol, W, combineSamples)
      acc[s, ] <- aggCol(cmB$sense - cmA$sense)
    }
    control_v <- colMeans(acc, na.rm = TRUE)
  } else if (control == "background") {
    bg <- backgroundGR
    if (length(bg) > nWin) {
      bg <- bg[sample.int(length(bg), nWin)]
      message("motifAnchoredDelta(): subsampled backgroundGR to ", nWin, " sites")
    }
    bmw <- makeWindows(bg)
    if (length(bmw$win)) {
      cmA <- .conditionMatrices(countSignalListA, bmw$win, scoreCol, W, combineSamples)
      cmB <- .conditionMatrices(countSignalListB, bmw$win, scoreCol, W, combineSamples)
      control_v <- aggCol(cmB$sense - cmA$sense)
    }
  }

  out <- data.frame(row.names = NULL,
    distance = distance,
    deltaMean = deltaMean,
    deltaLo = deltaLo,
    deltaHi = deltaHi,
    signalA = signalA,
    signalB = signalB,
    deltaAntisense = deltaAntisense,
    control = control_v,
    n = nWin)
  attr(out, "window") <- window
  attr(out, "nOccurrences") <- nWin
  attr(out, "agg") <- agg
  out
}
