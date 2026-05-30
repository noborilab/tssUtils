#' Shannon entropy of a score vector.
#'
#' @param scores Numeric vector of non-negative scores.
#' @param minScore Scores strictly below this are dropped before computing
#'   the distribution.
#'
#' @return A single numeric. `NA` if no scores pass `minScore` or if all
#'   passing scores sum to zero.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
shannonEntropy <- function(scores, minScore = 0) {
  i <- !is.na(scores) & scores >= minScore
  if (!any(i)) return(NA_real_)
  scores <- scores[i]
  s <- sum(scores)
  if (s <= 0) return(NA_real_)
  p <- scores / s
  -sum(p * log(p))
}

#' Simpson diversity of a score vector.
#'
#' @param scores Numeric vector of non-negative scores.
#' @param minScore Scores strictly below this are dropped before computing
#'   the distribution.
#'
#' @return A single numeric in `[0, 1)`. `NA` if no scores pass `minScore` or
#'   if all passing scores sum to zero.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
simpsonDiversity <- function(scores, minScore = 0) {
  i <- !is.na(scores) & scores >= minScore
  if (!any(i)) return(NA_real_)
  scores <- scores[i]
  s <- sum(scores)
  if (s <= 0) return(NA_real_)
  p <- scores / s
  1 - sum(p^2)
}

#' Per-TSS shape statistics across one or more samples.
#'
#' For each TSS, and within each sample, this computes the position of the
#' highest-scoring base, the maximum score itself, the Shannon entropy, the
#' Simpson diversity, and the positions at whatever score percentiles you ask
#' for. The percentile coordinates come from [calcPctiles()].
#'
#' @param TSS A `GRanges` of TSSs (with `mcols$name`).
#' @param signalList Named list of stranded coverage `GRanges` (one per
#'   sample, e.g. produced by [readSignal()]).
#' @param percentiles Numeric vector of cumulative-score percentiles to
#'   report. The default mirrors the analysis pipeline.
#' @param minScore Threshold passed to `shannonEntropy()` and
#'   `simpsonDiversity()`.
#'
#' @return A list of matrices, all of dimensions `length(TSS) ×
#'   length(signalList)`: `maxPos`, `maxScore`, `shannon`, `simpson`, plus
#'   one matrix per percentile, named `pct1`, `pct2`, ... in the order of
#'   `percentiles`.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
tssShape <- function(TSS, signalList,
  percentiles = c(1e-99, 0.1, 0.2, 0.5, 0.8, 0.9, 1),
  minScore = 0.5) {

  if (is.null(mcols(TSS)$name))
    stop("TSS must have an mcol named 'name'")

  nT <- length(TSS)
  nS <- length(signalList)
  sampleNames <- names(signalList)
  if (is.null(sampleNames)) sampleNames <- paste0("sample", seq_len(nS))

  rowDimnames <- list(TSS$name, sampleNames)
  mkMat <- function() matrix(NA_real_, nrow = nT, ncol = nS, dimnames = rowDimnames)

  out <- list(
    maxPos = mkMat(),
    maxScore = mkMat(),
    shannon = mkMat(),
    simpson = mkMat()
  )
  pctNames <- paste0("pct", seq_along(percentiles))
  for (nm in pctNames) out[[nm]] <- mkMat()

  tssStarts <- start(TSS)
  tssEnds <- end(TSS)

  for (i in seq_len(nS)) {
    bw <- signalList[[i]]
    ovs <- findOverlaps(TSS, bw)
    if (!length(ovs)) next
    bwPos <- start(bw)
    bwScore <- bw$score
    qH <- queryHits(ovs)
    sH <- subjectHits(ovs)
    uH <- sort(unique(qH))

    out$maxPos[uH, i] <- tapply(sH, qH,
      function(x) bwPos[x][which.max(bwScore[x])])
    out$maxScore[uH, i] <- tapply(sH, qH,
      function(x) bwScore[x][which.max(bwScore[x])])
    out$shannon[uH, i] <- tapply(sH, qH,
      function(x) shannonEntropy(bwScore[x], minScore))
    out$simpson[uH, i] <- tapply(sH, qH,
      function(x) simpsonDiversity(bwScore[x], minScore))

    pctMat <- mapply(function(qi, sx)
        calcPctiles(tssStarts[qi], tssEnds[qi], bwPos[sx], bwScore[sx], percentiles),
      as.list(uH), unname(tapply(sH, qH, list)))
    if (is.null(dim(pctMat))) pctMat <- matrix(pctMat, nrow = length(percentiles))
    for (j in seq_along(percentiles)) {
      out[[pctNames[j]]][uH, i] <- pctMat[j, ]
    }
  }

  out
}

#' File-based wrapper around [tssShape()].
#'
#' Loads the paired stranded bigWig files for a set of samples and then runs
#' [tssShape()] over the resulting signal list, so you do not have to assemble
#' that list yourself.
#'
#' @param TSS A `GRanges` of TSSs.
#' @param bwDir Directory containing the bigWigs.
#' @param samples Character vector of sample names. For each `s`, the files
#'   `<bwDir>/<s><posSuffix>` and `<bwDir>/<s><negSuffix>` are read.
#' @param posSuffix,negSuffix File-name suffixes for positive- and
#'   negative-strand bigWigs.
#' @param ... Further arguments passed to [tssShape()].
#'
#' @return The list of matrices returned by [tssShape()].
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
tssShapeFromBigWig <- function(TSS, bwDir, samples,
  posSuffix = ".rpm.pos.bw", negSuffix = ".rpm.neg.bw", ...) {
  signalList <- lapply(samples, function(s)
    readSignal(file.path(bwDir, paste0(s, posSuffix)),
      file.path(bwDir, paste0(s, negSuffix)),
      nScoreIsNegative = FALSE))
  names(signalList) <- samples
  tssShape(TSS, signalList, ...)
}

#' Aggregate per-sample TSS shape into cross-sample peak coordinates.
#'
#' Produces a `GRanges` whose ranges match the input `TSS`, where the
#' `thickStart` and `thickEnd` mcols give either the peak position from the
#' single sample with the highest score (`method = "max"`), or the median peak
#' position taken across all the samples (`method = "median"`).
#'
#' @param shape List of matrices returned by [tssShape()] (must contain
#'   `maxPos` and `maxScore`).
#' @param TSS A `GRanges` of TSSs (rows aligned with `shape$maxPos`).
#' @param method Either `"max"` (per-TSS sample with highest score) or
#'   `"median"` (across-sample median peak position).
#'
#' @return A `GRanges` mirroring `TSS` with integer `thickStart`/`thickEnd`
#'   mcols.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
aggregateTSSShape <- function(shape, TSS, method = c("max", "median")) {
  method <- match.arg(method)
  maxPos <- as.matrix(shape$maxPos)
  if (method == "max") {
    maxScore <- as.matrix(shape$maxScore)
    pos <- vapply(seq_len(nrow(maxPos)), function(i) {
      ms <- maxScore[i, ]
      if (all(is.na(ms))) return(NA_real_)
      maxPos[i, which.max(ms)]
    }, numeric(1))
  } else {
    pos <- matrixStats::rowMedians(maxPos, na.rm = TRUE)
  }
  out <- TSS
  out$thickStart <- as.integer(round(pos))
  out$thickEnd <- out$thickStart
  out
}

#' Build a thick-coordinate `GRanges` for BED12-style export.
#'
#' A small convenience reshaper that sets the range coordinates, along with the
#' thick-start and thick-end mcols, on a `GRanges` from plain numeric vectors.
#'
#' @param TSS A `GRanges` whose ranges will be replaced.
#' @param start,end Numeric vectors of new range start and end coordinates.
#' @param thickStart,thickEnd Numeric vectors of thick-start and thick-end
#'   coordinates (BED12 convention).
#'
#' @return A `GRanges` with updated ranges and `thickStart`/`thickEnd` mcols.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
tssThickBED <- function(TSS, start, end, thickStart, thickEnd) {
  out <- TSS
  start(out) <- as.integer(start)
  end(out) <- as.integer(end)
  out$thickStart <- as.integer(thickStart)
  out$thickEnd <- as.integer(thickEnd)
  out
}
