#' Read paired stranded or single unstranded bigWig/bedGraph files.
#'
#' @param p Filename for the bigWig or bedGraph file to read. If stranded, then
#'   this file represents the positive strand.
#' @param n Filename for the negative strand bigWig or bedGraph file to read.
#' @param nScoreIsNegative If `TRUE`, then ensure the negative strand scores
#'   are negative. If `FALSE`, ensure they are positive. If `NULL`, leave them
#'   as they are in the original file.
#' @param seqLvls If desired, ensure the final `GRanges` object uses these
#'   seqlevels.
#' @param windows If desired, only read signal overlapping ranges within this
#'   `GRanges` object.
#'
#' @return A `GRanges` object.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
readSignal <- function(p, n = NULL, nScoreIsNegative = NULL, seqLvls = NULL, windows = NULL) {
  p <- if (!is.null(windows)) import(p, selection = windows) else import(p)
  if (!is.null(n)) {
    n <- if (!is.null(windows)) import(n, selection = windows) else import(n)
    strand(p) <- "+"
    strand(n) <- "-"
  }
  if (isTRUE(nScoreIsNegative)) {
    n$score <- -abs(n$score)
  } else if (isFALSE(nScoreIsNegative)) {
    n$score <- abs(n$score)
  }
  if (!is.null(n)) {
    x <- sort(c(p, n))
  } else {
    x <- p
  }
  if (!is.null(seqLvls)) {
    seqlevels(x) <- as.character(seqLvls)
    x <- sort(x)
  }
  x
}

#' Read signal in windows from a `GRanges` object or a bedGraph/bigWig file
#'
#' @param signal Either a `GRanges` object or a filename.
#' @param windows A `GRanges` object containing equal width ranges.
#' @param scoreCol The name of the column containing scores if the input is a `GRanges` object.
#' @param trimSignal Lower and upper fractions of the score to trim.
#' @param keepTopPctile Only keep windows which are within the top percentile value.
#' @param rowNorm Normalize window scores between 0 and 1.
#' @param windowAverage Average signal values across all windows.
#' @param signalFactor Value to multiply scores by.
#' @param ... Additional arguments for [genomation::ScoreMatrix()].
#'
#' @return If `windowAverage = TRUE`, a `data.frame`, otherwise a `matrix`.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
readWindowsUnstranded <- function(signal, windows, scoreCol = "score",
  trimSignal = c(0, 1), keepTopPctile = 0, rowNorm = FALSE, windowAverage = TRUE,
  signalFactor = 1, ...) {

  sig <- ScoreMatrix(signal, windows, strand.aware = TRUE, weight.col = scoreCol, ...)

  sig <- sig@.Data

  sig <- sig * signalFactor
  sig[sig < 0 | is.na(sig)] <- 0

  if (trimSignal[1] != 0 || trimSignal[2] != 1) {
    qL <- stats::quantile(as.numeric(sig), p = trimSignal, na.rm = TRUE)
    sig[sig < qL[1]] <- qL[1]
    sig[sig > qL[2]] <- qL[2]
  }

  if (keepTopPctile > 0) {
    thresh <- stats::quantile(rowMaxs(sig, na.rm = TRUE) -
      rowMins(sig, na.rm = TRUE), p = keepTopPctile, na.rm = TRUE)
    sig <- sig[(rowMaxs(sig, na.rm = TRUE) - rowMins(sig, na.rm = TRUE)) >= thresh, ]
  }

  if (rowNorm) {
    s_i <- rowMaxs(sig, na.rm = TRUE)
    b_i <- rowMins(sig, na.rm = TRUE)
    sig[s_i == b_i, ] <- 0
    for (i in seq_len(nrow(sig))) {
      if (s_i[i] != b_i[i] && s_i[i] > 0) {
        sig[i, ] <- (sig[i, ] - b_i[i]) / (s_i[i] - b_i[i])
      }
    }
  }

  if (windowAverage) {
    sig_avg <- colMeans(sig, na.rm = TRUE)
    if (rowNorm) {
      sig_avg <- (sig_avg - min(sig_avg, na.rm = TRUE)) /
        (max(sig_avg, na.rm = TRUE) - min(sig_avg, na.rm = TRUE))
    }

    data.frame(row.names = NULL,
      Position = seq_along(sig_avg),
      Signal = sig_avg
    )
  } else {
    as.matrix(sig)
  }
}

#' Read stranded signal in windows from a `GRanges` object or paired bedGraph/bigWig files
#'
#' Given strand-aware input signal, this computes the per-window sense and
#' antisense signal matrices, and will optionally average them across windows
#' for you.
#'
#' @param signal Either a stranded `GRanges` object, or a length-2 character
#'   vector of filenames (positive strand first, then negative strand).
#' @param windows A stranded `GRanges` object containing equal width ranges.
#' @param scoreCol The name of the column containing scores if the input is a `GRanges` object.
#' @param trimSignal Lower and upper fractions of the score to trim.
#' @param keepTopPctile Only keep windows which are within the top percentile value.
#' @param rowNorm Normalize window scores between 0 and 1.
#' @param windowAverage Average signal values across all windows.
#' @param mergeSignalMatrix If `windowAverage = FALSE`, whether to return a
#'   single sense-minus-antisense matrix or a list with `Positive` and
#'   `Negative` matrices.
#' @param signalFactor Length-2 numeric vector of factors to multiply the
#'   positive- and negative-strand scores by.
#' @param negStrandValues If `TRUE`, antisense scores are returned as negative
#'   numbers.
#' @param ... Additional arguments for [genomation::ScoreMatrix()].
#'
#' @return If `windowAverage = TRUE`, a `data.frame` with columns `Position`,
#'   `SignalPositive`, `SignalNegative`, `SignalMerged`. Otherwise, a `matrix`
#'   (when `mergeSignalMatrix = TRUE`) or a list of two matrices.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
readWindowsStranded <- function(signal, windows, scoreCol = "score",
  trimSignal = c(0, 1), keepTopPctile = 0, rowNorm = FALSE, windowAverage = TRUE,
  mergeSignalMatrix = TRUE, signalFactor = c(1, 1), negStrandValues = FALSE, ...) {

  windows$order <- seq_along(windows)

  pos_win <- windows[as.character(strand(windows)) == "+"]
  neg_win <- windows[as.character(strand(windows)) == "-"]
  if (is(signal, "GRanges")) {
    pos_sig <- signal[as.character(strand(signal)) == "+"]
    neg_sig <- signal[as.character(strand(signal)) == "-"]
  } else if (is.character(signal) && length(signal) == 2) {
    pos_sig <- signal[1]
    neg_sig <- signal[2]
  } else {
    stop("signal should be a GRanges object or two filenames")
  }

  if (!length(pos_win) && !length(neg_win))
    stop("Couldn't find any windows on +/- strands")

  if (length(pos_win)) {
    pos_pos <- ScoreMatrix(pos_sig, pos_win, strand.aware = TRUE, weight.col = scoreCol, ...)
    neg_pos <- ScoreMatrix(neg_sig, pos_win, strand.aware = TRUE, weight.col = scoreCol, ...)
  } else {
    pos_pos <- matrix(nrow = 0, ncol = width(windows)[1])
    neg_pos <- matrix(nrow = 0, ncol = width(windows)[1])
  }
  if (length(neg_win)) {
    pos_neg <- ScoreMatrix(pos_sig, neg_win, strand.aware = TRUE, weight.col = scoreCol, ...)
    neg_neg <- ScoreMatrix(neg_sig, neg_win, strand.aware = TRUE, weight.col = scoreCol, ...)
  } else {
    pos_neg <- matrix(nrow = 0, ncol = width(windows)[1])
    neg_neg <- matrix(nrow = 0, ncol = width(windows)[1])
  }

  if (nrow(pos_pos) && nrow(neg_pos)) {
    good_pos <- dimnames(pos_pos)[[1]][dimnames(pos_pos)[[1]] %in% dimnames(neg_pos)[[1]]]
  } else {
    good_pos <- character()
  }
  if (nrow(pos_neg) && nrow(neg_neg)) {
    good_neg <- dimnames(pos_neg)[[1]][dimnames(pos_neg)[[1]] %in% dimnames(neg_neg)[[1]]]
  } else {
    good_neg <- character()
  }

  if (!length(good_pos) && !length(good_neg))
    stop("Could not recover any ranges with signal")

  if (length(good_pos)) {
    pos_pos <- pos_pos[good_pos, ]@.Data * signalFactor[1]
    neg_pos <- neg_pos[good_pos, ]@.Data * signalFactor[2]
    if (negStrandValues) neg_pos <- -neg_pos
  }
  if (length(good_neg)) {
    pos_neg <- pos_neg[good_neg, ]@.Data * signalFactor[1]
    neg_neg <- neg_neg[good_neg, ]@.Data * signalFactor[2]
    if (negStrandValues) neg_neg <- -neg_neg
  }

  pos_pos[pos_pos < 0 | is.na(pos_pos)] <- 0
  pos_neg[pos_neg < 0 | is.na(pos_neg)] <- 0
  neg_pos[neg_pos < 0 | is.na(neg_pos)] <- 0
  neg_neg[neg_neg < 0 | is.na(neg_neg)] <- 0

  if (trimSignal[1] != 0 || trimSignal[2] != 1) {
    qL <- stats::quantile(c(pos_pos, pos_neg, neg_pos, neg_neg), p = trimSignal, na.rm = TRUE)
    pos_pos[pos_pos < qL[1]] <- qL[1]
    pos_neg[pos_neg < qL[1]] <- qL[1]
    neg_pos[neg_pos < qL[1]] <- qL[1]
    neg_neg[neg_neg < qL[1]] <- qL[1]
    pos_pos[pos_pos > qL[2]] <- qL[2]
    pos_neg[pos_neg > qL[2]] <- qL[2]
    neg_pos[neg_pos > qL[2]] <- qL[2]
    neg_neg[neg_neg > qL[2]] <- qL[2]
  }

  if (keepTopPctile > 0) {
    thresh_sense <- stats::quantile(
      c(rowMaxs(pos_pos, na.rm = TRUE) -
        rowMins(pos_pos, na.rm = TRUE),
        rowMaxs(neg_neg, na.rm = TRUE) -
          rowMins(neg_neg, na.rm = TRUE)),
      p = keepTopPctile, na.rm = TRUE)
    thresh_antisense <- stats::quantile(
      c(rowMaxs(neg_pos, na.rm = TRUE) -
        rowMins(neg_pos, na.rm = TRUE),
        rowMaxs(pos_neg, na.rm = TRUE) -
          rowMins(pos_neg, na.rm = TRUE)),
      p = keepTopPctile, na.rm = TRUE)
    pos_pos <- pos_pos[(rowMaxs(pos_pos, na.rm = TRUE) -
      rowMins(pos_pos, na.rm = TRUE)) >= thresh_sense, ]
    neg_neg <- neg_neg[(rowMaxs(neg_neg, na.rm = TRUE) -
      rowMins(neg_neg, na.rm = TRUE)) >= thresh_sense, ]
    neg_pos <- neg_pos[(rowMaxs(neg_pos, na.rm = TRUE) -
      rowMins(neg_pos, na.rm = TRUE)) >= thresh_antisense, ]
    pos_neg <- pos_neg[(rowMaxs(pos_neg, na.rm = TRUE) -
      rowMins(pos_neg, na.rm = TRUE)) >= thresh_antisense, ]
  }

  if (rowNorm) {

    for (i in seq_len(nrow(pos_pos))) {
      s_i <- max(pos_pos[i, ], na.rm = TRUE)
      b_i <- min(pos_pos[i, ], na.rm = TRUE)
      if (s_i == b_i) pos_pos[i, ] <- 0
      else if (s_i > 0) pos_pos[i, ] <- (pos_pos[i, ] - b_i) / (s_i - b_i)
    }

    for (i in seq_len(nrow(neg_pos))) {
      s_i <- max(neg_pos[i, ], na.rm = TRUE)
      b_i <- min(neg_pos[i, ], na.rm = TRUE)
      if (s_i == b_i) neg_pos[i, ] <- 0
      else if (s_i > 0) neg_pos[i, ] <- (neg_pos[i, ] - b_i) / (s_i - b_i)
    }

    for (i in seq_len(nrow(pos_neg))) {
      s_i <- max(pos_neg[i, ], na.rm = TRUE)
      b_i <- min(pos_neg[i, ], na.rm = TRUE)
      if (s_i == b_i) pos_neg[i, ] <- 0
      else if (s_i > 0) pos_neg[i, ] <- (pos_neg[i, ] - b_i) / (s_i - b_i)
    }

    for (i in seq_len(nrow(neg_neg))) {
      s_i <- max(neg_neg[i, ], na.rm = TRUE)
      b_i <- min(neg_neg[i, ], na.rm = TRUE)
      if (s_i == b_i) neg_neg[i, ] <- 0
      else if (s_i > 0) neg_neg[i, ] <- (neg_neg[i, ] - b_i) / (s_i - b_i)
    }

  }

  if (windowAverage) {

    pos_avg <- colMeans(rbind(pos_pos, neg_neg), na.rm = TRUE)
    neg_avg <- colMeans(rbind(pos_neg, neg_pos), na.rm = TRUE)

    if (rowNorm) {
      if (any(pos_avg > 0))
        pos_avg <- (pos_avg - min(pos_avg, na.rm = TRUE)) /
          (max(pos_avg, na.rm = TRUE) - min(pos_avg, na.rm = TRUE))
      if (any(neg_avg > 0))
        neg_avg <- (neg_avg - min(neg_avg, na.rm = TRUE)) /
          (max(neg_avg, na.rm = TRUE) - min(neg_avg, na.rm = TRUE))
    }

    data.frame(row.names = NULL,
      Position = seq_along(pos_avg),
      SignalPositive = pos_avg,
      SignalNegative = neg_avg,
      SignalMerged = pos_avg - neg_avg
    )

  } else {

    good_pos <- pos_win$order[as.integer(good_pos)]
    good_neg <- neg_win$order[as.integer(good_neg)]

    sense_m <- rbind(pos_pos, neg_neg)[order(c(good_pos, good_neg)), ]
    antisense_m <- rbind(neg_pos, pos_neg)[order(c(good_pos, good_neg)), ]

    if (!mergeSignalMatrix) {
      list(Positive = sense_m, Negative = antisense_m)
    } else {
      sense_m - antisense_m
    }

  }
}
