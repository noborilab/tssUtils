#' Swap strands in a GRanges object.
#'
#' @param x A `GRanges` object.
#'
#' @return A sorted `GRanges` object with swapped strands for each range.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
swapStrand <- function(x) {
  xRC <- x
  strand(xRC) <- structure(c("+", "-", "*"), names = c("-", "+", "*"))[as.character(strand(x))]
  sort(xRC)
}

#' Determine score percentile positions in a range.
#'
#' @param rangeStart Starting position of the target range.
#' @param rangeEnd Ending position of the target range.
#' @param targetPos Corresponding position coordinates for the input scores.
#' @param targetScores input scores from which to calculate percentile coordinates.
#' @param pctiles A numeric vector of percentiles.
#'
#' @return A numeric vector of coordinate percentiles.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
calcPctiles <- function(rangeStart, rangeEnd, targetPos, targetScores, pctiles) {
  x <- rep(0, 1 + (rangeEnd - rangeStart))
  x[1 + (targetPos - rangeStart)] <- targetScores
  x <- cumsum(x) / sum(x)
  sapply(pctiles, function(y) (rangeStart - 1) + which(x >= y)[1])
}

