#' Write a BEDPE file from two `GRanges`
#'
#' Each pair `(x[i], y[i])` becomes a single BEDPE row. Whichever of the two
#' intervals sits furthest to the left is written as `chrom1/start1/end1`, and
#' the other as `chrom2/start2/end2`.
#'
#' @param x The first `GRanges` object. It must be in the same order as the second.
#' @param y The second `GRanges` object. It must be in the same order as the first.
#' @param f Output filename.
#' @param names Names to use for ranges. If `NULL`, names of the form
#'   `"INT_<i>"` are generated.
#'
#' @return A `data.frame` representation of the BEDPE file, invisibly.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
writeBEDPE <- function(x, y, f, names = NULL) {
  stopifnot(length(x) == length(y))
  stopifnot(all(as.character(seqnames(x)) == as.character(seqnames(y))))
  xLeft <- start(x) <= start(y)
  start1 <- ifelse(xLeft, start(x), start(y))
  end1 <- ifelse(xLeft, end(x), end(y))
  start2 <- ifelse(xLeft, start(y), start(x))
  end2 <- ifelse(xLeft, end(y), end(x))
  z <- data.frame(row.names = NULL,
    chrom1 = as.character(seqnames(x)),
    start1 = start1,
    end1 = end1,
    chrom2 = as.character(seqnames(x)),
    start2 = start2,
    end2 = end2
  )
  ord <- order(z$chrom1, z$start1, z$start2)
  z <- z[ord, ]
  if (is.null(names)) {
    z$name <- paste0("INT_", seq_len(nrow(z)))
  } else {
    z$name <- names[ord]
  }
  write.table(z, f, quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  invisible(z)
}
