#' Correlate TSS expression with enhancer expression.
#'
#' For each TSS-enhancer pair on the same chromosome, computes a Pearson
#' correlation across samples on log2(1 + CPM)-transformed values. Pairs with
#' `PCC < minPCC` are dropped, as are self-pairs (`Dist == 0`). Distances are
#' computed via [GenomicRanges::distance()].
#'
#' @param quantTSS A numeric matrix or `data.frame` of TSS quantifications
#'   (rows aligned with `TSS`).
#' @param quantEnh A numeric matrix or `data.frame` of enhancer
#'   quantifications (rows aligned with `enh`).
#' @param TSS A `GRanges` of TSSs.
#' @param enh A `GRanges` of enhancers.
#' @param minPCC Minimum Pearson correlation to retain.
#' @param anno Optional annotation `data.frame` from [annotateTSS()]; required
#'   when `codingOnly = TRUE`.
#' @param codingOnly If `TRUE`, restrict to TSSs with
#'   `anno$TSSTypeCoding == "pcTSS"`.
#'
#' @return A `data.frame` with columns `TSS`, `Enhancer`, `PCC`, `Dist`.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
correlateTSSEnhancers <- function(quantTSS, quantEnh, TSS, enh, minPCC = 0.5,
  anno = NULL, codingOnly = FALSE) {

  if (codingOnly) {
    if (is.null(anno)) stop("anno is required when codingOnly = TRUE")
    keep <- anno$TSSTypeCoding == "pcTSS"
    keep[is.na(keep)] <- FALSE
    keepNames <- anno$TSS[keep]
    sel <- rownames(quantTSS) %in% keepNames
    quantTSS <- quantTSS[sel, , drop = FALSE]
    TSS <- TSS[sel]
  }

  if (!nrow(quantTSS) || !nrow(quantEnh)) {
    return(data.frame(TSS = character(), Enhancer = character(),
      PCC = numeric(), Dist = integer(), stringsAsFactors = FALSE))
  }

  m <- suppressWarnings(stats::cor(t(log2(1 + as.matrix(quantTSS))),
    t(log2(1 + as.matrix(quantEnh)))))

  keepIdx <- which(!is.na(m) & m >= minPCC, arr.ind = TRUE)
  if (!nrow(keepIdx)) {
    return(data.frame(TSS = character(), Enhancer = character(),
      PCC = numeric(), Dist = integer(), stringsAsFactors = FALSE))
  }

  out <- data.frame(row.names = NULL,
    TSS = rownames(m)[keepIdx[, 1]],
    Enhancer = colnames(m)[keepIdx[, 2]],
    PCC = m[keepIdx],
    stringsAsFactors = FALSE
  )
  names(TSS) <- TSS$name
  names(enh) <- enh$name
  out$Dist <- distance(enh[out$Enhancer], TSS[out$TSS])
  out <- out[!is.na(out$Dist) & out$Dist > 0, , drop = FALSE]
  rownames(out) <- NULL
  out
}
