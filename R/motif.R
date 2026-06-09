# Motif scanning around promoters. See R/shapeActivity.R for the analyses that
# consume the output of scanPromoterMotifs().

#' Normalize a motif input to a named list of universalmotif objects.
#'
#' @param pwms A `universalmotif` object, a list of them, or a universalmotif
#'   list.
#'
#' @return A list of `universalmotif` objects with unique, non-blank names.
#'
#' @keywords internal
#' @noRd
.prepMotifs <- function(pwms) {
  motifs <- universalmotif::convert_motifs(pwms)
  if (methods::is(motifs, "universalmotif")) motifs <- list(motifs)
  ids <- vapply(motifs, function(m) as.character(m["name"]), character(1))
  bad <- is.na(ids) | !nzchar(ids)
  if (any(bad)) ids[bad] <- paste0("motif", which(bad))
  if (anyDuplicated(ids)) ids <- make.unique(ids, sep = "_")
  for (k in seq_along(motifs)) motifs[[k]]["name"] <- ids[k]
  names(motifs) <- ids
  motifs
}

#' Build per-promoter motif matrices from an occurrence GRanges.
#'
#' @param occ Occurrence `GRanges` (possibly zero-row) with mcols `promoter`,
#'   `motif`, `score`.
#' @param promoterNames Character vector of all promoter names (matrix rows).
#' @param motifIds Character vector of all motif ids (matrix columns).
#'
#' @return A list with `presence`, `count`, and `bestScore` matrices.
#'
#' @keywords internal
#' @noRd
.motifMatricesFromOccurrences <- function(occ, promoterNames, motifIds) {
  nP <- length(promoterNames)
  nM <- length(motifIds)
  count <- matrix(0L, nP, nM, dimnames = list(promoterNames, motifIds))
  bestScore <- matrix(NA_real_, nP, nM, dimnames = list(promoterNames, motifIds))
  if (length(occ)) {
    pf <- factor(occ$promoter, levels = promoterNames)
    mf <- factor(occ$motif, levels = motifIds)
    ct <- table(pf, mf)
    count[] <- as.integer(ct)
    agg <- tapply(occ$score, list(pf, mf), max)
    bestScore[] <- as.numeric(agg)
  }
  presence <- matrix((count > 0L) * 1L, nP, nM,
    dimnames = list(promoterNames, motifIds))
  list(presence = presence, count = count, bestScore = bestScore)
}

#' Scan a strand-oriented window around each promoter for motif occurrences.
#'
#' For each promoter, a window is defined relative to the dominant TSS position
#' (in TSS-strand orientation), its genomic sequence is fetched, and a motif
#' library is scanned against it with [universalmotif::scan_sequences()]. The
#' function returns both the per-occurrence hits (for the position-resolved
#' meta-profile in [motifAnchoredDelta()]) and the per-promoter
#' presence/count/best-score matrices (for the regression in
#' [shapeMotifActivity()]), so a single scan feeds both analyses.
#'
#' The anchor position is best taken from [aggregateTSSShape()], e.g.
#' `dominantPos = aggregateTSSShape(shape, TSS)$thickStart`, which keeps motif
#' scanning independent of any particular signal list.
#'
#' @param TSS A `GRanges` of TSSs (with `mcols$name`).
#' @param genome A `BSgenome` object, an [Rsamtools::FaFile()], or a character
#'   path to an indexed FASTA (its `.fai` index must exist).
#' @param pwms A `universalmotif` object, a list of them, or a universalmotif
#'   list. Conversion to a log-odds PWM, both-strand scanning, and thresholding
#'   are all handled by [universalmotif::scan_sequences()].
#' @param window Integer length-2 vector `c(upstream, downstream)` of offsets
#'   relative to the dominant TSS, in TSS-strand orientation (negative =
#'   upstream). The default `c(-50L, 10L)` spans the canonical core-promoter
#'   window.
#' @param threshold,thresholdType Passed to [universalmotif::scan_sequences()]
#'   as `threshold` and `threshold.type`. The default
#'   `threshold = 1e-4`, `thresholdType = "pvalue"` is a per-motif significance
#'   cutoff; use `thresholdType = "logodds"` with `threshold = 0.8` for the
#'   "80% of maximal log-odds" behaviour.
#' @param dominantPos Optional integer vector (one per promoter) giving the
#'   genomic position that anchors each window. Defaults to `start(TSS)`.
#' @param ignoreStrand If `FALSE` (default), both genomic strands are scanned
#'   and each hit is labelled sense or antisense relative to its TSS. If `TRUE`,
#'   only the oriented sense sequence is scanned.
#'
#' @return A list of class `"promoterMotifScan"` with two elements:
#'   `occurrences`, a `GRanges` with one row per hit (genomic coordinates and
#'   strand, plus mcols `promoter`, `motif`, `score`, `offset`, `relStrand`,
#'   `tssStrand`); and `matrices`, a list of `presence`, `count`, and
#'   `bestScore` matrices of dimension `length(TSS) × nMotif`. The `offset` mcol
#'   is the signed distance, in TSS orientation (negative = upstream), of the
#'   hit's leftmost base in the oriented window from the dominant TSS; for an
#'   antisense hit that leftmost base is the motif's 3' end on its own strand.
#'
#' @seealso [shapeMotifActivity()], [motifAnchoredDelta()],
#'   [aggregateTSSShape()].
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
scanPromoterMotifs <- function(TSS, genome, pwms,
  window = c(-50L, 10L), threshold = 1e-4, thresholdType = "pvalue",
  dominantPos = NULL, ignoreStrand = FALSE) {

  if (!requireNamespace("Biostrings", quietly = TRUE))
    stop("Package 'Biostrings' is required for scanPromoterMotifs(); install it.",
      call. = FALSE)
  if (!requireNamespace("universalmotif", quietly = TRUE))
    stop("Package 'universalmotif' is required for scanPromoterMotifs(); install it.",
      call. = FALSE)

  if (is.null(mcols(TSS)$name)) {
    TSS$name <- paste0("TSS_", seq_along(TSS))
  }
  names(TSS) <- TSS$name
  if (anyDuplicated(TSS$name))
    warning("TSS$name has duplicated values; per-promoter matrices will collapse them")

  if (length(window) != 2L)
    stop("window must be a length-2 integer vector c(upstream, downstream)")
  window <- as.integer(window)
  # window holds TSS-oriented offsets (negative = upstream), so for the
  # strand-aware promoters() call below we need upstream = -window[1] >= 0 and
  # downstream = window[2] + 1 >= 0, and a non-empty span.
  if (window[1] > window[2])
    stop("window must be c(upstream, downstream) with upstream <= downstream")
  if (window[1] > 0L || window[2] < -1L)
    stop("window offsets are relative to the TSS (negative = upstream); ",
      "window[1] must be <= 0 and window[2] >= -1")

  nT <- length(TSS)
  anchor <- if (is.null(dominantPos)) start(TSS) else as.integer(dominantPos)
  if (length(anchor) != nT)
    stop("dominantPos must have one entry per TSS")

  # Resolve a character path into an FaFile.
  if (is.character(genome)) {
    if (!requireNamespace("Rsamtools", quietly = TRUE))
      stop("Package 'Rsamtools' is required to read a FASTA path; install it, ",
        "or pass a BSgenome / FaFile object.", call. = FALSE)
    fai <- paste0(genome, ".fai")
    if (!file.exists(fai))
      stop("FASTA index not found: ", fai,
        " (build one with Rsamtools::indexFa())", call. = FALSE)
    genome <- Rsamtools::FaFile(genome)
  }

  # Build the per-promoter anchor ranges, then strand-aware windows. Drop any
  # NA anchors first, since IRanges() cannot hold NA positions.
  keepIdx <- which(!is.na(anchor))
  nDroppedNA <- nT - length(keepIdx)

  anchors <- GRanges(seqnames(TSS)[keepIdx],
    IRanges(anchor[keepIdx], anchor[keepIdx]), strand = strand(TSS)[keepIdx])
  mcols(anchors)$idx <- keepIdx

  win <- promoters(anchors, upstream = -window[1], downstream = window[2] + 1L)

  # Drop windows that fall off a chromosome end, when seqlengths are known.
  sl <- tryCatch(GenomeInfoDb::seqlengths(genome), error = function(e) NULL)
  nDroppedOOB <- 0L
  if (!is.null(sl) && length(sl)) {
    lenForWin <- sl[as.character(seqnames(win))]
    inBounds <- start(win) >= 1L &
      (is.na(lenForWin) | end(win) <= lenForWin)
    nDroppedOOB <- sum(!inBounds)
    win <- win[inBounds]
  }

  motifs <- .prepMotifs(pwms)
  motifIds <- names(motifs)

  # Motifs longer than the window cannot match; keep them as (zero) matrix
  # columns but exclude them from the scan.
  winWidth <- window[2] - window[1] + 1L
  motLen <- vapply(motifs, function(m) ncol(m["motif"]), integer(1))
  tooLong <- motLen > winWidth
  if (any(tooLong))
    message("scanPromoterMotifs(): ", sum(tooLong),
      " motif(s) longer than the ", winWidth,
      "-bp window; reported with zero occurrences")
  scanMotifs <- motifs[!tooLong]

  emptyOcc <- function() {
    gr <- GRanges(seqinfo = seqinfo(TSS))
    mcols(gr)$promoter <- character(0)
    mcols(gr)$motif <- character(0)
    mcols(gr)$score <- numeric(0)
    mcols(gr)$offset <- integer(0)
    mcols(gr)$relStrand <- character(0)
    mcols(gr)$tssStrand <- character(0)
    gr
  }

  if (nDroppedNA || nDroppedOOB)
    message("scanPromoterMotifs(): scanning ", length(win), "/", nT,
      " promoters (", nDroppedNA, " with NA anchor, ",
      nDroppedOOB, " off chromosome end)")

  if (!length(win) || !length(scanMotifs)) {
    out <- list(occurrences = emptyOcc(),
      matrices = .motifMatricesFromOccurrences(emptyOcc(), TSS$name, motifIds))
    class(out) <- "promoterMotifScan"
    return(out)
  }

  # Fetch oriented sequences; names are stable integer ids into `win`.
  seqs <- tryCatch(
    Biostrings::getSeq(genome, win),
    error = function(e)
      stop("Could not fetch sequence for the promoter windows: ", conditionMessage(e),
        "\nThis is often a seqlevels-style mismatch between TSS (",
        paste(utils::head(unique(as.character(seqnames(win))), 3), collapse = ", "),
        ") and the genome; align them with GenomeInfoDb::seqlevelsStyle().",
        call. = FALSE))
  names(seqs) <- as.character(seq_along(win))

  hits <- universalmotif::scan_sequences(scanMotifs, seqs,
    threshold = threshold, threshold.type = thresholdType,
    RC = !isTRUE(ignoreStrand), return.granges = TRUE)

  if (!length(hits)) {
    out <- list(occurrences = emptyOcc(),
      matrices = .motifMatricesFromOccurrences(emptyOcc(), TSS$name, motifIds))
    class(out) <- "promoterMotifScan"
    return(out)
  }

  occ <- .hitsToOccurrences(hits, win, TSS, window)

  out <- list(occurrences = occ,
    matrices = .motifMatricesFromOccurrences(occ, TSS$name, motifIds))
  class(out) <- "promoterMotifScan"
  out
}

#' Map oriented-frame scan hits back to genomic occurrences.
#'
#' @param hits The `GRanges` returned by [universalmotif::scan_sequences()] with
#'   `return.granges = TRUE`; seqnames are integer ids into `win`.
#' @param win The genomic promoter windows `GRanges`.
#' @param TSS The original TSS `GRanges`.
#' @param window The integer length-2 window offsets.
#'
#' @return A genomic occurrence `GRanges`.
#'
#' @keywords internal
#' @noRd
.hitsToOccurrences <- function(hits, win, TSS, window) {
  winIdx <- as.integer(as.character(seqnames(hits)))
  origIdx <- mcols(win)$idx[winIdx]
  tssStrand <- as.character(strand(win))[winIdx]

  s <- start(hits)
  e <- end(hits)
  hitStrand <- as.character(strand(hits))
  relStrand <- ifelse(hitStrand == "-", "antisense", "sense")

  winStart <- start(win)[winIdx]
  winEnd <- end(win)[winIdx]

  isMinus <- tssStrand == "-"
  gStart <- ifelse(isMinus, winEnd - e + 1L, winStart + s - 1L)
  gEnd <- ifelse(isMinus, winEnd - s + 1L, winStart + e - 1L)

  flip <- structure(c("-", "+", "*"), names = c("+", "-", "*"))
  gStrand <- ifelse(relStrand == "sense", tssStrand, flip[tssStrand])

  offset <- window[1] + (s - 1L)

  occ <- GRanges(seqnames(TSS)[origIdx],
    IRanges(gStart, gEnd), strand = gStrand,
    seqinfo = seqinfo(TSS))
  occ$promoter <- TSS$name[origIdx]
  occ$motif <- as.character(mcols(hits)$motif)
  occ$score <- as.numeric(mcols(hits)$score)
  occ$offset <- as.integer(offset)
  occ$relStrand <- relStrand
  occ$tssStrand <- tssStrand
  occ
}
