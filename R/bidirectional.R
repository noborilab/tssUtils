#' Find divergent promoters.
#'
#' @param pcTSS A `GRanges` object containing protein-coding TSSs.
#' @param ncTSS A `GRanges` object containing non-coding TSSs.
#' @param maxDist The maximum distance between antisense non-coding and
#'   protein-coding TSSs allowed for divergent promoters.
#' @param returnMerged Whether to return a set of unstranded ranges covering
#'   the divergent promoters, or a list containing separate `GRanges` for
#'   divergent protein-coding and non-coding TSSs.
#'
#' @return If `returnMerged = TRUE`, a single unstranded `GRanges` object,
#'   otherwise a list containing separate `GRanges` objects for protein-coding
#'   and non-coding TSSs.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
findDivergent <- function(pcTSS, ncTSS, maxDist = 500, returnMerged = TRUE) {
  if (is.null(mcols(ncTSS)$name)) ncTSS$name <- paste0("ncTSS_", 1:length(ncTSS))
  if (is.null(mcols(pcTSS)$name)) pcTSS$name <- paste0("pcTSS_", 1:length(pcTSS))

  names(pcTSS) <- pcTSS$name
  names(ncTSS) <- ncTSS$name

  ncTSSrc <- swapStrand(ncTSS)
  pcTSSrc <- swapStrand(pcTSS)

  ncTSSrc1 <- sort(resize(ncTSSrc, 1, "end"))
  pcTSSrc1 <- sort(resize(pcTSSrc, 1, "end"))

  ncTSS1 <- sort(resize(ncTSS, 1, "end"))
  pcTSS1 <- sort(resize(pcTSS, 1, "end"))

  pcTSS$DivTSS <- NA_character_
  pcTSS[pcTSS1$name]$DivTSS <- ncTSSrc1$name[follow(pcTSS1, ncTSSrc1)]
  pcTSS$Dist <- NA_integer_
  pcTSS$Dist[!is.na(pcTSS$DivTSS)] <- distance(pcTSS[!is.na(pcTSS$DivTSS)],
    ncTSSrc[pcTSS$DivTSS[!is.na(pcTSS$DivTSS)]])

  ncTSS$DivTSS <- NA_character_
  ncTSS[ncTSS1$name]$DivTSS <- pcTSSrc1$name[follow(ncTSS1, pcTSSrc1)]
  ncTSS$Dist <- NA_integer_
  ncTSS$Dist[!is.na(ncTSS$DivTSS)] <- distance(ncTSS[!is.na(ncTSS$DivTSS)],
    pcTSSrc[ncTSS$DivTSS[!is.na(ncTSS$DivTSS)]])

  pcTSSdiv <- pcTSS[!is.na(pcTSS$Dist) & pcTSS$Dist <= maxDist]
  ncTSSdiv <- ncTSS[!is.na(ncTSS$Dist) & ncTSS$Dist <= maxDist]

  if (returnMerged) {
    pcDiv <- promoters(resize(pcTSSdiv, 1, "end"), pcTSSdiv$Dist + width(pcTSSdiv), 1)
    ncDiv <- promoters(resize(ncTSSdiv, 1, "end"), ncTSSdiv$Dist + width(ncTSSdiv), 1)
    div <- sort(c(pcDiv, ncDiv))
    strand(div) <- "*"
    div <- sort(reduce(div))
    div$name <- paste0("DIV", as.character(seqnames(div)), "_",
      as.integer(start(div) + ((end(div) - start(div)) / 2)))
  } else {
    div <- list(pcTSS = pcTSSdiv, ncTSS = ncTSSdiv)
  }
  div
}

#' Find bidirectional promoters.
#'
#' @param TSS A `GRanges` object of TSSs.
#' @param maxDist The maximum distance between antisense TSSs to consider
#'   them bidirectional.
#'
#' @return An unstranded `GRanges` object representing bidirectional promoters.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
findBidirectional <- function(TSS, maxDist = 500) {
  if (is.null(mcols(TSS)$name)) TSS$name <- paste0("TSS_", seq_along(TSS))
  names(TSS) <- TSS$name

  TSS1 <- sort(resize(TSS, 1, "end"))
  TSSrc <- swapStrand(TSS)
  TSSrc1 <- sort(resize(TSSrc, 1, "end"))

  TSS$bdTSS <- NA_character_
  TSS[TSS1$name]$bdTSS <- TSSrc1$name[follow(TSS1, TSSrc1)]
  TSS$Dist <- NA_integer_
  TSS$Dist[!is.na(TSS$bdTSS)] <- distance(TSS[!is.na(TSS$bdTSS)], TSSrc[TSS$bdTSS[!is.na(TSS$bdTSS)]])

  bdTSS <- TSS[!is.na(TSS$Dist) & TSS$Dist <= maxDist]
  bdTSS <- resize(bdTSS, 1, fix = "end")
  bdTSS$Width <- distance(bdTSS, TSS1[bdTSS$bdTSS], ignore.strand = TRUE)

  bdGR <- promoters(bdTSS, bdTSS$Width, 1)
  strand(bdGR) <- "*"
  bdGR <- sort(reduce(sort(bdGR)))
  bdGR$name <- paste0("BD", as.character(seqnames(bdGR)), "_",
    as.integer(start(bdGR) + ((end(bdGR) - start(bdGR)) / 2)))

  bdGR
}

