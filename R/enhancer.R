#' Define enhancer regions from annotated TSSs.
#'
#' Builds two kinds of enhancer-like region and hands them back together with
#' their merged union:
#'
#' - **Bidirectional ncTSS regions**: pairs of antisense ncTSSs whose ends fall
#'   within `divergentDist` bp of each other, with the span between the pair
#'   taken as the region.
#' - **Unidirectional intergenic ncTSS regions**: ncTSSs that were classified as
#'   `Intergenic` (or as one of `intergenicTypes`) and are not divergent, each
#'   extended `intergenicFlank` bp upstream of its 5' end.
#'
#' Any region that ends up overlapping a pcTSS is dropped, since it is more
#' likely to belong to a known promoter than to an enhancer.
#'
#' @param TSS A `GRanges` of called TSSs (with `mcols$name` matching `anno$TSS`).
#' @param anno Annotation `data.frame` produced by [annotateTSS()].
#' @param exprFilt Logical vector aligned with `anno`/`TSS` indicating which
#'   TSSs pass the expression filter (e.g. from [filterByExpression()]).
#' @param divergentDist Maximum distance between paired antisense ncTSSs that
#'   counts as bidirectional.
#' @param intergenicFlank Upstream extension applied to unidirectional
#'   intergenic ncTSSs.
#' @param intergenicTypes Feature types treated as intergenic-equivalent for
#'   the unidirectional rule.
#'
#' @return A list with three `GRanges` slots: `bdGR` (bidirectional regions),
#'   `intTSS` (unidirectional intergenic regions), and `enh` (the sorted union
#'   of both, with `name` mcol).
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
defineEnhancers <- function(TSS, anno, exprFilt, divergentDist = 500,
  intergenicFlank = 500,
  intergenicTypes = c("transposable_element", "transposable_element_gene",
    "antisense_lncRNA", "antisense_RNA", "lnc_RNA",
    "miRNA_primary_transcript", "ncRNA", "pseudogenic_transcript",
    "transcript_region")) {

  isPC <- anno$TSSTypeCoding == "pcTSS"
  isNC <- !isPC

  pcTSS <- TSS[isPC & exprFilt]
  ncTSS <- TSS[isNC & exprFilt]

  if (length(ncTSS) > 1) {
    ncTSSrc <- swapStrand(ncTSS)
    ncTSS1e <- sort(resize(ncTSS, 1, "end"))
    ncTSSrc1e <- sort(resize(ncTSSrc, 1, "end"))
    names(ncTSS) <- ncTSS$name
    names(ncTSSrc) <- ncTSSrc$name
    names(ncTSS1e) <- ncTSS1e$name
    names(ncTSSrc1e) <- ncTSSrc1e$name

    divHit <- ncTSSrc1e$name[follow(ncTSS1e, ncTSSrc1e)]
    ncTSS$DivTSS <- NA_character_
    ncTSS[ncTSS1e$name]$DivTSS <- divHit
    ncTSS$Dist <- NA_integer_
    haveDiv <- !is.na(ncTSS$DivTSS)
    ncTSS$Dist[haveDiv] <- distance(ncTSS[haveDiv], ncTSSrc[ncTSS$DivTSS[haveDiv]])

    bdSeed <- ncTSS[!is.na(ncTSS$Dist) & ncTSS$Dist <= divergentDist]
  } else {
    bdSeed <- ncTSS[0]
  }

  if (length(bdSeed)) {
    bdSeed <- resize(bdSeed, 1, fix = "end")
    bdSeed$Width <- distance(bdSeed,
      sort(resize(ncTSS, 1, "end"))[bdSeed$DivTSS], ignore.strand = TRUE)
    bdGR <- flank(bdSeed, bdSeed$Width)
    strand(bdGR) <- "*"
    if (length(pcTSS))
      bdGR <- bdGR[!overlapsAny(bdGR, pcTSS, ignore.strand = TRUE)]
    bdGR <- reduce(bdGR)
    bdGR$name <- paste0("BD", as.character(seqnames(bdGR)), "_",
      as.integer(start(bdGR) + ((end(bdGR) - start(bdGR)) / 2)))
  } else {
    bdGR <- GRanges()
    bdGR$name <- character()
  }

  intCandidates <- (anno$TSSTypeLocation == "Intergenic" |
    anno$FeatureType %in% intergenicTypes) & !anno$IsDivergent
  intTSS <- TSS[intCandidates & exprFilt]
  if (length(intTSS) && length(bdGR))
    intTSS <- intTSS[!overlapsAny(resize(intTSS, 1, "start"), bdGR, ignore.strand = TRUE)]
  if (length(intTSS)) {
    intTSS <- reduce(flank(resize(intTSS, 1, fix = "end"), intergenicFlank))
    strand(intTSS) <- "*"
    intTSS$name <- paste0("UD", as.character(seqnames(intTSS)), "_",
      as.integer(start(intTSS) + ((end(intTSS) - start(intTSS)) / 2)))
    if (length(pcTSS))
      intTSS <- intTSS[!overlapsAny(intTSS, pcTSS, ignore.strand = TRUE)]
  } else {
    intTSS <- GRanges()
    intTSS$name <- character()
  }

  enh <- sort(reduce(sort(c(bdGR, intTSS))))
  if (length(enh)) {
    enh$name <- paste0("ENH", as.character(seqnames(enh)), "_",
      as.integer(start(enh) + ((end(enh) - start(enh)) / 2)))
    names(enh) <- enh$name
  }

  list(bdGR = bdGR, intTSS = intTSS, enh = enh)
}

#' Sum TSS quantification over enhancer regions.
#'
#' @param enh A `GRanges` of enhancer regions (e.g. from [defineEnhancers()]).
#' @param TSS A `GRanges` of TSSs aligned to the rows of `quant`.
#' @param quant A numeric matrix or `data.frame` of TSS quantifications
#'   (rows aligned with `TSS`, samples in columns).
#' @param ignoreStrand Whether to ignore strand when overlapping enhancers and
#'   TSSs.
#'
#' @return A numeric matrix of dimensions `length(enh) × ncol(quant)`.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
quantifyEnhancers <- function(enh, TSS, quant, ignoreStrand = TRUE) {
  OVs <- findOverlaps(enh, TSS, ignore.strand = ignoreStrand)
  out <- matrix(0, nrow = length(enh), ncol = ncol(quant),
    dimnames = list(enh$name, colnames(quant)))
  if (!length(OVs)) return(out)
  qH <- queryHits(OVs)
  sH <- subjectHits(OVs)
  uH <- sort(unique(qH))
  for (i in seq_len(ncol(quant))) {
    out[uH, i] <- tapply(sH, qH, function(x) sum(quant[x, i], na.rm = TRUE))
  }
  out
}

#' Compute strand-resolved enhancer statistics.
#'
#' For each enhancer this sums the positive- and negative-strand TSS
#' quantifications separately, and from those works out a Pearson correlation
#' between the two strand profiles across samples, the median strand ratio, the
#' sample where total expression peaks, and (if you pass any) the overlap with
#' peak annotations.
#'
#' @param enh A `GRanges` of enhancer regions.
#' @param TSS A `GRanges` of TSSs aligned to the rows of `quant`.
#' @param quant A numeric matrix or `data.frame` of TSS quantifications.
#' @param exprFilt Logical vector aligned with `quant`/`TSS`.
#' @param peakAnnotations Optional named list of `GRanges`. For each entry
#'   `nm`, an `Overlapping<nm>` column is added.
#'
#' @return A `data.frame` with columns `Enhancer`, `Chr`, `Start`, `End`,
#'   `IsBidirectional`, `StrandCorr`, `MedianStrandRatio`, `ExprPeakSample`,
#'   plus one `Overlapping<name>` column per `peakAnnotations` entry.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
enhancerStrandStats <- function(enh, TSS, quant, exprFilt,
  peakAnnotations = list()) {

  enhP <- enh; strand(enhP) <- "+"
  enhN <- enh; strand(enhN) <- "-"
  TSSf <- TSS[exprFilt]
  quantF <- quant[exprFilt, , drop = FALSE]

  enhQuantP <- matrix(0, nrow = length(enh), ncol = ncol(quant),
    dimnames = list(enh$name, colnames(quant)))
  enhQuantN <- enhQuantP

  ovsP <- findOverlaps(enhP, TSSf, ignore.strand = FALSE)
  ovsN <- findOverlaps(enhN, TSSf, ignore.strand = FALSE)

  for (i in seq_len(ncol(quant))) {
    if (length(ovsP)) {
      qH <- queryHits(ovsP); sH <- subjectHits(ovsP)
      uH <- sort(unique(qH))
      enhQuantP[uH, i] <- tapply(sH, qH, function(x) sum(quantF[x, i], na.rm = TRUE))
    }
    if (length(ovsN)) {
      qH <- queryHits(ovsN); sH <- subjectHits(ovsN)
      uH <- sort(unique(qH))
      enhQuantN[uH, i] <- tapply(sH, qH, function(x) sum(quantF[x, i], na.rm = TRUE))
    }
  }

  strandRatio <- matrix(NA_real_, nrow = length(enh), ncol = ncol(quant),
    dimnames = list(enh$name, colnames(quant)))
  for (i in seq_len(nrow(strandRatio))) {
    iMax <- pmax(enhQuantP[i, ], enhQuantN[i, ])
    iMin <- pmin(enhQuantP[i, ], enhQuantN[i, ])
    strandRatio[i, ] <- (iMax + 1) / (iMin + 1)
  }
  strandRatio[enhQuantP < 2 & enhQuantN < 2] <- NA_real_

  enhCor <- suppressWarnings(stats::cor(t(enhQuantP), t(enhQuantN)))
  pairCor <- diag(enhCor)

  isBidir <- logical(length(enh))
  if (length(ovsP) && length(ovsN)) {
    bothQH <- intersect(queryHits(ovsP), queryHits(ovsN))
    isBidir[bothQH] <- TRUE
  }

  medRatio <- matrixStats::rowMedians(strandRatio, na.rm = TRUE)
  medRatio[!isBidir] <- NA_real_

  totalQuant <- enhQuantP + enhQuantN
  peakSample <- colnames(totalQuant)[apply(totalQuant, 1, which.max)]
  if (!length(peakSample)) peakSample <- character(length(enh))

  out <- data.frame(row.names = NULL,
    Enhancer = enh$name,
    Chr = as.character(seqnames(enh)),
    Start = start(enh),
    End = end(enh),
    IsBidirectional = isBidir,
    StrandCorr = pairCor,
    MedianStrandRatio = medRatio,
    ExprPeakSample = peakSample,
    stringsAsFactors = FALSE
  )

  for (nm in names(peakAnnotations)) {
    col <- paste0("Overlapping", nm)
    out[[col]] <- NA_character_
    pk <- peakAnnotations[[nm]]
    if (is.null(mcols(pk)$name))
      mcols(pk)$name <- paste0(nm, "_", seq_along(pk))
    OVs <- findOverlaps(enh, pk)
    OVs <- OVs[!duplicated(queryHits(OVs))]
    out[[col]][queryHits(OVs)] <- mcols(pk)$name[subjectHits(OVs)]
  }

  out
}

#' Narrow each bidirectional enhancer to its inner-peak coordinates.
#'
#' For each region in `bdGR`, this finds the highest-expression `+` strand and
#' `-` strand TSS overlapping it (using `tssMaxPos` to pick the peak position
#' within each), and then returns a `GRanges` running from the negative-strand
#' peak at its start to the positive-strand peak at its end.
#'
#' @param bdGR A `GRanges` of bidirectional enhancer regions.
#' @param TSS A `GRanges` of TSSs aligned to `quant` and `tssMaxPos`.
#' @param exprFilt Logical vector aligned with `TSS`.
#' @param tssMaxPos Numeric matrix of per-sample peak positions for each TSS
#'   (rows aligned with `TSS`, samples in columns).
#' @param quant Numeric matrix or `data.frame` of TSS quantifications.
#'
#' @return A sorted `GRanges` covering each bidirectional region's inner
#'   peak-to-peak span; metadata column `name` carries the original
#'   `bdGR$name`.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
bidirectionalNarrow <- function(bdGR, TSS, exprFilt, tssMaxPos, quant) {
  if (!length(bdGR)) return(bdGR)

  medMaxPos <- matrixStats::rowMedians(as.matrix(tssMaxPos), na.rm = TRUE)
  tss1 <- TSS
  tss1$thickPos <- round(medMaxPos)
  ok <- !is.na(tss1$thickPos)
  tss1 <- tss1[ok & exprFilt[ok]]
  if (!length(tss1)) return(bdGR[0])
  start(tss1) <- end(tss1) <- tss1$thickPos

  quantOK <- quant[ok & exprFilt[ok], , drop = FALSE]
  tss1$MedianCPM <- matrixStats::rowMedians(as.matrix(quantOK), na.rm = TRUE)

  bdOvs <- findOverlaps(bdGR, tss1, ignore.strand = TRUE)
  if (!length(bdOvs)) return(bdGR[0])

  tss1$BD <- NA_character_
  tss1$BD[subjectHits(bdOvs)] <- bdGR$name[queryHits(bdOvs)]
  tss1 <- tss1[order(tss1$MedianCPM, decreasing = TRUE)]
  tss1p <- tss1[as.character(strand(tss1)) == "+" & !is.na(tss1$BD)]
  tss1m <- tss1[as.character(strand(tss1)) == "-" & !is.na(tss1$BD)]
  tss1p <- tss1p[!duplicated(tss1p$BD)]
  tss1m <- tss1m[!duplicated(tss1m$BD)]
  shared <- intersect(tss1p$BD, tss1m$BD)
  tss1p <- tss1p[tss1p$BD %in% shared]
  tss1m <- tss1m[tss1m$BD %in% shared]

  out <- bdGR
  out$TSSp <- NA_integer_
  out$TSSm <- NA_integer_
  ovsP <- findOverlaps(bdGR, tss1p, ignore.strand = TRUE)
  ovsM <- findOverlaps(bdGR, tss1m, ignore.strand = TRUE)
  out$TSSp[queryHits(ovsP)] <- start(tss1p)[subjectHits(ovsP)]
  out$TSSm[queryHits(ovsM)] <- start(tss1m)[subjectHits(ovsM)]

  out <- out[!is.na(out$TSSp) & !is.na(out$TSSm) & out$TSSm < out$TSSp]
  if (!length(out)) return(out)
  start(out) <- out$TSSm
  end(out) <- out$TSSp
  sort(out)
}
