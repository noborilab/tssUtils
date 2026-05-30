#' Prepare reference data from a TxDb for TSS annotation.
#'
#' Pulls the transcripts and their TSS positions out of a TxDb (along with the
#' reverse-complement variants, which the antisense rules need), and builds a
#' transcript-to-type lookup to go with them. The returned list is what
#' [annotateTSS()] expects to be handed.
#'
#' @param TxDb A `TxDb` object (e.g. loaded with [AnnotationDbi::loadDb()]).
#' @param seqnames Optional character vector of seqlevels to restrict the
#'   transcripts to. If `NULL`, all seqlevels are kept.
#' @param geneNameCol Name of the metadata column on `transcripts(TxDb)` that
#'   holds the canonical transcript identifier.
#' @param annoCol Column to extract from `TxDb` for transcript-type annotation.
#' @param annoKey Key column used to look up `annoCol` values.
#'
#' @return A named list with elements `feat`, `featTSS`, `featRC`, `featRCTSS`,
#'   `featAnn` (a `data.frame` indexed by transcript ID), `nameCol`, `typeCol`,
#'   and `defaultType` (`NA_character_`).
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
prepTxData <- function(TxDb, seqnames = NULL, geneNameCol = "tx_name",
  annoCol = "TXTYPE", annoKey = "TXNAME") {

  Tx <- GenomicFeatures::transcripts(TxDb)
  if (!is.null(seqnames)) {
    Tx <- Tx[as.character(seqnames(Tx)) %in% as.character(seqnames)]
    seqlevels(Tx) <- as.character(seqnames)
  }
  TxTSS <- promoters(Tx, upstream = 0, downstream = 1)
  TxRC <- swapStrand(Tx)
  TxRCTSS <- promoters(TxRC, upstream = 0, downstream = 1)

  TxAnn <- suppressMessages(AnnotationDbi::select(TxDb,
    keys = mcols(Tx)[[geneNameCol]], columns = annoCol, keytype = annoKey))
  rownames(TxAnn) <- TxAnn[[annoKey]]

  names(Tx) <- mcols(Tx)[[geneNameCol]]
  names(TxTSS) <- mcols(TxTSS)[[geneNameCol]]
  names(TxRC) <- mcols(TxRC)[[geneNameCol]]
  names(TxRCTSS) <- mcols(TxRCTSS)[[geneNameCol]]

  list(
    feat = Tx,
    featTSS = TxTSS,
    featRC = TxRC,
    featRCTSS = TxRCTSS,
    featAnn = TxAnn,
    nameCol = geneNameCol,
    typeCol = annoCol,
    defaultType = NA_character_
  )
}

#' Prepare reference data from a `GRanges` of features for TSS annotation.
#'
#' The generic counterpart to [prepTxData()], meant for features that do not
#' come from a TxDb (transposable elements loaded from a BED file, say). It
#' produces the same list shape, so that [annotateTSS()] can treat every
#' reference set in the same way without having to know where it came from.
#'
#' @param features A `GRanges` of features.
#' @param nameCol Name of the metadata column on `features` that holds the
#'   canonical feature identifier.
#' @param defaultType Single string used as the `FeatureType` for every
#'   feature (e.g. `"transposable_element"`).
#'
#' @return A list with the same shape as [prepTxData()], with `featAnn = NULL`.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
prepFeatureData <- function(features, nameCol = "name",
  defaultType = "feature") {

  if (is.null(mcols(features)[[nameCol]]))
    stop("features must have an mcol named '", nameCol, "'")

  Feat <- features
  FeatTSS <- promoters(Feat, upstream = 0, downstream = 1)
  FeatRC <- swapStrand(Feat)
  FeatRCTSS <- promoters(FeatRC, upstream = 0, downstream = 1)

  names(Feat) <- mcols(Feat)[[nameCol]]
  names(FeatTSS) <- mcols(FeatTSS)[[nameCol]]
  names(FeatRC) <- mcols(FeatRC)[[nameCol]]
  names(FeatRCTSS) <- mcols(FeatRCTSS)[[nameCol]]

  list(
    feat = Feat,
    featTSS = FeatTSS,
    featRC = FeatRC,
    featRCTSS = FeatRCTSS,
    featAnn = NULL,
    nameCol = nameCol,
    typeCol = NA_character_,
    defaultType = defaultType
  )
}

#' Filter quantification by minimum CPM across samples.
#'
#' @param quant A numeric matrix or `data.frame` with TSSs in rows and
#'   samples in columns.
#' @param minCpm Minimum CPM value to count a sample as expressing.
#' @param minSamples Minimum number of samples that must satisfy the CPM
#'   threshold.
#'
#' @return A logical vector of length `nrow(quant)` indicating which TSSs pass.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
filterByExpression <- function(quant, minCpm = 0.5, minSamples = 3) {
  apply(quant, 1, function(x) sum(x >= minCpm, na.rm = TRUE) >= minSamples)
}

# Internal: vectorised overlap between x[i] and y[y_ind[i]] when y_ind contains
# names of y. Returns a logical of length(x).
.checkOverlapVec <- function(x, y, y_ind) {
  i <- !is.na(y_ind) & y_ind %in% names(y)
  o <- logical(length(x))
  o[i] <- as.logical(poverlaps(x[i], y[y_ind[i]]))
  o
}

# Internal: returns TRUE for entries of x that fall strictly upstream of
# y[y_ind[i]] in a strand-aware sense.
.checkIfUpstreamVec <- function(x, y, y_ind) {
  i <- !is.na(y_ind) & y_ind %in% names(y)
  o <- logical(length(x))
  o[i] <- ifelse(as.character(strand(x[i])) == "+",
    end(x[i]) < start(y[y_ind[i]]),
    start(x[i]) > end(y[y_ind[i]]))
  o
}

# Internal: build the empty Anno data.frame.
.initAnno <- function(TSS, peakAnnotations) {
  Anno <- data.frame(row.names = TSS$name,
    TSS = TSS$name,
    Chr = as.character(seqnames(TSS)),
    Start = start(TSS),
    End = end(TSS),
    Strand = as.character(strand(TSS)),
    TSSTypeCoding = "ncTSS",
    TSSTypeLocation = NA_character_,
    FeatureID = NA_character_,
    FeatureSymbol = NA_character_,
    FeatureType = NA_character_,
    TSSFeatOrientation = NA_character_,
    DistanceToFeature = NA_integer_,
    stringsAsFactors = FALSE
  )
  for (nm in names(peakAnnotations)) {
    col <- paste0("Overlapping", nm)
    Anno[[col]] <- NA_character_
    pk <- peakAnnotations[[nm]]
    if (is.null(mcols(pk)$name))
      mcols(pk)$name <- paste0(nm, "_", seq_along(pk))
    OVs <- findOverlaps(TSS, pk)
    OVs <- OVs[!duplicated(queryHits(OVs))]
    Anno[[col]][queryHits(OVs)] <- mcols(pk)$name[subjectHits(OVs)]
  }
  Anno
}

# Rule: sense Tx promoter (5'UTR-aware).
.classifyTxPromoter <- function(TSS, Anno, txData, params) {
  TxTSS <- txData$featTSS
  TxAnn <- txData$featAnn
  Tx5utr <- params$tx5utr

  if (!is.null(Tx5utr)) {
    Overlapping5utr <- findOverlaps(TSS, Tx5utr)
    overlap5 <- rep(NA_character_, length(TSS))
    overlap5[queryHits(Overlapping5utr)] <- mcols(Tx5utr)$name[subjectHits(Overlapping5utr)]
  } else {
    overlap5 <- rep(NA_character_, length(TSS))
  }

  Nearest <- distanceToNearest(TSS, TxTSS, ignore.strand = FALSE)
  nearestId <- rep(NA_character_, length(TSS))
  nearestDist <- rep(NA_integer_, length(TSS))
  nearestId[queryHits(Nearest)] <- names(TxTSS)[subjectHits(Nearest)]
  nearestDist[queryHits(Nearest)] <- mcols(Nearest)$distance

  overlapsUtrOfNearest <- if (!is.null(Tx5utr))
    .checkOverlapVec(TSS, Tx5utr, nearestId) else logical(length(TSS))
  isUpstreamNearest <- .checkIfUpstreamVec(TSS, TxTSS, nearestId)

  ok <- overlapsUtrOfNearest |
    (isUpstreamNearest & !is.na(nearestDist) & nearestDist <= params$promoterUpstream) |
    (!is.na(nearestDist) & nearestDist == 0)

  okFix <- !ok & !is.na(overlap5)
  if (any(okFix)) {
    nearestId[okFix] <- overlap5[okFix]
    nearestDist[okFix] <- distance(TSS[okFix], TxTSS[overlap5[okFix]])
    isUpstreamNearest[okFix] <- FALSE
    ok[okFix] <- TRUE
  }

  ok[is.na(ok)] <- FALSE

  hitNames <- TSS$name[ok]
  if (length(hitNames)) {
    nearestType <- TxAnn[nearestId[ok], txData$typeCol]
    isPC <- !is.na(nearestType) & nearestType %in% params$mRNATypes
    Anno[hitNames, "TSSTypeCoding"][isPC] <- "pcTSS"
    Anno[hitNames, "TSSTypeLocation"] <- "Promoter"
    Anno[hitNames, "FeatureID"] <- nearestId[ok]
    Anno[hitNames, "TSSFeatOrientation"] <- "Sense"
    Anno[hitNames, "FeatureType"] <- nearestType
    d <- nearestDist[ok]
    d <- ifelse(isUpstreamNearest[ok], -d, d)
    Anno[hitNames, "DistanceToFeature"] <- d
  }

  list(remaining = TSS[!ok], Anno = Anno)
}

# Rule: sense feature promoter (e.g. TE promoter). Generic over any
# featureSet shape.
.classifyFeatPromoter <- function(TSS, Anno, featData, params) {
  FeatTSS <- featData$featTSS
  Feat <- featData$feat

  Nearest <- distanceToNearest(TSS, FeatTSS, ignore.strand = FALSE)
  nearestId <- rep(NA_character_, length(TSS))
  nearestDist <- rep(NA_integer_, length(TSS))
  nearestId[queryHits(Nearest)] <- names(FeatTSS)[subjectHits(Nearest)]
  nearestDist[queryHits(Nearest)] <- mcols(Nearest)$distance

  overlapsNearest <- .checkOverlapVec(TSS, FeatTSS, nearestId)
  isUpstreamNearest <- .checkIfUpstreamVec(TSS, Feat, nearestId)

  ok <- overlapsNearest |
    (isUpstreamNearest & !is.na(nearestDist) & nearestDist <= params$promoterUpstream)
  ok[is.na(ok)] <- FALSE

  hitNames <- TSS$name[ok]
  if (length(hitNames)) {
    Anno[hitNames, "TSSTypeCoding"] <- params$featCoding
    Anno[hitNames, "TSSTypeLocation"] <- "Promoter"
    Anno[hitNames, "FeatureID"] <- nearestId[ok]
    Anno[hitNames, "TSSFeatOrientation"] <- "Sense"
    Anno[hitNames, "FeatureType"] <- featData$defaultType
    d <- nearestDist[ok]
    d <- ifelse(isUpstreamNearest[ok], -d, d)
    Anno[hitNames, "DistanceToFeature"] <- d
  }

  list(remaining = TSS[!ok], Anno = Anno)
}

# Rule: any other Tx promoter (close enough to, or overlapping, the nearest TSS).
.classifyOtherTxPromoter <- function(TSS, Anno, txData, params) {
  TxTSS <- txData$featTSS
  TxAnn <- txData$featAnn

  Nearest <- distanceToNearest(TSS, TxTSS, ignore.strand = FALSE)
  nearestId <- rep(NA_character_, length(TSS))
  nearestDist <- rep(NA_integer_, length(TSS))
  nearestId[queryHits(Nearest)] <- names(TxTSS)[subjectHits(Nearest)]
  nearestDist[queryHits(Nearest)] <- mcols(Nearest)$distance

  overlapsNearest <- .checkOverlapVec(TSS, TxTSS, nearestId)
  isUpstreamNearest <- .checkIfUpstreamVec(TSS, TxTSS, nearestId)

  ok <- overlapsNearest | (!is.na(nearestDist) & nearestDist <= params$promoterUpstream)
  ok[is.na(ok)] <- FALSE

  hitNames <- TSS$name[ok]
  if (length(hitNames)) {
    Anno[hitNames, "TSSTypeLocation"] <- "Promoter"
    Anno[hitNames, "FeatureID"] <- nearestId[ok]
    Anno[hitNames, "TSSFeatOrientation"] <- "Sense"
    Anno[hitNames, "FeatureType"] <- TxAnn[nearestId[ok], txData$typeCol]
    d <- nearestDist[ok]
    d <- ifelse(isUpstreamNearest[ok], -d, d)
    Anno[hitNames, "DistanceToFeature"] <- d
  }

  list(remaining = TSS[!ok], Anno = Anno)
}

# Rule: intragenic, sense (the TSS overlaps a feature body).
.classifyIntragenic <- function(TSS, Anno, featData, params) {
  Feat <- featData$feat
  FeatTSS <- featData$featTSS
  FeatAnn <- featData$featAnn
  typeCol <- featData$typeCol

  OVs <- findOverlaps(TSS, Feat)
  OVs <- OVs[!duplicated(queryHits(OVs))]
  if (!length(OVs)) return(list(remaining = TSS, Anno = Anno))

  hit <- queryHits(OVs)
  sub <- subjectHits(OVs)
  hitNames <- TSS$name[hit]
  featNames <- names(Feat)[sub]
  d <- distance(TSS[hit], FeatTSS[sub])

  Anno[hitNames, "TSSTypeLocation"] <- "Intragenic"
  Anno[hitNames, "FeatureID"] <- featNames
  Anno[hitNames, "TSSFeatOrientation"] <- "Sense"
  if (!is.null(FeatAnn) && !is.na(typeCol)) {
    Anno[hitNames, "FeatureType"] <- FeatAnn[featNames, typeCol]
  } else {
    Anno[hitNames, "FeatureType"] <- featData$defaultType
  }
  Anno[hitNames, "DistanceToFeature"] <- d

  ok <- logical(length(TSS))
  ok[hit] <- TRUE
  list(remaining = TSS[!ok], Anno = Anno)
}

# Rule: antisense ncTSS against a feature set's reverse-complement view.
.classifyAntisense <- function(TSS, Anno, featData, params, codingLabel) {
  FeatRC <- featData$featRC
  FeatRCTSS <- featData$featRCTSS
  FeatAnn <- featData$featAnn
  typeCol <- featData$typeCol
  defaultType <- featData$defaultType

  Nearest <- distanceToNearest(TSS, FeatRCTSS, ignore.strand = FALSE)
  OVs <- findOverlaps(TSS, FeatRC)
  OVs <- OVs[!duplicated(queryHits(OVs))]

  overlapId <- rep(NA_character_, length(TSS))
  overlapId[queryHits(OVs)] <- names(FeatRC)[subjectHits(OVs)]

  nearestId <- overlapId
  fillIdx <- is.na(nearestId)
  if (length(Nearest)) {
    nearestFromQ <- rep(NA_character_, length(TSS))
    nearestFromQ[queryHits(Nearest)] <- names(FeatRCTSS)[subjectHits(Nearest)]
    nearestId[fillIdx] <- nearestFromQ[fillIdx]
  }

  nearestDist <- rep(NA_integer_, length(TSS))
  haveId <- !is.na(nearestId) & nearestId %in% names(FeatRCTSS)
  if (any(haveId)) {
    nearestDist[haveId] <- distance(TSS[haveId], FeatRCTSS[nearestId[haveId]])
  }

  isUpstreamNearest <- .checkIfUpstreamVec(TSS, FeatRCTSS, nearestId)

  ok <- !is.na(overlapId) |
    (isUpstreamNearest & !is.na(nearestDist) & nearestDist <= params$promoterUpstream)
  ok[is.na(ok)] <- FALSE

  hitNames <- TSS$name[ok]
  if (length(hitNames)) {
    Anno[hitNames, "TSSTypeCoding"] <- codingLabel
    Anno[hitNames, "TSSTypeLocation"] <- "Intragenic"
    Anno[hitNames, "TSSFeatOrientation"] <- "Antisense"
    Anno[hitNames, "FeatureID"] <- nearestId[ok]
    if (!is.null(FeatAnn) && !is.na(typeCol)) {
      Anno[hitNames, "FeatureType"] <- FeatAnn[nearestId[ok], typeCol]
    } else {
      Anno[hitNames, "FeatureType"] <- defaultType
    }
    d <- nearestDist[ok]
    d <- ifelse(isUpstreamNearest[ok], -d, d)
    Anno[hitNames, "DistanceToFeature"] <- d
  }

  list(remaining = TSS[!ok], Anno = Anno)
}

# Rule: divergence column population using findDivergent().
.classifyDivergent <- function(Anno, divergentDist) {
  Anno$IsDivergent <- FALSE
  Anno$DivergentDistance <- NA_integer_
  Anno$DivergentTSS <- NA_character_

  isPC <- Anno$TSSTypeCoding == "pcTSS"
  isNC <- !isPC
  if (!any(isPC) || !any(isNC)) return(Anno)

  pcGR <- GRanges(Anno$Chr[isPC], IRanges(Anno$Start[isPC], Anno$End[isPC]),
    Anno$Strand[isPC], name = Anno$TSS[isPC])
  ncGR <- GRanges(Anno$Chr[isNC], IRanges(Anno$Start[isNC], Anno$End[isNC]),
    Anno$Strand[isNC], name = Anno$TSS[isNC])

  div <- findDivergent(pcGR, ncGR, maxDist = divergentDist, returnMerged = FALSE)
  pcDiv <- div$pcTSS
  ncDiv <- div$ncTSS

  if (length(pcDiv)) {
    Anno[pcDiv$name, "IsDivergent"] <- TRUE
    Anno[pcDiv$name, "DivergentDistance"] <- pcDiv$Dist
    Anno[pcDiv$name, "DivergentTSS"] <- pcDiv$DivTSS
  }
  if (length(ncDiv)) {
    Anno[ncDiv$name, "IsDivergent"] <- TRUE
    Anno[ncDiv$name, "DivergentDistance"] <- ncDiv$Dist
    Anno[ncDiv$name, "DivergentTSS"] <- ncDiv$DivTSS
  }

  Anno
}

#' Annotate TSSs against transcript and feature reference data.
#'
#' Each TSS is classified by looking at how it overlaps, and how close it sits
#' to, the transcript and feature sets you supply. The rules are applied in a
#' fixed order of precedence, from highest to lowest: sense Tx promoter (which
#' is 5'UTR-aware), sense feature promoter, any other sense Tx promoter,
#' intragenic sense Tx, intragenic sense feature, antisense Tx, antisense
#' feature, and finally intergenic (the catch-all for anything that is left
#' over). Once everything has been classified, divergent protein-coding /
#' non-coding TSS pairs are flagged with [findDivergent()].
#'
#' @param TSS A `GRanges` of called TSSs. Must have `mcols(TSS)$name` set; it
#'   will be set to `paste0("TSS_", seq_along(TSS))` if missing.
#' @param txData Output of [prepTxData()].
#' @param featureSets Named list of [prepFeatureData()] outputs (e.g.
#'   `list(TE = teData)`). The list order is the order classification rules
#'   are applied. May be empty.
#' @param tx5utr Optional `GRanges` of 5'UTR ranges (must have `mcols$name`
#'   matching transcript IDs in `txData$feat`).
#' @param peakAnnotations Optional named list of `GRanges` (e.g. enhancer or
#'   ATAC peak annotations). For each entry `nm`, an `Overlapping<nm>` column
#'   is added to the result.
#' @param geneSymbols Optional `data.frame` with columns `Gene` and `Symbol`
#'   giving a transcript-ID-stripped → symbol lookup.
#' @param geneFromTxFn Function applied to `FeatureID` to derive a gene ID
#'   that can be looked up in `geneSymbols$Gene`. Defaults to stripping a
#'   trailing `.<digits>`.
#' @param teGenes Optional character vector of gene IDs (post-`geneFromTxFn`)
#'   that should be reclassified as `transposable_element_gene` /
#'   `teTSS` even when they were classified as `pcTSS` based on `mRNATypes`.
#' @param featureCoding Named character vector mapping each feature-set name
#'   to the `TSSTypeCoding` value to assign on a sense promoter / antisense
#'   match (e.g. `c(TE = "teTSS")`). Defaults to `"teTSS"` for any unnamed
#'   entries.
#' @param mRNATypes Character vector of `featAnn[, txData$typeCol]` values
#'   that mark a sense Tx promoter as `pcTSS`. Other values stay `ncTSS`.
#' @param promoterUpstream,promoterDownstream Distance thresholds (bp) for
#'   the promoter rules. Sense Tx promoter accepts overlap with the
#'   transcript 5'UTR, distance 0, or upstream within `promoterUpstream`.
#'   Other promoter rules use `promoterUpstream` only. `promoterDownstream`
#'   is reserved for future extensions.
#' @param divergentDist Maximum distance between paired antisense pcTSS/ncTSS
#'   to flag as divergent.
#'
#' @return A `data.frame` with one row per `TSS` and columns: `TSS`, `Chr`,
#'   `Start`, `End`, `Strand`, `TSSTypeCoding`, `TSSTypeLocation`,
#'   `FeatureID`, `FeatureSymbol`, `FeatureType`, `TSSFeatOrientation`,
#'   `DistanceToFeature`, `IsDivergent`, `DivergentDistance`, `DivergentTSS`,
#'   plus one `Overlapping<name>` column per `peakAnnotations` entry.
#'
#' @author Benjamin Jean-Marie Tremblay, \email{benjamin.tremblay@tsl.ac.uk}
#' @export
annotateTSS <- function(TSS, txData, featureSets = list(), tx5utr = NULL,
  peakAnnotations = list(), geneSymbols = NULL,
  geneFromTxFn = function(x) sub("\\.\\d+$", "", x),
  teGenes = character(), featureCoding = character(),
  mRNATypes = "mRNA", promoterUpstream = 200, promoterDownstream = 200,
  divergentDist = 500) {

  if (is.null(mcols(TSS)$name)) TSS$name <- paste0("TSS_", seq_along(TSS))
  names(TSS) <- TSS$name

  Anno <- .initAnno(TSS, peakAnnotations)

  params <- list(
    tx5utr = tx5utr,
    mRNATypes = mRNATypes,
    promoterUpstream = promoterUpstream,
    promoterDownstream = promoterDownstream
  )

  step <- .classifyTxPromoter(TSS, Anno, txData, params)
  TSS <- step$remaining; Anno <- step$Anno

  for (nm in names(featureSets)) {
    fd <- featureSets[[nm]]
    coding <- if (nm %in% names(featureCoding)) featureCoding[[nm]] else "teTSS"
    fparams <- c(params, list(featCoding = coding))
    step <- .classifyFeatPromoter(TSS, Anno, fd, fparams)
    TSS <- step$remaining; Anno <- step$Anno
  }

  step <- .classifyOtherTxPromoter(TSS, Anno, txData, params)
  TSS <- step$remaining; Anno <- step$Anno

  step <- .classifyIntragenic(TSS, Anno, txData, params)
  TSS <- step$remaining; Anno <- step$Anno

  for (nm in names(featureSets)) {
    fd <- featureSets[[nm]]
    step <- .classifyIntragenic(TSS, Anno, fd, params)
    TSS <- step$remaining; Anno <- step$Anno
  }

  step <- .classifyAntisense(TSS, Anno, txData, params, codingLabel = "ncTSS")
  TSS <- step$remaining; Anno <- step$Anno

  for (nm in names(featureSets)) {
    fd <- featureSets[[nm]]
    coding <- if (nm %in% names(featureCoding)) featureCoding[[nm]] else "teTSS"
    step <- .classifyAntisense(TSS, Anno, fd, params, codingLabel = coding)
    TSS <- step$remaining; Anno <- step$Anno
  }

  Anno$TSSTypeLocation[is.na(Anno$TSSTypeLocation)] <- "Intergenic"

  if (!is.null(geneSymbols) && all(c("Gene", "Symbol") %in% colnames(geneSymbols))) {
    sym <- structure(geneSymbols$Symbol, names = geneSymbols$Gene)
    Anno$FeatureSymbol <- sym[geneFromTxFn(Anno$FeatureID)]
  }

  if (length(teGenes)) {
    isTEG <- !is.na(Anno$FeatureID) & geneFromTxFn(Anno$FeatureID) %in% teGenes
    Anno$FeatureType[isTEG] <- "transposable_element_gene"
    Anno$TSSTypeCoding[isTEG & Anno$TSSTypeLocation == "Promoter"] <- "teTSS"
  }

  Anno <- .classifyDivergent(Anno, divergentDist)

  Anno
}
