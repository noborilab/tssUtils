suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
})

mkTSS <- function(seqs = "1", starts, strands = "+", names = NULL) {
  gr <- GRanges(seqs, IRanges(starts, starts), strand = strands)
  if (is.null(names)) names <- paste0("TSS_", seq_along(gr))
  gr$name <- names
  names(gr) <- gr$name
  gr
}

mkTxData <- function() {
  Tx <- GRanges("1",
    IRanges(start = c(1000, 2500, 5000, 8000, 12000),
      end =      c(2000, 3500, 6000, 9000, 13000)),
    strand = c("+", "+", "-", "+", "-"),
    tx_name = c("tx1.1", "tx2.1", "tx3.1", "tx4.1", "tx5.1"))
  TxAnn <- data.frame(row.names = Tx$tx_name,
    TXNAME = Tx$tx_name,
    TXTYPE = c("mRNA", "mRNA", "mRNA", "lnc_RNA", "mRNA"),
    stringsAsFactors = FALSE)
  TxTSS <- promoters(Tx, upstream = 0, downstream = 1)
  TxRC <- tssUtils::swapStrand(Tx)
  TxRCTSS <- promoters(TxRC, upstream = 0, downstream = 1)
  names(Tx) <- Tx$tx_name
  names(TxTSS) <- TxTSS$tx_name
  names(TxRC) <- TxRC$tx_name
  names(TxRCTSS) <- TxRCTSS$tx_name
  list(feat = Tx, featTSS = TxTSS, featRC = TxRC, featRCTSS = TxRCTSS,
    featAnn = TxAnn, nameCol = "tx_name", typeCol = "TXTYPE",
    defaultType = NA_character_)
}

mkFeatureData <- function() {
  Feat <- GRanges("1",
    IRanges(start = c(15000, 20000),
      end =      c(15500, 20500)),
    strand = c("+", "-"),
    name = c("TE1", "TE2"))
  FeatTSS <- promoters(Feat, upstream = 0, downstream = 1)
  FeatRC <- tssUtils::swapStrand(Feat)
  FeatRCTSS <- promoters(FeatRC, upstream = 0, downstream = 1)
  names(Feat) <- Feat$name
  names(FeatTSS) <- FeatTSS$name
  names(FeatRC) <- FeatRC$name
  names(FeatRCTSS) <- FeatRCTSS$name
  list(feat = Feat, featTSS = FeatTSS, featRC = FeatRC, featRCTSS = FeatRCTSS,
    featAnn = NULL, nameCol = "name", typeCol = NA_character_,
    defaultType = "transposable_element")
}

mkQuant <- function(rowNames, nSamples = 3, fill = 1) {
  m <- matrix(fill, nrow = length(rowNames), ncol = nSamples,
    dimnames = list(rowNames, paste0("s", seq_len(nSamples))))
  m
}

mkSignalList <- function(samples = c("s1", "s2"), positions = 100:200,
  scoreFn = function(p) 100 - abs(150 - p)) {
  out <- lapply(samples, function(nm) {
    gr <- GRanges("1", IRanges(positions, positions), strand = "+")
    gr$score <- scoreFn(positions)
    gr
  })
  names(out) <- samples
  out
}
