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

# Integer-count analogue of mkSignalList(), one stranded count GRanges per
# sample, with a triangular sense peak (and optional antisense reads).
mkCountSignalList <- function(depths = c(s1 = 20, s2 = 5),
  samples = names(depths), positions = 100:200, peak = 150,
  antisense = FALSE, seqlen = 5000) {
  si <- GenomeInfoDb::Seqinfo("1", seqlen)
  out <- lapply(samples, function(nm) {
    mult <- depths[[nm]]
    sc <- as.integer(round(mult * pmax(50 - abs(peak - positions), 0)) + 1L)
    gr <- GRanges("1", IRanges(positions, positions), strand = "+")
    gr$score <- sc
    GenomeInfoDb::seqinfo(gr) <- si
    if (isTRUE(antisense)) {
      ap <- positions[positions >= peak - 25 & positions <= peak + 25]
      asc <- as.integer(round(mult * pmax(20 - abs(peak - ap), 0)) + 1L)
      gra <- GRanges("1", IRanges(ap, ap), strand = "-")
      gra$score <- asc
      GenomeInfoDb::seqinfo(gra) <- si
      gr <- sort(c(gr, gra))
    }
    gr
  })
  names(out) <- samples
  out
}

# Write a tiny indexed FASTA to a tempfile and return the path. Plants a
# TATAAA at `plantAt` (1-based) on the forward strand.
mkFastaFixture <- function(len = 400, plantAt = c(50, 250), motif = "TATAAA",
  seed = 1) {
  set.seed(seed)
  g <- paste(sample(c("A", "C", "G", "T"), len, replace = TRUE), collapse = "")
  for (p in plantAt) substr(g, p, p + nchar(motif) - 1L) <- motif
  dna <- Biostrings::DNAStringSet(c("1" = g))
  fa <- tempfile(fileext = ".fa")
  Biostrings::writeXStringSet(dna, fa)
  Rsamtools::indexFa(fa)
  fa
}
