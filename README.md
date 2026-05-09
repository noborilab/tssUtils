# tssUtils

A collection of utilities for analyzing csRNA-seq data, or any type of
transcription start site (TSS) sequencing data. The package wraps the common
operations needed to go from a set of called TSSs and per-sample quantifications
to a fully annotated table, defined enhancer regions, TSS-enhancer correlation
links, and per-sample TSS shape statistics.

## Installation

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("bjmt/tssUtils")
```

`tssUtils` depends on the Bioconductor stack
(`GenomicRanges`, `IRanges`, `GenomeInfoDb`, `S4Vectors`, `rtracklayer`,
`genomation`, `AnnotationDbi`, `GenomicFeatures`). If any are missing,
install them first via `BiocManager::install()`.

## Quick start

```r
library(tssUtils)

# Reference data
TxDb   <- AnnotationDbi::loadDb("Araport11.TxDb")
TE     <- rtracklayer::import("Araport11_TEs.bed")
Tx5utr <- rtracklayer::import("Araport11_5UTR.bed")

txData <- prepTxData(TxDb, seqnames = as.character(1:5))
teData <- prepFeatureData(TE, defaultType = "transposable_element")

# Inputs
TSS   <- rtracklayer::import("tss.final.bed")
quant <- as.matrix(read.delim("tss.final.cpm.txt", row.names = 1))
quant <- quant[TSS$name, ]

# Annotate
anno <- annotateTSS(TSS, txData,
  featureSets   = list(TE = teData),
  tx5utr        = Tx5utr,
  featureCoding = c(TE = "teTSS"))

# Filter, define enhancers, quantify, correlate
exprFilt <- filterByExpression(quant, minCpm = 0.5, minSamples = 3)
enhRes   <- defineEnhancers(TSS, anno, exprFilt)
enhQuant <- quantifyEnhancers(enhRes$enh, TSS, quant)
enhStats <- enhancerStrandStats(enhRes$enh, TSS, quant, exprFilt)
corrDf   <- correlateTSSEnhancers(quant[exprFilt, ], enhQuant,
  TSS[exprFilt], enhRes$enh, minPCC = 0.5)
```

For the full pipeline, including per-sample TSS shape (Shannon entropy,
Simpson diversity, percentile-based thick BED12 coordinates) and writing the
output files, see the vignette:

```r
vignette("annotation", package = "tssUtils")
```

## Function overview

| Stage | Functions |
| --- | --- |
| Reference prep | `prepTxData`, `prepFeatureData` |
| Annotation | `annotateTSS`, `filterByExpression`, `findDivergent`, `findBidirectional` |
| Enhancers | `defineEnhancers`, `quantifyEnhancers`, `enhancerStrandStats`, `bidirectionalNarrow` |
| Correlation | `correlateTSSEnhancers`, `writeBEDPE` |
| Per-sample TSS shape | `tssShape`, `tssShapeFromBigWig`, `aggregateTSSShape`, `tssThickBED`, `shannonEntropy`, `simpsonDiversity` |
| Signal I/O | `readSignal`, `readWindowsUnstranded`, `readWindowsStranded` |
| Misc | `swapStrand`, `calcPctiles` |

## License

GPL (>= 3)
