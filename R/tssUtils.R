#' tssUtils: Utilities for working with TSS sequencing data
#'
#' @description
#' A collection of utilities for working with csRNA-seq, or with any TSS
#' sequencing experimental data.
#'
#' @name tssUtils-pkg
#' @aliases tssUtils-pkg
#'
#' @importFrom GenomeInfoDb seqlevels `seqlevels<-`
#' @importFrom genomation ScoreMatrix
#' @importFrom matrixStats rowMins rowMaxs rowMedians
#' @importFrom methods is
#' @importFrom rtracklayer import
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom utils write.table
#' @import GenomicRanges
#' @import IRanges
"_PACKAGE"
