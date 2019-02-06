#===============================================================================
# consensus_peaks.R
#===============================================================================

# Imports ======================================================================

#' @import BiocParallel
#' @import edgeR
#' @import GenomicAlignments
#' @import Rsamtools




# Functions ====================================================================

#' @title edgeR counts
#'
#' @description read counts for input into edgeR
#'
#' @param peaks Granges object representing consensus peaks
#' @param reads_file_paths list or named character vector giving paths to BAM
#'   files
#' @return DGEList object for use with edgeR
#' @export
edger_counts <- function(peaks, reads_file_paths) {
  summarized_experiments <- summarizeOverlaps(
    peaks,
    BamFileList(as.character(reads_file_paths)),
    mode = "IntersectionNotEmpty",
    ignore.strand = TRUE,
    BPPARAM = MulticoreParam(
      workers = min(length(reads_file_paths), multicoreWorkers())
    )
  )
  DGEList(
    assays(summarized_experiments)[["counts"]],
    group = names(reads_file_paths)
  )
}