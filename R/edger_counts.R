#===============================================================================
# edger_counts.R
#===============================================================================

# Imports ======================================================================

#' @import BiocParallel
#' @import edgeR
#' @import GenomicAlignments
#' @import Rsamtools




# Functions ====================================================================

#' @title extract group vector
#'
#' @description extract a group vector (for DGEList) from the read file paths
#'   data
#'
#' @param reads_file_paths list containing the reads file paths for each group
#' @return a group vector for DGEList
#' @export
extract_group_vector <- function(reads_file_paths) {
  unlist(
    lapply(
      names(reads_file_paths),
      function(x) rep(x, length(reads_file_paths[[x]]))
    )
  )
}

#' @title edgeR counts
#'
#' @description read counts for input into edgeR
#'
#' @param peaks Granges object representing consensus peaks
#' @param reads_file_paths list or named character vector giving paths to BAM
#'   files
#' @return DGEList object for use with edgeR
#' @export
edger_counts <- function(peaks, reads_file_paths, group, cores = NULL) {
  if (is.null(cores)) {
    BPPARAM <- MulticoreParam(
      workers = min(length(reads_file_paths), multicoreWorkers())
    )
  } else {
    BPPARAM <- MulticoreParam(
      workers = min(length(reads_file_paths), multicoreWorkers(), cores)
    )
  }
  summarized_experiment <- summarizeOverlaps(
    peaks,
    BamFileList(reads_file_paths),
    mode = "IntersectionNotEmpty",
    ignore.strand = TRUE,
    BPPARAM = BPPARAM
  )
  # DGEList(assays(summarized_experiment)[["counts"]], group = group)
  summarized_experiment
}
