#===============================================================================
# count_matrix.R
#===============================================================================

# Imports ======================================================================

#' @import BiocParallel
#' @import edgeR
#' @import GenomicAlignments
#' @import SummarizedExperiment




# Functions ====================================================================

#' @title DGEList of read counts
#'
#' @description read counts for input into edgeR
#'
#' @param peaks Granges object representing consensus peaks
#' @param reads_file_paths list or named character vector giving paths to BAM
#'   files
#' @return DGEList object representing read counts
#' @export
count_dgelist <- function(peaks, reads_file_paths, group, cores = NULL) {
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
  calcNormFactors(
    DGEList(
      assays(summarized_experiment)[["counts"]],
      group = group
    )
  )
}

#' @title apply transformations to read counts
#'
#' @description apply voom and remove batch effects
#'
#' @param count_dgelist DGEList representing read counts
#' @return numeric matrix representing transformed counts
#' @export
transform_counts <- function(count_dgelist) {
  count_elist <- voom(counts_dgelist)
  removeBatchEffect(
    count_elist[["E"]],
    batch = x[["batch"]],
    covariates = x[["tss_enrichment"]],
    design = count_elist[["design"]]
  )
}

