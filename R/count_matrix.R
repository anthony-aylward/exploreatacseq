#===============================================================================
# count_matrix.R
#===============================================================================

# Imports ======================================================================

#' @import BiocParallel
#' @import edgeR
#' @import GenomicAlignments
#' @import limma
#' @import Rsamtools
#' @import SummarizedExperiment




# Functions ====================================================================

#' @title DGEList of read counts
#'
#' @description generate a DGEList from input peaks and BAM files
#'
#' @param peaks Granges object representing consensus peaks
#' @param reads_file_paths list or named character vector giving paths to BAM
#'   files
#' @return DGEList object representing read counts
count_dgelist <- function(peaks, reads_file_paths, group, cores = 1) {
  summarized_experiment <- summarizeOverlaps(
    peaks,
    BamFileList(reads_file_paths),
    mode = "IntersectionNotEmpty",
    ignore.strand = TRUE,
    BPPARAM = MulticoreParam(
      workers = min(length(reads_file_paths), multicoreWorkers(), cores)
    )
  )
  calcNormFactors(
    DGEList(
      assays(summarized_experiment)[["counts"]],
      group = rep(1, times = length(group))
    )
  )
}

#' @title apply transformations to read counts
#'
#' @description apply voom to input read counts (DGEList) and remove batch
#'   effects
#'
#' @param count_dgelist DGEList representing read counts
#' @param batch factor or vector indicating batches
#' @param covariates matrix or vector of numeric covariates
#' @return numeric matrix representing transformed counts
#' @export
transform_counts <- function(
  count_dgelist,
  batch = NULL,
  covariates = NULL
) {
  count_elist <- voom(count_dgelist)
  removeBatchEffect(
    count_elist[["E"]],
    batch = batch,
    covariates = covariates,
    design = count_elist[["design"]]
  )
}

#' @title most variable peaks
#'
#' @description extract the most variable peaks from the count matrix
#'
#' @param count_matrix numeric matrix representing transformed counts
#' @param n max number of peaks to keep
#' @return numeric matrix representing counts for most variable peaks
most_variable_peaks <- function(count_matrix, n = 1e5) {
  count_matrix[
    rev(order(apply(count_matrix, 1, var)))[1:min(n, nrow(count_matrix))],
  ]
}

#' @title extract sample vector
#'
#' @description extract sample vector from count matrix
#'
#' @param count_matrix matrix of read counts
#' @return character vector indicating sample
extract_sample_vector <- function(count_matrix) {
  vapply(
    strsplit(colnames(count_matrix), split = ".", fixed = TRUE),
    function(x) x[[1]],
    character(length = 1)
  )
}

#' @title extract treatment vector
#'
#' @description extract treatment vector from count matrix
#'
#' @param count_matrix matrix of read counts
#' @return character vector indicating treatment
extract_treatment_vector <- function(count_matrix) {
  vapply(
    strsplit(colnames(count_matrix), split = ".", fixed = TRUE),
    function(x) x[[2]],
    character(length = 1)
  )
}
