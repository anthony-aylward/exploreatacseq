#===============================================================================
# consensus_peaks.R
#===============================================================================

# Imports ======================================================================

#' @importFrom BiocGenerics %in%
#' @importFrom GenomicRanges GRanges union width seqnames
#' @importFrom IRanges IRanges subsetByOverlaps
#' @importFrom S4Vectors Rle
#' @importFrom stats median
#' @importFrom utils combn read.table


# Functions ====================================================================

#' @title Read peaks from a file
#'
#' @description Read peaks from a BED, narrowPeak, or similar file into a
#'   GRanges object
#'
#' @param peaks_file_path character, path to the peaks file
#' @return GRanges object representing the peaks
read_peaks <- function(peaks_file_path) {
  peaks_df <- read.table(peaks_file_path, stringsAsFactors = FALSE)
  GRanges(
    seqnames = Rle(peaks_df[[1]]),
    ranges = IRanges(start = peaks_df[[2]], end = peaks_df[[3]])
  )
}

#' @title Unify GRanges
#'
#' @description Union of several GRanges objects
#'
#' @param x list of GRanges objects
#' @param ignore_strand logical, if TRUE strand information will be ignored
#' @return GRanges, the union of the GRanges in x
unify_granges <- function(x, ignore_strand = FALSE) {
  if (length(x) > 1) {
    Reduce(function(x, y) union(x, y, ignore.strand = ignore_strand), x)
  } else {
    x[[1]]
  }
}

#' @title Unify peaks
#'
#' @description Create a GRanges object representing the union of peaks from
#'  one or more input files
#'
#' @param peaks_file_paths character, paths to the files containing the peaks
#' @return GRanges, the union of peaks from the input files
unify_peaks <- function(peaks_file_paths) {
  unify_granges(lapply(peaks_file_paths, read_peaks))
}

#' @title Peaks by sample
#'
#' @description Get a GRanges object containing a set of peaks for each input
#'   sample
#'
#' @param peaks_paths_by_sample a list of character vectors containing paths of
#'   peaks files for each sample
#' @return list of GRanges objects representing the peaks for each sample
peaks_by_sample <- function(peaks_paths_by_sample) {
  lapply(peaks_paths_by_sample, unify_peaks)
}

#' @title Consensus peaks
#'
#' @description generate a set of consensus peaks from per-sample peaks
#'
#' @details The consensus set is the union of all peaks identified in any
#'   sample that overlap a peak in another sample.
#'
#' @param peaks list of GRanges objects representing the peaks for each sample
#' @return GRanges, the consensus peak set
consensus_peaks <- function(peaks) {
  unify_granges(
    lapply(
      combn(peaks, 2, simplify = FALSE),
      function(pair) {
        union(
          subsetByOverlaps(pair[[1]], pair[[2]], ignore.strand = TRUE),
          subsetByOverlaps(pair[[2]], pair[[1]], ignore.strand = TRUE),
          ignore.strand = TRUE
        )
      }
    )
  )
}

#' @title Filter peaks
#'
#' @description filter out peaks that are greater than 3kb in length or are
#'   non-autosomal
#'
#' @param peaks GRanges object representing peaks
#' @return GRanges, the filtered peaks
filter_peaks <- function(peaks) {
  peaks[
    (width(peaks) <= 3000)
    & (seqnames(peaks) %in% paste("chr", as.character(seq_len(22)), sep = ""))
  ]
}

#' @title Median Peak Length
#'
#' @description get the median peak length
#'
#' @param peaks Granges object representing peaks
#' @return integer, the median peak length
median_peak_length <- function(peaks) {
  median(width(peaks))
}
