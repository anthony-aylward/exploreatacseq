#!/usr/bin/Rscript
#===============================================================================
# consensus_peaks.R
#===============================================================================

# Imports ======================================================================

#' @import GenomicRanges
#' @import jsonlite
#' @import IRanges




# Functions ====================================================================

#' @title Read peaks from a file
#'
#' @description Read peaks from a BED or similar file into a GRanges object
#'
#' @param peaks_file_path character, path to the peaks file
#' @return GRanges object representing the peaks
#' @export
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
#' @export
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
#' @export
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
#' @export
peaks_by_sample <- function(peaks_paths_by_sample) {
  lapply(peaks_paths_by_sample, unify_peaks)
}

#' @title Consensus peaks
#'
#' @descriptions generate a set of consensus peaks from per-sample peaks
#'
#' @details The consensus set is the union of all peaks identified in any
#'   sample that overlap a peak in another sample.
#'
#' @param peaks list of GRanges objects representing the peaks for each sample
#' @return GRanges, the consensus peak set
#' @export
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

#' @title Consensus peaks from a JSON file
#'
#' @description As `consensus_peaks` but takes as input the path to a JSON file
#'   containing paths to peaks files instead of a list containing peaks
#'
#' @param json_file_path character, the path to a JSON file
#' @return GRanges, the consensus peak set
#' @export
consensus_peaks_from_json <- function(json_file_path) {
  consensus_peaks(peaks_by_sample(fromJSON(json_file_path)))
}
