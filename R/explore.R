#===============================================================================
# explore.R
#===============================================================================

# Imports ======================================================================

#' @import edgeR




# Functions ====================================================================

#' @title Explore ATAC-seq data
#'
#' @description perform an exploratory analysis
#'
#' @param json_file_path path to a JSON file providing data details
#' @param output_prefix character, a prefix for output files
#' @export
explore <- function(json_file_path, output_prefix, treatment_groups = list()) {
  x <- parse_json(json_file_path)
  peaks <- filter_peaks(
    consensus_peaks(peaks_by_sample(x[["peaks_paths_by_sample"]]))
  )
  cat(
    "median peak length: ",
    median_peak_length(peaks),
    "\n",
    file = paste(output_prefix, ".txt", sep = "")
  )
  counts <- transform_counts(
    count_dgelist(
      peaks,
      x[["reads_file_paths"]],
      x[["group"]],
      cores = as.integer(opt[["cores"]])
    ),
    batch = x[["batch"]],
    covariates = x[["tss_enrichment"]]
  )
  save(counts, file = paste(output_prefix, ".RData", sep = ""))
  treatment <- sapply(
    strsplit(colnames(counts), split = ".", fixed = TRUE),
    function(x) x[[2]]
  )
  for (group in treatment_groups) {
    pca <- two_principal_components(counts[,treatment %in% group])
    pdf(
      paste(
        output_prefix,
        "-",
        paste(group, sep = "-", collapse = "-"),
        ".pdf",
        sep = ""
      ),
      height = 14
    )
    plot_pca(pca, draw_lines = list(group))
    dev.off()
    png(
      paste(
        output_prefix,
        "-",
        paste(group, sep = "-", collapse = "-"),
        ".png",
        sep = ""
      ),
      height = 960
    )
    plot_pca(pca, draw_lines = list(group))
    dev.off()
  }
}