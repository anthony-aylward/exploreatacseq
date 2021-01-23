#===============================================================================
# explore.R
#===============================================================================

#' @import svglite




# Functions ====================================================================

#' @title preprocess ATAC-seq data
#'
#' @description preprocess multiple ATAC-seq datasets into a read count matrix
#' @param input_data list contatining data details or character path to a JSON
#'   file
#' @param n_peaks integer, number of the top most variable peaks to include in
#'   the analysis
#' @param cores integer, max number of cores to use
#' @return list giving the median peak length and the matrix of read counts
#' @export
preprocess <- function(input_data, n_peaks = 1e5, cores = 1) {
  if (typeof(input_data) == "list") {
    x <- input_data
  } else if (typeof(input_data) == "character") {
    x <- parse_json(input_data)
  } else {
    stop("input data must be either list or character")
  }
  peaks <- filter_peaks(
    consensus_peaks(peaks_by_sample(x[["peaks_paths_by_sample"]]))
  )
  counts <- most_variable_peaks(
    transform_counts(
      count_dgelist(
        peaks,
        x[["reads_file_paths"]],
        x[["group"]],
        cores = cores
      ),
      batch = x[["batch"]],
      covariates = x[["tss_enrichment"]]
    ),
    n = n_peaks
  )
  list(median_peak_length = median_peak_length(peaks), counts = counts)
}

#' @title generate PCA plots
#'
#' @description plot PCA of the read counts
#'
#' @param counts matrix of read counts
#' @param output_prefix prefix for output files
#' @param treatment character vector indicating treatment, if NULL, it will be
#'   inferred from the read count matrix
#' @param treatment_groups list of treatment groups
#' @param labels if TRUE, text labels will be added to points in all plots
#' @param palette_order ordering of color palette, eithier "categorical" or
#'   "sequential"
#' @param palette the color palette
#' @export
generate_pca_plots <- function(
  counts,
  output_prefix,
  treatment = NULL,
  treatment_order = NULL,
  treatment_groups = list(),
  labels = FALSE,
  palette_order = "categorical",
  palette = NULL
) {
  if (is.null(treatment)) treatment <- extract_treatment_vector(counts)
  if (is.null(treatment_order)) treatment_order <- unique(treatment)
  if (is.null(palette)) {
    palette_vector = setNames(
      exploreatacseq_color_palette(order = palette_order)[
        1:length(treatment_order)
      ],
      treatment_order
    )
  } else {
    palette_vector = setNames(
      palette[1:length(treatment_order)],
      treatment_order
    )
  }
  pca <- prcomp(counts, rank = 2)
  svglite(paste(output_prefix, "-pca.svg", sep = ""), height = 7, width = 7)
  plot_pca(pca, labels = labels, palette = palette_vector)
  dev.off()
  pdf(paste(output_prefix, "-pca.pdf", sep = ""))
  plot_pca(pca, labels = labels, palette = palette_vector)
  dev.off()
  png(paste(output_prefix, "-pca.png", sep = ""))
  plot_pca(pca, labels = labels, palette = palette_vector)
  dev.off()
  for (group in treatment_groups) {
    palette = palette_vector[group]
    pca <- prcomp(counts[,treatment %in% group], rank = 2)
    svglite(
      paste(
        output_prefix,
        "-",
        paste(group, sep = "-", collapse = "-"),
        "-pca.svg",
        sep = ""
      ),
      height = 7,
      width = 7
    )
    plot_pca(
      pca,
      draw_lines = list(group),
      labels = labels,
      palette = palette
    )
    dev.off()
    pdf(
      paste(
        output_prefix,
        "-",
        paste(group, sep = "-", collapse = "-"),
        "-pca.pdf",
        sep = ""
      )
    )
    plot_pca(
      pca,
      draw_lines = list(group),
      labels = labels,
      palette = palette
    )
    dev.off()
    png(
      paste(
        output_prefix,
        "-",
        paste(group, sep = "-", collapse = "-"),
        "-pca.png",
        sep = ""
      )
    )
    plot_pca(
      pca,
      draw_lines = list(group),
      labels = labels,
      palette = palette
    )
    dev.off()
  }
}

#' @title generate UMAP plots
#'
#' @description plot UMAP of the read counts
#'
#' @param counts matrix of read counts
#' @param output_prefix prefix for output files
#' @param sample character vector indicating sample, if NULL, it will be
#'   inferred from the read count matrix
#' @param treatment character vector indicating treatment, if NULL, it will be
#'   inferred from the read count matrix
#' @param treatment_groups list of treatment groups
#' @param labels if TRUE, text labels will be added to points in all plots
#' @param cores number of cores to use
#' @export
generate_umap_plots <- function(
  counts,
  output_prefix,
  sample = NULL,
  treatment = NULL,
  treatment_groups = list(),
  labels = FALSE,
  n_neighbors = 15,
  metric = "euclidean",
  n_pc = NULL,
  cores = 1
) {
  if (is.null(sample)) sample <- extract_sample_vector(counts)
  if (is.null(treatment)) treatment <- extract_treatment_vector(counts)
  if (is.null(n_pc)) {
    u <- umap(
      t(counts),
      n_neighbors = min(n_neighbors, ncol(counts) - 1),
      metric = metric,
      n_threads = cores
    )
  } else {
    pca <- prcomp(counts, rank = n_pc)
    u <- umap(
      pca[["rotation"]],
      n_neighbors = min(n_neighbors, ncol(counts) - 1),
      metric = metric,
      n_threads = cores
    )
  }
  rownames(u) <- paste(sample, treatment, sep = ".")
  pdf(paste(output_prefix, "-umap.pdf", sep = ""), height = 14)
  plot_umap(u, labels = labels)
  dev.off()
  png(paste(output_prefix, "-umap.png", sep = ""), height = 960)
  plot_umap(u, labels = labels)
  dev.off()
  for (group in treatment_groups) {
    u <- umap(
      t(counts[,treatment %in% group]),
      n_threads = cores,
      n_neighbors = min(15, sum(treatment %in% group) - 1)
    )
    rownames(u) <- paste(sample, treatment, sep = ".")[treatment %in% group]
    pdf(
      paste(
        output_prefix,
        "-",
        paste(group, sep = "-", collapse = "-"),
        "-umap.pdf",
        sep = ""
      ),
      height = 14
    )
    plot_umap(u, draw_lines = list(group), labels = labels)
    dev.off()
    png(
      paste(
        output_prefix,
        "-",
        paste(group, sep = "-", collapse = "-"),
        "-umap.png",
        sep = ""
      ),
      height = 960
    )
    plot_umap(u, draw_lines = list(group), labels = labels)
    dev.off()
  }
}

#' @title Explore ATAC-seq data
#'
#' @description perform an exploratory analysis
#'
#' @details
#' The analysis consists of the following steps:
#'
#' @param input_data list contatining data details or character path to a JSON
#'   file
#' @param output_prefix character, a prefix for output files
#' @param n_peaks integer, number of the top most variable peaks to include in
#'   the analysis
#' @param treatment_groups list providing groups of treatments to be compared
#' @param labels if TRUE, text labels will be added to points in all plots
#' @param write_counts logical, if TRUE the read count matrix will be written
#'   to disk as a TSV file
#' @param cores integer, max number of cores to use
#' @param pca logical, if true perform PCA analysis
#' @param umap logical, if true perform UMAP analysis
#' @param palette_order ordering of color palette, eithier "categorical" or
#'   "sequential"
#' @param palette the color palette
#' @export
explore <- function(
  input_data,
  output_prefix,
  treatment_groups = list(),
  labels = FALSE,
  n_peaks = 1e5,
  n_neighbors = 15,
  metric = "euclidean",
  n_pc = NULL,
  write_counts = FALSE,
  cores = 1,
  pca = TRUE,
  umap = FALSE,
  palette_order = "categorical",
  palette = NULL
) {
  preprocessed_data <- preprocess(
    input_data,
    n_peaks = n_peaks,
    cores = cores
  )
  cat(
    "median peak length: ",
    preprocessed_data[["median_peak_length"]],
    "\n",
    sep = "",
    file = paste(output_prefix, ".txt", sep = "")
  )
  if (write_counts) {
    write.table(
      preprocessed_data[["counts"]],
      file = paste(output_prefix, ".tsv", sep = ""),
      quote = FALSE,
      sep = "\t",
      row.names = FALSE
    )
  }
  if (pca) {
    generate_pca_plots(
      preprocessed_data[["counts"]],
      output_prefix,
      labels = labels,
      treatment_groups = treatment_groups,
      palette_order = palette_order,
      palette = palette
    )
  }
  if (umap) {
    generate_umap_plots(
      preprocessed_data[["counts"]],
      output_prefix,
      treatment_groups = treatment_groups,
      n_neighbors = n_neighbors,
      metric = metric,
      n_pc = n_pc,
      cores = cores
    )
  }
}
