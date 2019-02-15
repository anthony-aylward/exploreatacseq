#===============================================================================
# explore.R
#===============================================================================

#' @import uwot




# Functions ====================================================================

#' @title generate PCA plots
#'
#' @description plot PCA of the read counts
#'
#' @param counts matrix of read counts
#' @param output_prefix prefix for output files
#' @param treatment character vector indicating treatment
#' @param treatment_groups list of treatment groups
#' @export
generate_pca_plots <- function(
  counts,
  output_prefix,
  treatment,
  treatment_groups = list()
) {
  pca <- two_principal_components(counts)
  pdf(paste(output_prefix, "-pca.pdf", sep = ""), height = 14)
  plot_pca(pca)
  dev.off()
  png(paste(output_prefix, "-pca.png", sep = ""), height = 960)
  plot_pca(pca)
  dev.off()
  for (group in treatment_groups) {
    pca <- two_principal_components(counts[,treatment %in% group])
    pdf(
      paste(
        output_prefix,
        "-",
        paste(group, sep = "-", collapse = "-"),
        "-pca.pdf",
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
        "-pca.png",
        sep = ""
      ),
      height = 960
    )
    plot_pca(pca, draw_lines = list(group))
    dev.off()
  }
}

#' @title generate UMAP plots
#'
#' @description plot UMAP of the read counts
#'
#' @param counts matrix of read counts
#' @param output_prefix prefix for output files
#' @param sample character vector indicating sample
#' @param treatment character vector indicating treatment
#' @param treatment_groups list of treatment groups
#' @export
generate_umap_plots <- function(
  counts,
  output_prefix,
  sample,
  treatment,
  treatment_groups = list(),
  cores = 1
) {
  u <- umap(
    t(counts),
    n_threads = cores,
    n_neighbors = min(15, ncol(counts) - 1)
  )
  rownames(u) <- paste(sample, treatment, sep = ".")
  pdf(paste(output_prefix, "-umap.pdf", sep = ""), height = 14)
  plot_umap(u)
  dev.off()
  png(paste(output_prefix, "-umap.png", sep = ""), height = 960)
  plot_umap(u)
  dev.off()
  for (group in treatment_groups) {
    u <- umap(
      t(counts[,treatment %in% group]),
      n_threads = cores
      n_neighbors = min(15, sum(treatment %in% group) - 1)
    )
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
    plot_umap(u, draw_lines = list(group))
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
    plot_umap(u, draw_lines = list(group))
    dev.off()
  }
}

#' @title Explore ATAC-seq data
#'
#' @description perform an exploratory analysis
#'
#' @param json_file_path path to a JSON file providing data details
#' @param output_prefix character, a prefix for output files
#' @param treatment_groups list providing groups of treatments to be compared
#' @param cores integer, max number of cores to use
#' @export
explore <- function(
  json_file_path,
  output_prefix,
  treatment_groups = list(),
  cores = 1
) {
  x <- parse_json(json_file_path)
  peaks <- filter_peaks(
    consensus_peaks(peaks_by_sample(x[["peaks_paths_by_sample"]]))
  )
  cat(
    "median peak length: ",
    median_peak_length(peaks),
    "\n",
    sep = "",
    file = paste(output_prefix, ".txt", sep = "")
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
    )
  )
  write.table(
    counts,
    file = paste(output_prefix, ".tsv", sep = ""),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  treatment <- extract_treatment_vector(counts)
  generate_pca_plots(
    counts,
    output_prefix,
    treatment,
    treatment_groups = treatment_groups
  )
  generate_umap_plots(
    counts,
    output_prefix,
    extract_sample_vector(counts),
    treatment,
    treatment_groups = treatment_groups,
    cores = cores
  )
}
