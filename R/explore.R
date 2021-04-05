#===============================================================================
# explore.R
#===============================================================================

#' @importFrom grDevices dev.off pdf png
#' @importFrom stats prcomp setNames
#' @importFrom svglite svglite
#' @importFrom utils combn write.table




# Functions ====================================================================

#' @title preprocess ATAC-seq data
#'
#' @description preprocess multiple ATAC-seq datasets into a read count matrix
#' @param input_data list contatining data details or character path to a JSON
#'     file
#' @param n_peaks integer, number of the top most variable peaks to include in
#'     the analysis
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
    raw_counts <- get_raw_counts(peaks, x[["reads_file_paths"]], cores = cores)
    transformed_counts <- most_variable_peaks(
        transform_counts(
            count_dgelist(raw_counts, x[["group"]]),
            batch = x[["batch"]],
            covariates = x[["tss_enrichment"]]
        ),
        n = n_peaks
    )
    list(
        median_peak_length = median_peak_length(peaks),
        raw_counts = raw_counts,
        transformed_counts = transformed_counts
    )
}

#' @title generate PCA plots
#'
#' @description plot PCA of the read counts
#'
#' @param counts matrix of read counts
#' @param output_prefix prefix for output files
#' @param treatment character vector indicating treatment, if NULL, it will be
#'     inferred from the read count matrix
#' @param treatment_order order of treatments WRT the palette
#' @param treatment_groups list of treatment groups
#' @param labels if TRUE, text labels will be added to points in all plots
#' @param palette_order ordering of color palette, eithier "categorical" or
#'     "sequential"
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
            exploreatacseq_color_palette(order = palette_order)[seq_len(length(treatment_order))],
            treatment_order
        )
    } else {
        palette_vector = setNames(
            palette[seq_len(length(treatment_order))],
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
            paste(output_prefix, "-", paste(group, sep = "-", collapse = "-"), "-pca.svg", sep = ""),
            height = 7,
            width = 7
        )
        plot_pca(pca, draw_lines = list(group), labels = labels, palette = palette)
        dev.off()
        pdf(paste(output_prefix, "-", paste(group, sep = "-", collapse = "-"), "-pca.pdf", sep = ""))
        plot_pca(pca, draw_lines = list(group), labels = labels, palette = palette)
        dev.off()
        png(paste(output_prefix, "-", paste(group, sep = "-", collapse = "-"), "-pca.png", sep = ""))
        plot_pca(pca, draw_lines = list(group), labels = labels, palette = palette)
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
#'     file
#' @param output_prefix character, a prefix for output files
#' @param n_peaks integer, number of the top most variable peaks to include in
#'     the analysis
#' @param treatment_groups list providing groups of treatments to be compared
#' @param labels if TRUE, text labels will be added to points in all plots
#' @param write_raw_counts logical, if TRUE the raw read count matrix will be
#'     written to disk as a TSV file
#' @param write_transformed_counts logical, if TRUE the transformed read count
#'     matrix will be written to disk as a TSV file
#' @param write_counts logical, if TRUE both raw and transformed read count
#'     matrices will be written to disk as TSV files
#' @param cores integer, max number of cores to use
#' @param palette_order ordering of color palette, eithier "categorical" or
#'     "sequential"
#' @param palette the color palette
#' @export
explore <- function(
    input_data,
    output_prefix,
    treatment_groups = list(),
    labels = FALSE,
    n_peaks = 1e5,
    write_raw_counts = FALSE,
    write_transformed_counts = FALSE,
    write_counts = FALSE,
    cores = 1,
    palette_order = "categorical",
    palette = NULL
) {
    preprocessed_data <- preprocess(input_data, n_peaks = n_peaks, cores = cores)
    cat(
        "median peak length: ",
        preprocessed_data[["median_peak_length"]],
        "\n",
        sep = "",
        file = paste(output_prefix, ".txt", sep = "")
    )
    if (write_counts || write_raw_counts) {
        write.table(
            preprocessed_data[["raw_counts"]],
            file = paste(output_prefix, ".raw.tsv", sep = ""),
            quote = FALSE,
            sep = "\t"
        )
    }
    if (write_counts || write_transformed_counts) {
        write.table(
            preprocessed_data[["transformed_counts"]],
            file = paste(output_prefix, ".transformed.tsv", sep = ""),
            quote = FALSE,
            sep = "\t"
        )
    }
    generate_pca_plots(
        preprocessed_data[["transformed_counts"]],
        output_prefix,
        labels = labels,
        treatment_groups = treatment_groups,
        palette_order = palette_order,
        palette = palette
    )
}
