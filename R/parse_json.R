#===============================================================================
# parse_json.R
#===============================================================================

# Imports ======================================================================

#' @importFrom jsonlite fromJSON


# Functions ====================================================================

#' @title extract peaks paths by sample
#'
#' @description extract peaks paths from input
#'
#' @param input_list list derived from input JSON
#' @return list of peaks file paths by sample
extract_peaks_paths_by_sample <- function(input_list) {
    lapply(
        input_list,
        function(treatment_list) unlist(
            lapply(
                treatment_list,
                function(treatment) {
                    if (sort(names(treatment)) == c("peaks", "reads", "tssenrich")) { 
                        treatment[["peaks"]] 
                    } else {
                        vapply(treatment, function(rep) rep[["peaks"]], character(length = 1))
                    }
                }
            )
        )
    )
}

#' @title extract reads file paths
#'
#' @description extract reads file paths from input
#'
#' @param input_list list derived from input JSON
#' @return character vector of reads file paths
extract_reads_file_paths <- function(input_list) {
    unlist(
        lapply(
            input_list,
            function(treatment_list) lapply(
                treatment_list,
                function(treatment) {
                    if (sort(names(treatment)) == c("peaks", "reads", "tssenrich")) {
                        treatment[["reads"]]
                    } else {
                        vapply(treatment, function(rep) rep[["reads"]], character(length = 1))
                    }
                }
            )
        )
    )
}

#' @title extract TSS enrichments
#'
#' @description extract TSS enrichment vector from input
#'
#' @param input_list list derived from input JSON
#' @return numeric vector of TSS enrichments
extract_tss_enrichments <- function(input_list) {
    unlist(
        lapply(
            input_list,
            function(treatment_list) lapply(
                treatment_list,
                function(treatment) {
                    if (sort(names(treatment)) == c("peaks", "reads", "tssenrich")) {
                        treatment[["tssenrich"]]
                    } else {
                        vapply(treatment, function(rep) rep[["tssenrich"]], double(length = 1))
                    }
                }
            )
        )
    )
}

#' @title extract group
#'
#' @description extract group vector from input
#'
#' @param input_list list derived from input JSON
#' @return character vector of groups
extract_group_vector <- function(input_list) {
    unlist(
        lapply(
            input_list,
            function(treatment_list) lapply(
                names(treatment_list),
                function(name) rep(name, length(unlist(treatment_list[[name]]))/3)
            )
        )
    )
}

#' @title extract batch
#'
#' @description extract batch vector from input
#'
#' @param input_list list derived from input JSON
#' @return character vector of batch
extract_batch_vector <- function(input_list) {
    unlist(
        lapply(
            names(input_list),
            function(name) rep(name, length(unlist(input_list[[name]]))/3)
        )
    )
}

#' @title parse input data
#'
#' @description parse input data into useful vectors
#'
#' @details
#' The list returned by fromJSON has an awkward structure. This function
#' re-organizes the input data into a list of useful vectors.
#'
#' @param input_list list derived from input JSON
#' @return list of useful vectors
parse_input_data <- function(input_list) {
    list(
        peaks_paths_by_sample = extract_peaks_paths_by_sample(input_list),
        reads_file_paths = extract_reads_file_paths(input_list),
        group = extract_group_vector(input_list),
        batch = extract_batch_vector(input_list),
        tss_enrichment = extract_tss_enrichments(input_list)
    )
}

#' @title parse input JSON
#'
#' @description parse input JSON into useful vectors using `parse_input_data()`.
#'
#' @param json_file_path path to input JSON file
#' @return list of useful vectors
parse_json <- function(json_file_path) {
    parse_input_data(fromJSON(json_file_path))
}
