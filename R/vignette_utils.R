#' @importFrom stats setNames
#' @importFrom stringr str_replace

#' @title locate SCATEData files
#'
#' @description locate BAM files from the SCATEData package
#'
#' @return list of paths to SCATEData BAM files
#' @export
locate_SCATEData_files <- function() {
    bam_list <- list.files(
        paste0(system.file(package = "SCATEData"), "/extdata/"),
        full.names = TRUE,
        pattern = ".bam$"
    )
    bam_list <- bam_list[grepl(".bam$", bam_list)]
    names(bam_list) <- str_replace(basename(bam_list), ".bam$", "")
    bam_list
}

#' @title map SCATEData file paths to sample types
#'
#' @description map SCATEData file paths to sample types
#'
#' @return list mapping file paths to sample types, sorted by sample types
#' @export
map_paths_to_samples <- function() {
    data(SCATEData_sample)
    SCATEData_path <- locate_SCATEData_files()
    samples <- vapply(
        names(SCATEData_path),
        function(name) SCATEData_sample[[name]],
        character(length = 1)
    )
    sort(setNames(samples, unname(SCATEData_path)))
}


