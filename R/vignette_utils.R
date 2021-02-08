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
