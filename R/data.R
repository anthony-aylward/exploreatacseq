#' Table mapping SCATEData accession numbers to sample type
#'
#' The bioconductor package SCATEData includes 20 BAM files containing
#' single-cell ATAC-seq data. Each dataset is named by its NCBI Sequence Read
#' Archive accession number, "SRR1779746.bam". This table lists the sample type
#' of each dataset.
#'
#' @docType data
#'
#' @usage data(SCATEData_accession_sample)
#'
#' @format A data frame with 20 rows and 2 variables:
#' \describe{
#'     \item{accession}{SRA accession number}
#'     \item{sample}{Sample type}
#' }
#'
#' @keywords datasets
#'
#' @references Zhicheng Ji, Weiqiang Zhou, Wenpin Hou, Hongkai Ji,
#' \href{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02075-3}{Single-cell ATAC-seq Signal Extraction and Enhancement with SCATE},
#' Genome Biol 21, 161 (2020)
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/sra/}{NCBI Sequence Read Archive}
"SCATEData_accession_sample"

#' Vector mapping SCATEData accession numbers to sample type
#'
#' The bioconductor package SCATEData includes 20 BAM files containing
#' single-cell ATAC-seq data. Each dataset is named by its NCBI Sequence Read
#' Archive accession number, "SRR1779746.bam". This vector maps the accession
#' number to the sample type of each dataset.
#'
#' @docType data
#'
#' @usage data(SCATEData_sample)
#'
#' @format A character vector with 20 values
#'
#' @keywords datasets
#'
#' @references Zhicheng Ji, Weiqiang Zhou, Wenpin Hou, Hongkai Ji,
#' \href{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-02075-3}{Single-cell ATAC-seq Signal Extraction and Enhancement with SCATE},
#' Genome Biol 21, 161 (2020)
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/sra/}{NCBI Sequence Read Archive}
"SCATEData_sample"