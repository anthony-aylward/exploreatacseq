#===============================================================================
# principal_components_analysis.R
#===============================================================================

# Imports ======================================================================

#' @import RColorBrewer




# Functions ====================================================================

#' @title two principal components
#'
#' @description compute two principal components for the count matrix
#'
#' @param count_matrix matrix of read counts
#' @return matrix with two columns giving the principal component coordinates
#' @export
two_principal_components <- function(count_matrix) {
  prcomp(count_matrix, rank = 2)[["rotation"]]
}

#' @title coordinates by treatment
#'
#' @description organize PCA coordinates by treatment
#'
#' @param pca two-column matrix of principal component coordinates
#' @return list of two-column matrices, one per treatment
#' @export
coordinates_by_treatment <- function(pca) {
  treatment <- sapply(
    strsplit(rownames(pca), split = ".", fixed = TRUE),
    function(x) x[[2]]
  )
  lapply(unique(treatment), function(x) pca[treatment == x,])
}

#' @title plot principal components
#'
#' @description generate a PCA plot
#'
#' @param pca two-column matrix of principal component coordinates
#' @export
plot_pca <- function(pca) {
  plot(pca[[1]], pca[[2]], col = "white")
  coord_by_treat <- coordinates_by_treatment(pca)
  n_treatments <- length(coord_by_treat)
  palette_size <- min(3, n_treatments)
  palette <- brewer.pal(palette_size, "Set1")[c(2, 1, 3:palette_size)]
  plot(pca, col = "white")
  for (i in n_treatments) {
    points(coord_by_treat[[i]], col = palette[[i]], pch = 19)
  }
}