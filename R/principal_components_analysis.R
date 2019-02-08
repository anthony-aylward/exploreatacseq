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
  lapply(unique(treatment), function(x) pca[treatment == x,, drop = FALSE])
}

#' @title plot principal components
#'
#' @description generate a PCA plot
#'
#' @param pca two-column matrix of principal component coordinates
#' @param draw_lines list of treatment groups to draw lines through
#' @export
plot_pca <- function(pca, draw_lines = list()) {
  coord_by_treat <- coordinates_by_treatment(pca)
  palette <- brewer.pal(9, "Set1")[c(2, 1, 3:5, 7:9)]
  plot(pca[,1], pca[,2], col = "white", xlab = "PC1", ylab = "PC2")

  for (group in draw_ines) {
    for (index in 1:length(group) - 1) {
      start_treatment = group[index]
      end_treatment = group[index + 1]
      start_samples = strsplit(
        rownames(coord_by_treat[[start_treatment]]),
        split = ".",
        fixed = TRUE
      )
      end_samples = strsplit(
        rownames(coord_by_treat[[start_treatment]]),
        split = ".",
        fixed = TRUE
      )
    }
  }

  for (i in 1:length(coord_by_treat)) {
    coord <- coord_by_treat[[i]]
    points(coord[,1], coord[,2], col = palette[[i]], pch = 19)
  }
}