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
  prcomp(count_matrix, rank = 2)
}

#' @title coordinates by treatment
#'
#' @description organize PCA coordinates by treatment
#'
#' @param coord two-column matrix of principal component coordinates
#' @return list of two-column matrices, one per treatment
#' @export
coordinates_by_treatment <- function(coord) {
  treatment <- sapply(
    strsplit(rownames(coord), split = ".", fixed = TRUE),
    function(x) x[[2]]
  )
  treat_uniq <- unique(treatment)
  setNames(
    lapply(treat_uniq, function(x) coord[treatment == x,, drop = FALSE]),
    treat_uniq
  )
}

#' @title plot principal components
#'
#' @description generate a PCA plot
#'
#' @param pca list with class "prcomp"
#' @param draw_lines list of treatment groups to draw lines through
#' @export
plot_pca <- function(pca, draw_lines = list()) {
  coord_by_treat <- coordinates_by_treatment(pca[["rotation"]])
  percent_of_variance <- round(summary(pca)[["importance"]][2, 1:2] * 100)

  draw_line <- function(sample, start_treatment, end_treatment) {
    lines(
      c(
        coord_by_treat[[start_treatment]][
          paste(sample, start_treatment, sep = "."), 1
        ],
        coord_by_treat[[end_treatment]][
          paste(sample, end_treatment, sep = "."), 1
        ]
      ),
      c(
        coord_by_treat[[start_treatment]][
          paste(sample, start_treatment, sep = "."), 2
        ],
        coord_by_treat[[end_treatment]][
          paste(sample, end_treatment, sep = "."), 2
        ]
      ),
      lwd = 4,
      col = "lightgray"
    )
  }

  palette <- brewer.pal(9, "Set1")[c(2, 1, 3:5, 7:9)]
  par(mfcol = c(2, 1))
  plot(
    pca[["rotation"]][,1],
    pca[["rotation"]][,2],
    col = "white",
    xlab = paste("PC1 [", percent_of_variance[["PC1"]], "%]", sep = ""),
    ylab = paste("PC2 [", percent_of_variance[["PC2"]], "%]", sep = "")
  )
  for (group in draw_lines) {
    for (i in 1:(length(group) - 1)) {
      start_treatment = group[i]
      end_treatment = group[i + 1]
      start_samples = sapply(
        strsplit(
          rownames(coord_by_treat[[start_treatment]]),
          split = ".",
          fixed = TRUE
        ),
        function(x) x[[1]]
      )
      end_samples = sapply(
        strsplit(
          rownames(coord_by_treat[[end_treatment]]),
          split = ".",
          fixed = TRUE
        ),
        function(x) x[[1]]
      )
      samples = intersect(start_samples, end_samples)
      for (sample in samples) {
        draw_line(sample, start_treatment, end_treatment)
      }
      if (length(group) > i + 1) {
        bridge_treatment = group[i + 2]
        bridge_samples = sapply(
          strsplit(
            rownames(coord_by_treat[[bridge_treatment]]),
            split = ".",
            fixed = TRUE
          ),
          function(x) x[[1]]
        )
        samples = setdiff(intersect(start_samples, bridge_samples), samples)
        for (sample in samples) {
          draw_line(sample, start_treatment, bridge_treatment)
        }
      }
    }
  }
  n_treatments <- length(coord_by_treat)
  for (i in 1:n_treatments) {
    coord <- coord_by_treat[[i]]
    points(coord[,1], coord[,2], col = palette[[i]], pch = 19, cex = 2)
    text(coord[,1], coord[,2], labels = rownames(coord))
  }
  plot(0:1, 0:1, col = "white", xaxt = "n", yaxt = "n", bty = "n", ann = FALSE)
  legend(
    0,
    1,
    legend = names(coord_by_treat),
    col = palette[1:n_treatments],
    pch = 19,
    cex = 2
  )
}
