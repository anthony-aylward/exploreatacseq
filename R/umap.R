#===============================================================================
# umap.R
#===============================================================================

#' @import RColorBrewer


# Functions ====================================================================

#' @title plot umap coordinates
#'
#' @description generate a UMAP plot
#'
#' @param umap_matrix matrix of umap coordinates with appropriate rownames
#' @param draw_lines list of treatment groups to draw lines through
#' @param palette color palette for points
#' @param labels if TRUE, text labels will be added to the points
#' @export
plot_umap <- function(
  umap_matrix,
  draw_lines = list(),
  labels = FALSE,
  palette = exploreatacseq_color_palette()
) {
  coord_by_treat <- coordinates_by_treatment(umap_matrix)

  draw_line <- function(sample, start_treatment, end_treatment, rep = NULL) {
    if (is.null(rep)) {
      lines(
        c(
          coord_by_treat[[start_treatment]][paste(sample, start_treatment, sep = "."), 1],
          coord_by_treat[[end_treatment]][paste(sample, end_treatment, sep = "."), 1]
        ),
        c(
          coord_by_treat[[start_treatment]][paste(sample, start_treatment, sep = "."), 2],
          coord_by_treat[[end_treatment]][paste(sample, end_treatment, sep = "."), 2]
        ),
        lwd = 4,
        col = "lightgray"
      )
    } else {
      lines(
        c(
          coord_by_treat[[start_treatment]][paste(sample, start_treatment, rep, sep = "."), 1],
          coord_by_treat[[end_treatment]][paste(sample, end_treatment, rep, sep = "."), 1]
        ),
        c(
          coord_by_treat[[start_treatment]][paste(sample, start_treatment, rep, sep = "."), 2],
          coord_by_treat[[end_treatment]][paste(sample, end_treatment, rep, sep = "."), 2]
        ),
        lwd = 4,
        col = "lightgray"
      )
    }
  }

  layout(matrix(c(3, 4, 1, 2), 2, 2, byrow = FALSE), widths = c(1, 2), heights = c(2, 1))
  par(mai = c(0.65, 0.65, 0.1, 0.1), omi = c(0.1, 0.1, 0.1, 0.1))

  plot(umap_matrix[,1], umap_matrix[,2], col = "white", xaxt = "n", yaxt = "n", ann = FALSE)
  for (group in draw_lines) {
    for (i in 1:(length(group) - 1)) {
      start_treatment = group[i]
      end_treatment = group[i + 1]
      start_samples = lapply(
        strsplit(rownames(coord_by_treat[[start_treatment]]), split = ".", fixed = TRUE),
        function(x) if (length(x) == 3) x[c(1,3)] else x[[1]]
      )
      end_samples = lapply(
        strsplit(rownames(coord_by_treat[[end_treatment]]), split = ".", fixed = TRUE),
        function(x) if (length(x) == 3) x[c(1,3)] else x[[1]]
      )
      samples = intersect(start_samples, end_samples)
      for (sample in samples) {
        if (length(sample) == 2) {
          draw_line(sample[[1]], start_treatment, end_treatment, rep = sample[[2]])
        } else {
          draw_line(sample, start_treatment, end_treatment)
        }
      }
      if (length(group) > i + 1) {
        bridge_treatment = group[i + 2]
        bridge_samples = lapply(
          strsplit(rownames(coord_by_treat[[bridge_treatment]]), split = ".", fixed = TRUE),
          function(x) if (length(x) == 3) x[c(1,3)] else x[[1]]
        )
        samples = setdiff(intersect(start_samples, bridge_samples), samples)
        for (sample in samples) {
          if (length(sample) == 2) {
            draw_line(sample[[1]], start_treatment, bridge_treatment, rep = sample[[2]])
          } else {
            draw_line(sample, start_treatment, bridge_treatment)
          }
        }
      }
    }
  }

  n_treatments <- length(coord_by_treat)
  for (i in 1:n_treatments) {
    coord <- coord_by_treat[[i]]
    points(coord[,1], coord[,2], col = palette[[i]], pch = 19, cex = 2)
    if (labels) text(coord[,1], coord[,2], labels = rownames(coord), pos = 1)
  }
  grp <- unlist(lapply(names(coord_by_treat), function(x) rep(x, nrow(coord_by_treat[[x]]))))
  pc <- lapply(c(1, 2), function(y) unlist(lapply(coord_by_treat, function(x) as.numeric(x[,y]))))
  box_colors <- palette[1:n_treatments][order(vapply(coord_by_treat, function(x) median(x[,1])))]
  by_median <- reorder(grp, pc[[1]], median)
  
  boxplot(pc[[1]] ~ by_median, horizontal = TRUE, las = 1, col = box_colors, yaxt="n", ann=FALSE)
  title(xlab = paste("PC1 [", percent_of_variance[["PC1"]], "%]", sep = ""))
  
  boxplot(pc[[2]] ~ by_median, col =  box_colors, xaxt="n", ann=FALSE)
  title(ylab = paste("PC2 [", percent_of_variance[["PC2"]], "%]", sep = ""))

  plot(0:1, 0:1, col = "white", xaxt = "n", yaxt = "n", bty = "n", ann = FALSE)
  legend(0, 1, legend = names(coord_by_treat), col = palette[1:n_treatments], pch = 19)
}
