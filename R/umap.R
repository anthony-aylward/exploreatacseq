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
#' @export
plot_umap <- function(umap_matrix, draw_lines = list()) {
  coord_by_treat <- coordinates_by_treatment(umap_matrix)

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
  plot(umap_matrix[,1], umap_matrix[,2], col = "white", ann = FALSE)
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
    text(
      coord[,1],
      coord[,2],
      labels = sapply(
        strsplit(rownames(coord), split = ".", fixed = TRUE),
        function(x) x[[1]]
      ),
      pos = 1
    )
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
