#===============================================================================
# umap.R
#===============================================================================

# Imports ======================================================================

#' @import uwot
#' @import vizier




# Functions ====================================================================

#' @title preprocess count data for UMAP
#'
#' @description take the transpose of the count data and add a treatment column
#'
#' @param x count matrix
#' @param treatment character vector indicating treatment for each experiment
#' @return matrix for input to `umap()`
#' @export
preprocess_for_umap <- function(x, treatment) {
  cbind(as.data.frame(t(x)), treatment)
}

#' @title reduce dimensionality with umap
#'
#' @description apply the `umap()` function from `uwot`
#'
#' @param x input matrix or data.frame
#' @param n_threads number of threads to use
#' @return the result of the `umap()` call
#' @export
umap_reduce <- function(x, n_threads = 1) {
  umap(x, n_threads = n_threads, n_sgd_threads = n_threads, approx_pow = TRUE)
}

#' @title embed image
#'
#' @description plot an embedding from umap
#'
#' @param X input matrix
#' @param Y input umap
#' @export
embed_img <- function(X, Y, k = 15, ...) {
  args <- list(...)
  args[["coords"]] <- Y
  args[["x"]] <- X

  do.call(embed_plot, args)
}

#' @title plot a umap
#'
#' @description plot an embedding from umap
#'
#' @param x input matrix
#' @param x_umap input umap
#' @export
plot_umap <- function(x, x_umap) {
  embed_img(
    x,
    x_umap,
    pc_axes = TRUE,
    equal_axes = TRUE,
    alpha_scale = 0.5,
    cex = 1
  )
}
