#' UpSet plot for gene sets
#'
#' Convenience wrapper around \code{aplotExtra::upset_plot()} so the plot is
#' directly available from TCMDATA. This is particularly useful for the output
#' of [get_ml_gene_sets()].
#'
#' @param list A named list of gene sets.
#' @param nintersects Number of intersections to show. \code{NULL} shows all.
#' @param order.intersect.by One of \code{"size"} or \code{"name"}.
#' @param order.set.by One of \code{"size"} or \code{"name"}.
#' @param color.intersect.by Color scheme for intersection bars.
#' @param color.set.by Color scheme for set bars.
#' @param remove_empty_intersects Whether to remove empty intersections.
#'
#' @return An UpSet plot object returned by \code{aplotExtra}.
#'
#' @examples
#' \dontrun{
#'   gene_sets <- list(
#'     LASSO = c("TP53", "BRCA1", "EGFR"),
#'     RF = c("TP53", "AKT1"),
#'     XGBOOST = c("EGFR", "AKT1", "MTOR")
#'   )
#'   upsetplot(gene_sets)
#' }
#' @export
upsetplot <- function(list,
                      nintersects = NULL,
                      order.intersect.by = c("size", "name"),
                      order.set.by = c("size", "name"),
                      color.intersect.by = "none",
                      color.set.by = "none",
                      remove_empty_intersects = TRUE) {
  if (!requireNamespace("aplotExtra", quietly = TRUE)) {
    stop("Package 'aplotExtra' is required for upsetplot(). Please install it first.",
         call. = FALSE)
  }
  if (!is.list(list)) {
    stop("'list' must be a named list of sets.", call. = FALSE)
  }
  if (length(list) < 2) {
    stop("Need at least 2 sets for upsetplot().", call. = FALSE)
  }

  nm <- names(list)
  if (is.null(nm)) {
    names(list) <- paste("Set", seq_along(list), sep = "_")
  } else if (any(!nzchar(nm))) {
    names(list)[!nzchar(nm)] <- paste("Set", which(!nzchar(nm)), sep = "_")
  }

  aplotExtra::upset_plot(
    list = list,
    nintersects = nintersects,
    order.intersect.by = match.arg(order.intersect.by),
    order.set.by = match.arg(order.set.by),
    color.intersect.by = color.intersect.by,
    color.set.by = color.set.by,
    remove_empty_intersects = remove_empty_intersects
  )
}
