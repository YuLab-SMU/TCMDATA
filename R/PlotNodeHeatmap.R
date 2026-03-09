#' Plot PPI nodes metrics Heatmap
#'
#' Create a Z-score normalized heatmap of PPI node topological metrics using
#' \code{ComplexHeatmap}.  Provides convenient parameters for the most common
#' styling adjustments while still accepting any additional argument via
#' \code{...}.
#'
#' @param data A data frame containing the node names and metric values.
#' @param id_col Character. The column name serving as row identifiers.
#'   Default is \code{"name"}.
#' @param select_cols Character vector. Columns to include in the heatmap.
#'   If \code{NULL} (default), all numeric columns except \code{id_col} are
#'   used.
#' @param colors Character vector of length 3. Colours mapped to Z-score
#'   values \code{-2}, \code{0}, and \code{2}. Default is
#'   \code{c("#2166AC", "white", "#B2182B")}.
#' @param cluster_rows Logical. Perform hierarchical clustering on rows?
#'   Default \code{TRUE}.
#' @param cluster_cols Logical. Perform hierarchical clustering on columns?
#'   Default \code{FALSE}.
#' @param row_fontsize Numeric. Font size for row (gene) names.
#'   Default \code{12}.
#' @param row_fontface Character. Font face for row names. One of
#'   \code{"plain"}, \code{"italic"}, \code{"bold"}, or
#'   \code{"bold.italic"}. Default \code{"italic"}.
#' @param col_fontsize Numeric. Font size for column (metric) names.
#'   Default \code{10}.
#' @param col_fontface Character. Font face for column names. Default
#'   \code{"bold"}.
#' @param col_rotation Numeric. Rotation angle (degrees) for column labels.
#'   Default \code{45}.
#' @param show_row_names Logical. Show row names? Default \code{TRUE}.
#' @param show_column_names Logical. Show column names? Default \code{TRUE}.
#' @param legend_title Character. Title for the colour legend. Default
#'   \code{"Z-score"}.
#' @param border_color Character. Cell border colour. Default \code{"white"}.
#'   Set to \code{NA} to remove borders.
#' @param border_width Numeric. Cell border line width. Default \code{1}.
#' @param ... Additional arguments passed to
#'   \code{\link[ComplexHeatmap]{Heatmap}}.
#'
#' @importFrom grid gpar unit
#' @return A \code{Heatmap} object that can be printed or combined with other
#'   heatmaps via \code{+} or \code{\%v\%}.
#'
#' @examples
#' data(demo_ppi)
#' ppi <- compute_nodeinfo(demo_ppi)
#' rk_res <- rank_ppi_nodes(ppi)[[2]]
#' selected_cols <- colnames(rk_res)[c(2, 3, 4, 6, 9, 12, 14, 15, 16)]
#' p1 <- plot_node_heatmap(rk_res, select_cols = selected_cols)
#' print(p1)
#'
#' @export
plot_node_heatmap <- function(data,
                              id_col = "name",
                              select_cols = NULL,
                              colors = c("#2166AC", "white", "#B2182B"),
                              cluster_rows = TRUE,
                              cluster_cols = FALSE,
                              row_fontsize = 12,
                              row_fontface = "italic",
                              col_fontsize = 10,
                              col_fontface = "bold",
                              col_rotation = 45,
                              show_row_names = TRUE,
                              show_column_names = TRUE,
                              legend_title = "Z-score",
                              border_color = "white",
                              border_width = 1,
                              ...) {

  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required. ",
         "Install via: BiocManager::install(\"ComplexHeatmap\")")
  }

  if (!id_col %in% colnames(data)) {
    stop("Column '", id_col, "' not found in the input data.")
  }

  plot_data <- data
  rownames(plot_data) <- plot_data[[id_col]]
  plot_data[[id_col]] <- NULL

  if (!is.null(select_cols)) {
    missing <- setdiff(select_cols, colnames(plot_data))
    if (length(missing) > 0) {
      stop("Columns not found: ", paste(missing, collapse = ", "))
    }
    plot_data <- plot_data[, select_cols, drop = FALSE]
  }

  is_num <- vapply(plot_data, is.numeric, logical(1))
  plot_data <- plot_data[, is_num, drop = FALSE]

  if (ncol(plot_data) == 0) {
    stop("No numeric columns left to plot after selection.")
  }

  mat_scaled <- scale(plot_data)
  mat_scaled[is.na(mat_scaled)] <- 0

  col_fun <- circlize::colorRamp2(c(-2, 0, 2), colors)

  rect_gp <- if (is.na(border_color)) {
    gpar(col = NA)
  } else {
    gpar(col = border_color, lwd = border_width)
  }

  p <- ComplexHeatmap::Heatmap(
    matrix            = mat_scaled,
    name              = legend_title,
    col               = col_fun,
    cluster_rows      = cluster_rows,
    cluster_columns   = cluster_cols,
    show_row_names    = show_row_names,
    show_column_names = show_column_names,
    rect_gp           = rect_gp,
    row_names_gp      = gpar(fontsize = row_fontsize, fontface = row_fontface),
    column_names_gp   = gpar(fontsize = col_fontsize, fontface = col_fontface),
    column_names_rot  = col_rotation,
    row_dend_width    = unit(2, "cm"),
    ...
  )

  return(p)
}

#' @rdname plot_node_heatmap
#' @export
PlotNodeHeatmap <- function(...) {
  warning("PlotNodeHeatmap is deprecated. Please use plot_node_heatmap instead.")
  plot_node_heatmap(...)
}
