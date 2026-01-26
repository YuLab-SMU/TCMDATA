#' Plot PPI nodes metrics Heatmap
#'
#' @param data A data frame containing the nodes names and metric values.
#' @param id_col Character. The column name serving as row names. Default is "name".
#' @param select_cols Character Vector. The specific columns to include in the heatmap. Default all columns except `id_col` are used.
#' @param colors Vector of length 3. Colors for low (-2), zero (0), and high (2) values. Default is c("#2166AC", "white", "#B2182B").
#' @param cluster_rows Logical. Whether to perform hierarchical clustering on rows. Default is TRUE.
#' @param cluster_cols Logical. Whether to perform hierarchical clustering on columns. Default is FALSE.
#' @param ... Additional arguments passed to `ComplexHeatmap::Heatmap`.
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar unit
#' @return A Heatmap object.
#' @examples
#' data(demo_ppi)
#' ppi <- compute_nodeinfo(demo_ppi)
#' rk_res <- rank_ppi_nodes(ppi)[[2]]
#' selected_cols <- colnames(rk_res)[c(2, 3, 4, 6, 9, 12, 14, 15, 16)]
#' p1 <- PlotNodeHeatmap(rk_res, select_cols = selected_cols)
#' print(p1)
#' @export
PlotNodeHeatmap <- function(data, 
                            id_col = "name", 
                            select_cols = NULL, 
                            colors = c("#2166AC", "white", "#B2182B"), 
                            cluster_rows = TRUE, 
                            cluster_cols = FALSE,
                            ...) {
  
  # Check for ComplexHeatmap availability
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("Package 'ComplexHeatmap' is required but not installed. Please install it from Bioconductor using: BiocManager::install(\"ComplexHeatmap\")")
  }
  
  if (!id_col %in% colnames(data)) {
    stop(paste("Error: The column", id_col, "does not exist in the input data."))
  }
  
  plot_data <- data
  rownames(plot_data) <- plot_data[[id_col]]
  plot_data[[id_col]] <- NULL 
  
  if (!is.null(select_cols)) {
    missing <- setdiff(select_cols, colnames(plot_data))
    if (length(missing) > 0) {
      stop(paste("Error: The following columns were not found:", paste(missing, collapse=", ")))
    }
    plot_data <- plot_data[, select_cols, drop = FALSE]
  }
  
  is_num <- sapply(plot_data, is.numeric)
  plot_data <- plot_data[, is_num, drop = FALSE]
  
  if (ncol(plot_data) == 0) {
    stop("Error: No numeric columns left to plot after selection.")
  }
  
  # scale data and deal with NA
  mat_scaled <- scale(plot_data)
  mat_scaled[is.na(mat_scaled)] <- 0
  
  # Color Mapping
  col_fun <- circlize::colorRamp2(c(-2, 0, 2), colors)
  
  p <- ComplexHeatmap::Heatmap(matrix = mat_scaled,
                               name = "Z-score",
                               col = col_fun,
                               cluster_rows = cluster_rows,
                               cluster_columns = cluster_cols,
                               rect_gp = gpar(col = "white", lwd = 1), 
                               row_names_gp = gpar(fontsize = 12),     
                               column_names_gp = gpar(fontsize = 10, fontface = "italic"), 
                               column_names_rot = 45,
                               row_dend_width = unit(2, "cm"),
                               ...)
  
  return(p)
}
