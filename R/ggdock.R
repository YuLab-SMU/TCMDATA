#' Plot molecule–target docking affinity heatmaps
#'
#' @param dock_data A numeric matrix or data frame.
#' @param order Reordering method for factor levels (optional).
#' @param type Type of visualization, either `"dot"` or `"tile"`.
#' @param point_size Numeric. Point size used in dot plots. Default is 3.
#' @param base_size Numeric. Base font size for the theme. Default is 7.
#' @param angle Numeric. Rotation angle for x-axis text labels. Default is 50.
#' @param hjust,vjust Numeric. Horizontal and vertical justification for x-axis labels.
#' @param palette Character. Name of the continuous color palette to use via
#'   \pkg{paletteer}. If \code{NULL}, a built-in low-saturation publication
#'   palette is used.
#' @param label Logical. If `TRUE`, print affinity numbers on the plot. Default `FALSE`.
#' @param label_digits Integer. Number of digits to show for the affinity text. Default `2`.
#' @param label_size Numeric. Text size for affinity labels. Default `3.5`.
#' @param label_color Character. Text color for affinity labels. Default `"black"`.
#' @param label_family Character. Font family for labels (e.g., `"sans"`, `"Times"`). Default `"sans"`.
#' @param label_fontface Character. Font face for labels: `"plain"`, `"bold"`, `"italic"`, `"bold.italic"`. Default `"plain"`.
#' @param ... Additional parameters passed to `geom_point()` or `geom_tile()`.
#'
#' @return A `ggplot` object.
#'
#' @importFrom ggplot2 aes coord_fixed element_blank geom_point geom_text
#' @importFrom ggplot2 geom_tile ggplot labs scale_color_gradient2
#' @importFrom ggplot2 scale_fill_gradient2 theme
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate all_of
#' @importFrom stats median
#' @importFrom rlang .data
#' @importFrom magrittr %>%
#'
#' @export

ggdock <- function(
    dock_data,
    order = NULL,
    type = "dot",
    point_size = 3,
    base_size = 7,
    angle = 50,
    hjust = 1,
    vjust = 1,
    palette = NULL,
    label = FALSE,
    label_digits = 2,
    label_size = 3,
    label_color = "black",
    label_family = "sans",
    label_fontface = "plain",
    ...){

  if (!is.null(palette) && !requireNamespace("paletteer", quietly = TRUE)) {
    stop("Package 'paletteer' is required for ggdock(). Please install it with: install.packages('paletteer')")
  }
  if (!requireNamespace("forcats", quietly = TRUE)) {
    stop("Package 'forcats' is required for ggdock(). Please install it with: install.packages('forcats')")
  }

  type <- match.arg(type, c("dot", "tile"))

  if (is.matrix(dock_data)) dock_data <- as.data.frame(dock_data)

  ## convert data structure
  df_long <- dock_data %>%
    tibble::rownames_to_column("target") %>%
    tidyr::pivot_longer(
      cols = !dplyr::all_of("target"),
      names_to = "molecule",
      values_to = "affinity")

  ## get order function (eg. max, min, median)
  order_fun <- get_order(order)

  if (!is.null(order_fun)) {
    df_long <- df_long |>
      dplyr::mutate(
        target   = forcats::fct_reorder(.data$target,   .data$affinity,
                                        .fun = order_fun, .desc = TRUE),
        molecule = forcats::fct_reorder(.data$molecule, .data$affinity,
                                        .fun = order_fun, .desc = TRUE)
      )
  }

  if (type == "dot") {
    p <- ggplot(data = df_long, aes(x = .data[["molecule"]], y = .data[["target"]])) +
      geom_point(aes(color = .data[["affinity"]]), size = point_size, ...) +
      labs(x = NULL, y = NULL, color = "Binding affinity\n(kcal/mol)") +
      .theme_tcm_pub(base_size = base_size, grid = "both") +
      theme(
        axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust),
        axis.line = element_blank(),
        axis.ticks = element_blank()
      ) +
      coord_fixed()
    if (is.null(palette)) {
      p <- p + scale_color_gradient2(
        low = "#4E79A7", mid = "white", high = "#B85C5C",
        midpoint = stats::median(df_long$affinity, na.rm = TRUE),
        name = "Binding affinity\n(kcal/mol)"
      )
    } else {
      p <- p + paletteer::scale_color_paletteer_c(
        palette, name = "Binding affinity\n(kcal/mol)"
      )
    }
  }
  else {
    p <- ggplot(data = df_long, aes(x = .data[["molecule"]], y = .data[["target"]], fill = .data[["affinity"]])) +
      geom_tile(...) +
      labs(x = NULL, y = NULL, fill = "Binding affinity\n(kcal/mol)") +
      .theme_tcm_pub(base_size = base_size, grid = "none") +
      theme(
        axis.text.x = element_text(angle = angle, hjust = hjust, vjust = vjust),
        axis.line = element_blank(),
        axis.ticks = element_blank()
      ) +
      coord_fixed()
    if (is.null(palette)) {
      p <- p + scale_fill_gradient2(
        low = "#4E79A7", mid = "white", high = "#B85C5C",
        midpoint = stats::median(df_long$affinity, na.rm = TRUE),
        name = "Binding affinity\n(kcal/mol)"
      )
    } else {
      p <- p + paletteer::scale_fill_paletteer_c(
        palette, name = "Binding affinity\n(kcal/mol)"
      )
    }
  }

  if (isTRUE(label)) {
    p <- p + geom_text(
      aes(label = sprintf(paste0("%.", label_digits, "f"), .data[["affinity"]])),
      size = label_size,
      color = label_color,
      family = label_family,
      fontface = label_fontface
    )
  }

  return(p)
}


#' Retrieve ordering function for reordering factors

#' @param order Ordering option. Can be one of:
#' @return A function object or `NULL` if no ordering is requested.
#'
#' @importFrom stats median
#' @keywords internal

get_order <- function(order) {
  if (is.null(order) || identical(order, "none")) return(NULL)
  if (is.function(order)) return(order)
  if (is.character(order)) {
    m <- c("median" = stats::median, "mean" = base::mean,
           "max" = base::max, "min" = base::min)
    if (order %in% names(m)) return(m[[order]])
    stop("`order` only supports: NULL/'none'/'median'/'mean'/'max'/'min' or function.")
  }
  stop("`order` should be NULL/'none'/character/functions.")
}
