#' Normalize a numeric vector to the range \eqn{[0, 1]}.
#'
#' A simple and robust 0–1 scaler. Handles NA, Inf, or constant values by
#' returning 0.5 for those entries.
#'
#' @param x Numeric vector.
#' @param na_rm Logical; whether to ignore non-finite values when computing
#'   the range. Defaults to \code{TRUE}.
#'
#' @return Numeric vector which has been scaled.
#'
#' @keywords internal
norm01 <- function(x, na_rm = TRUE) {
  x <- as.numeric(x)
  if (na_rm) {
    finite_idx <- is.finite(x)
    if (!any(finite_idx)) return(rep(0.5, length(x)))
    rng <- range(x[finite_idx], na.rm = TRUE)
  } else {
    rng <- range(x, na.rm = TRUE)
  }
  if (diff(rng) == 0) {
    return(rep(0.5, length(x)))
  } else {
    out <- (x - rng[1]) / diff(rng)
    out[!is.finite(out)] <- 0.5
    return(out)
  }
}


#' Extract normalized centrality profile for a node
#'
#' This function extracts selected metrics for a given node and normalizes
#' them to the range \eqn{[0, 1]} for use in plots such as radar charts.
#'
#' @param tab A data frame containing at least a \code{name} column and the specified metric columns.
#' @param node_name Character; node name to extract.
#' @param metrics Character vector; names of metric columns to use. Defaults to common centrality measures.
#'
#' @return A data frame with two columns: \code{metric} and \code{value},
#'   where \code{value} is the normalized score of the given node for each
#'   metric.
#'
#' @export
get_node_profile <- function(tab,
                             node_name,
                             metrics = c("degree",
                                         "betweenness",
                                         "closeness",
                                         "MCC",
                                         "MNC",
                                         "DMNC",
                                         "coreness",
                                         "EPC")) {
  stopifnot("name" %in% names(tab))

  metrics <- intersect(metrics, names(tab))
  if (length(metrics) == 0L) {
    stop("No requested metrics found in table.")
  }

  M_raw  <- tab[, metrics, drop = FALSE]
  M_norm <- as.data.frame(lapply(M_raw, norm01))

  idx <- which(tab$name == node_name)
  if (length(idx) == 0L) {
    stop("Node name not found in table.")
  }

  vals <- M_norm[idx, , drop = FALSE]

  df_profile <- data.frame(
    metric = metrics,
    value  = as.numeric(vals[1, ]),
    stringsAsFactors = FALSE
  )

  return(df_profile)
}


#' Plot a radar chart for node centrality profile
#'
#' This function draws a radar-style polygon plot for a single node,
#' showing its values on multiple centrality metrics.
#'
#' @param profile_df A data frame containing at least one column of metric
#'   names and one column of numeric values.
#' @param category_col Character string; column name in \code{profile_df}
#'   for metric labels. Defaults to \code{"metric"}.
#' @param value_col Character string; column name in \code{profile_df}
#'   for metric values (usually scaled to \eqn{[0, 1]}). Defaults to \code{"value"}.
#' @param fill_color Polygon fill color. "#A3BEDD", "#A7CBA9", and "#D59390" are recommended.
#' @param line_color Polygon border color.
#' @param title Character; plot title. Default is NULL.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot geom_polygon geom_segment geom_text coord_equal
#' @importFrom ggplot2 theme_void ggtitle theme element_text margin
#' @importFrom dplyr mutate bind_rows
#' @importFrom rlang sym .data
#'
#' @export
radar_plot <- function(profile_df,
                       category_col = "metric",
                       value_col = "value",
                       fill_color = "#B6B0D4",
                       line_color = "#6D65A6",
                       title = NULL) {
  data <- profile_df
  angle <- NULL

  if (!(category_col %in% colnames(data))){
    category_col <- colnames(data)[1]
  }

  if (!(value_col %in% colnames(data))){
    value_col <- colnames(data)[2]
  }

  n <- nrow(data)
  angles <- (0:(n - 1)) / n * 2 * pi

  data <- data %>%
    dplyr::mutate(
      angle = angles,
      x = sin(angle) * !!rlang::sym(value_col),
      y = cos(angle) * !!rlang::sym(value_col))

  df_polygon <- rbind(data, data[1, ])

  max_value <- ceiling(max(data[[value_col]]) * 10) / 10
  radii <- seq(0.2, max_value, length.out = 4)

  ring_polygons <- lapply(radii, function(r) {
    x <- sin(angles) * r
    y <- cos(angles) * r
    data.frame(
      x = c(x, x[1]),
      y = c(y, y[1]),
      r = r
    )
  }) %>%
    dplyr::bind_rows(.id = "ring")

  df_segments <- data %>%
    dplyr::mutate(
      xend = sin(angle) * max_value,
      yend = cos(angle) * max_value)

  axis_labels <- data.frame(
    x = radii,
    y = 0,
    label = sprintf("%.1f", radii))

  # radar plot
  p <- ggplot() +
    geom_polygon(
      data = ring_polygons,
      aes(.data$x, .data$y, group = .data$ring),
      color = "grey80", fill = NA
    ) +
    geom_segment(
      data = df_segments,
      aes(x = 0, y = 0, xend = .data$xend, yend = .data$yend),
      color = "grey85"
    ) +
    geom_polygon(
      data = df_polygon,
      aes(.data$x, .data$y),
      fill = fill_color, color = line_color,
      linewidth = 0.8, alpha = 0.7
    ) +
    geom_text(
      data = data,
      aes(
        x = sin(.data$angle) * (max_value * 1.25),
        y = cos(.data$angle) * (max_value * 1.25),
        label = !!rlang::sym(category_col)
      ),
      size = 3
    ) +
    geom_text(
      data = axis_labels,
      aes(x = .data$x, y = .data$y, label = .data$label),
      color = "grey40", size = 3, vjust = -0.5
    ) +
    coord_equal(clip = "off") +
    theme_void() +
    theme(
      plot.margin = margin(10, 30, 10, 10))

  if (!is.null(title)){
    p <- p + ggtitle(title) +
      theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
  }

  return(p)
}
