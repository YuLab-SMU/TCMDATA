#' Lollipop Plot for Enrichment Results
#'
#' @param enrich_obj An enrichment result object from clusterProfiler.
#' @param x Character. The variable used for the x-axis. Default is `"RichFactor"`.
#' @param top_n Integer. Number of top enriched terms to display. Default is 10.
#' @param orderBy Character. Variable used to order the y-axis terms. Default is `"x"`.
#' @param text.col Character. Colors for text. Default is `black`.
#' @param text.size Numeric. Font size for axis text and title. Default is 8.
#' @param text.width Numeric. Font width for axis text and title. Default is 35.
#' @param palette Character. Color palette name from `RColorBrewer` to use for dot color. Default is `"RdBu"`.
#' @param line.col Character. Color of the segment lines. Default is `"grey60"`.
#' @param line.type Character. Line type for segments. Default is `"solid"`.
#' @param line.size Numeric. Line width for segments. Default is 0.9.
#' @param plot_title Character. Optional plot title. Default is `NULL`.
#' @param show_count Logical. Whether to display the count value as a text label next to each dot. Default is TRUE
#' @param ... Additional arguments passed to internal helper functions.

#' @return A `ggplot` object showing a lollipop-style enrichment plot.
#'
#' @import ggplot2
#' @importFrom ggrepel geom_text_repel
#' @importFrom RColorBrewer brewer.pal.info brewer.pal
#' @importFrom enrichplot dotplot
#' @importFrom rlang .data
#'
#' @export

gglollipop <- function(enrich_obj,
                               x = "RichFactor",
                               top_n = 10,
                               orderBy = NULL,
                               text.col = "black",
                               text.size = 8,
                               text.width = 35,
                               palette = "RdBu",
                               line.col = "grey60",
                               line.type = "solid",
                               line.size = 0.9,
                               plot_title = NULL,
                               show_count = TRUE,
                               ...) {

  if (is.null(orderBy)) {
    orderBy <- "x"
  }

  p <- enrichplot::dotplot(enrich_obj,
                           showCategory = top_n,
                           x = x,
                           orderBy = orderBy,
                           ...)

  df_plot <- p$data
  n_color <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]

  ## remove former color_scale
  p$scales$scales <- Filter(function(s) {
    !"ScaleContinuous" %in% class(s) || s$aesthetics != "fill"
  }, p$scales$scales)

  p <- p +
    geom_segment(data = df_plot,
                 aes(x = 0, xend = .data[[x]],
                     y = .data[["Description"]], yend = .data[["Description"]]),
                 color = line.col,
                 linetype = line.type,
                 linewidth = line.size,
                 inherit.aes = FALSE) +
    scale_fill_gradientn(
      colours = RColorBrewer::brewer.pal(n_color, palette),
      trans = "log10",
      guide = guide_colorbar(reverse = TRUE)
    )

  p$layers <- rev(p$layers)

  if (show_count){
    p <- p + geom_text_repel(data = df_plot,
              aes(x = .data[[x]],
                  y = .data[["Description"]],
                  label = .data[["Count"]]),
              size = text.size / 2,
              color = text.col,
              hjust = 0,
              direction = "x",
              nudge_x = 0.025 * max(df_plot[[x]], na.rm = TRUE),
              segment.color = NA,
              inherit.aes = FALSE)
  }

  p <- p +
    theme_bw() +
    theme(
      axis.text.x = element_text(size = text.size),
      axis.text.y = element_text(size = text.size),
      axis.title = element_text(size = text.size + 1),
      plot.title = element_text(size = text.size + 2, face = "bold", hjust = 0.5)
    ) +
    ggtitle(plot_title)

  return(p)
}
