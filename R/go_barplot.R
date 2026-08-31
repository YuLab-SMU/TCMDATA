#' Bar Plot for GO Enrichment Results
#'
#' Create a bar plot for GO enrichment analysis results, with bars colored
#' by ontology category (BP, CC, MF). Based on barplot method from enrichplot.
#'
#' @param enrich_obj An enrichResult object from \code{clusterProfiler::enrichGO()}.
#' @param x Character. The variable to be plotted on the x-axis. Default is "Count".
#' @param order Logical. Whether to sort `x`. Default is TRUE.
#' @param top_n Integer. Number of top terms to display per ontology category. Default is 10.
#' @param colors Named character vector for BP, CC, MF colors. Default is \code{c(BP = "#E64B35", CC = "#4DBBD5", MF = "#00A087")}.
#' @param label.size Numeric. Font size for count labels on bars. Default is 3.2.
#' @param show_count Logical. Whether to show count labels on bars. Default is TRUE.
#' @param label.bold Logical. Whether count labels should be bold. Default is FALSE.
#' @param base_size Numeric. Base font size. Default is 7.
#' @param base_family Character. Font family used throughout the plot.
#'   Default is `"Arial"`.
#' @param x.angle Numeric. Rotation angle of x-axis text (pathway descriptions). Default is 60.
#' @param legend.position Character. Position of the legend. Default is "top".
#' @param plot_title Character or NULL. Title of the plot. Default is NULL.
#' @param ... Additional parameters if neccessary.
#'
#' @return A ggplot2 object.
#'
#' @import ggplot2
#' @import enrichplot
#' @importFrom graphics barplot
#' @importFrom dplyr group_by slice_max ungroup mutate arrange desc
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## demo data
#' herbs <- c("灵芝")
#' lz <- search_herb(herb = herbs, type = "Herb_cn_name")
#' set.seed(2025)
#' g <- sample(lz$target, 200)
#'
#' ## GO enrichment analysis
#' library(clusterProfiler)
#' x <- enrichGO(g, ont="all", OrgDb='org.Hs.eg.db', keyType="SYMBOL")
#' p1 <- go_barplot(x)
#' print(p1)
#'
#' }
go_barplot <- function(enrich_obj,
                       x = "Count",
                       order = TRUE,
                       top_n = 10,
                       colors = NULL,
                       label.size = 3.2,
                       show_count = TRUE,
                       label.bold = FALSE,
                       base_size = 7,
                       base_family = "Arial",
                       x.angle = 60,
                       legend.position = "top",
                       plot_title = NULL,
                       ...) {

  base_family <- .resolve_tcm_font_family(base_family)
  # Check input
  if (!inherits(enrich_obj, "enrichResult")) {
    stop("enrich_obj must be an enrichResult object", call. = FALSE)
  }

  # Default colors
  if (is.null(colors)) {
    colors <- c("BP" = "#E64B35", "CC" = "#4DBBD5", "MF" = "#00A087")
  }

  label_mapping <- c(
    "BP" = "Biological Process",
    "CC" = "Cellular Component",
    "MF" = "Molecular Function"
  )

  # Use barplot method from enrichplot
  p <- graphics::barplot(enrich_obj,
                         x = x,
                         showCategory = top_n,
                         split = "ONTOLOGY",
                         order = order,
                         ...)

  # Get data from barplot object and reorder
  df <- p$data
  if (!"ONTOLOGY" %in% colnames(df)) {
    stop("No ONTOLOGY column found. Please use enrichGO(ont='ALL').")
  }

  xcol <- x

  # Reorder
  df <- df |>
    mutate(ONTOLOGY = factor(.data$ONTOLOGY, levels = c("BP", "CC", "MF"))) |>
    arrange(.data$ONTOLOGY, desc(.data[[xcol]]))

  desired_order <- df$Description

  # update
  p <- p + df

  p <- p +
    aes(fill = .data$ONTOLOGY) +
    scale_fill_manual(values = colors, labels = label_mapping, name = NULL) +
    coord_flip() +
    scale_y_discrete(limits = desired_order)

  # Create color mapping matching the plot order
  label_colors <- sapply(desired_order, function(desc) {
    ont <- as.character(df$ONTOLOGY[df$Description == desc][1])
    colors[ont]
  })

  # Add count labels on top of bars
  if (show_count) {
    p <- p + geom_text(
      aes(label = .data[[xcol]]),
      vjust = -0.3,
      size = label.size,
      family = base_family,
      fontface = if (label.bold) "bold" else "plain"
    )
  }

  # update layout
  p <- p +
    .theme_tcm_pub(
      base_size = base_size,
      base_family = base_family,
      grid = "none"
    ) +
    theme(
      axis.text.x = element_text(angle = x.angle, hjust = 1, vjust = 1,
                                 color = label_colors, face  = if (label.bold) "bold" else "plain"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.y = element_line(color = "black"),
      axis.line.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = legend.position,
      plot.margin = margin(t = 5, r = 20, b = 5, l = 60, unit = "pt")
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1)))

  if (!is.null(plot_title)){
    p <- p + ggtitle(plot_title)
  }

  return(p)
}
