#' TCM Sankey Plot for herb-molecule-target
#'
#' @param data A data frame containing at least three columns: \code{herb}, \code{molecule}, and \code{target}.
#' @param axis_order Character vector specifying the order of axes in the Sankey diagram. Default is \code{c("herb", "molecule", "target")}.
#' @param herb_cols Character vector defining the base color palette for the herb layer.
#' @param mol_cols Character vector defining the base color palette for the molecule.
#' @param target_cols Character vector defining the base color palette for the target.
#' @param plot_font Character string specifying the font family used for text. Default is \code{"sans"}.
#' @param font_face Character string specifying the font face. Default is \code{"plain"}.
#' @param target_fontface Character string specifying the font face for target labels (rightmost axis). Default is \code{"italic"}.
#' @param font_size Numeric value controlling the size of node labels. Default is \code{3.5}.
#' @param width Numeric value controlling the width of both nodes and flows. Default is \code{0.05}.
#' @param alpha Numeric value controlling the transparency of the flows. Default is \code{0.3}.
#' @param knot.pos Numeric value (between 0 and 1) determining the curvature position of flow lines. Default is \code{0.3}.
#' @param ... Additional arguments passed to `tcm_sankey` when using the deprecated alias `TCM_sankey`.
#' @import ggplot2
#' @importFrom dplyr count mutate group_by ungroup arrange select filter lead case_when all_of desc
#' @importFrom yulab.utils str_wrap
#' @importFrom stats setNames
#' @importFrom rlang .data
#'
#' @return A \code{ggplot} object representing the herb–compound–target Sankey diagram.
#' @export
tcm_sankey <- function(
    data,
    axis_order = c("herb", "molecule", "target"),
    herb_cols = c("#4E79A7", "#59A14F", "#B07AA1", "#9C755F",
                  "#76B7B2", "#BAB0AC", "#A0CBE8", "#8CD17D"),
    mol_cols = c("#D8A24A", "#86BCB6", "#A0CBE8", "#FFBE7D",
                 "#B7D59A", "#E6A0A0", "#C5B0D5", "#B85C5C", "#4E79A7"),
    target_cols = c("#C49C94", "#D7A1B8", "#C9C982", "#9CCB9A",
                    "#D7B7D8", "#B7D59A", "#E6D48A", "#B85C5C",
                    "#D8A24A", "#E15759"),
    plot_font = "Arial",
    font_face = "plain",
    target_fontface = "italic",
    font_size = 2.2,
    width = 0.05,
    alpha = 0.3,
    knot.pos = 0.3)
{
  if (!requireNamespace("ggalluvial", quietly = TRUE)) {
    stop("Package 'ggalluvial' is required for tcm_sankey(). Please install it.")
  }

  if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
    stop("Package 'RColorBrewer' is required for tcm_sankey(). Please install it with: install.packages('RColorBrewer')")
  }

  set.seed(2025)

  herb_y_high <- mol_y_high <- target_y_high <- 1

  ## prepare input data
  dfForLodes <- data %>%
    count(.data$herb, .data$molecule, .data$target, name = "Freq") %>%
    group_by(.data$molecule) %>%
    mutate(MolTotal = sum(.data$Freq)) %>%
    ungroup() %>%
    arrange(dplyr::desc(.data$MolTotal)) %>%
    select(!dplyr::all_of("MolTotal"))

  ## convert data type
  sankeyData <- ggalluvial::to_lodes_form(dfForLodes, key = "axis", axes = axis_order) %>%
    mutate(y_pos = dplyr::case_when(
      .data$axis == "herb" ~ herb_y_high,
      .data$axis == "molecule" ~ mol_y_high,
      .data$axis == "target" ~ target_y_high),
      alluvium = as.character(.data$alluvium),
      stratum = as.character(.data$stratum))

  sankeyData <- insert_spacer_nodes(sankeyData)

  ## make colors
  herb_colors <- make_colors(items = dfForLodes$herb, colors = herb_cols, insert = RColorBrewer::brewer.pal(12, "Paired"))
  mol_colors <- make_colors(items = dfForLodes$molecule, colors = mol_cols, insert = RColorBrewer::brewer.pal(8, "Set1"))
  target_cols <- make_colors(items = dfForLodes$target,colors = target_cols, insert = RColorBrewer::brewer.pal(8, "Set2"))
  spacer_strata <- grep("spacer_", unique(sankeyData$stratum), value = TRUE)
  spacer_colors <- stats::setNames(rep("transparent", length(spacer_strata)), spacer_strata)
  nodeColors <- c(herb_colors, mol_colors, target_cols, spacer_colors)

  sankeyData <- sankeyData %>%
    dplyr::mutate(axis = factor(.data$axis, levels = axis_order)) %>%
    dplyr::mutate(node_color = nodeColors[as.character(.data$stratum)]) %>%
    dplyr::group_by(.data$alluvium) %>%
    dplyr::mutate(
      to_node_name = dplyr::lead(as.character(.data$stratum), order_by = .data$axis),
      to_node_name = ifelse(is.na(.data$to_node_name), as.character(.data$stratum), .data$to_node_name)) %>%
    ungroup() %>%
    mutate(flow_color = nodeColors[.data$to_node_name])

  ## plot sankey
  sankeyPlot <- ggplot(data = sankeyData,
                       aes(x = .data$axis, stratum = .data$stratum, alluvium = .data$alluvium,
                           y = .data$y_pos, label = .data$stratum)) +
    ggalluvial::geom_stratum(aes(fill = .data$node_color), color = NA, width = width) +
    ggalluvial::geom_flow(aes(fill = .data$flow_color), alpha = alpha, width = width,
              knot.pos = knot.pos, color = "transparent") +

    # herb
    geom_text(stat = ggalluvial::StatStratum,
              data = function(x) dplyr::filter(x, .data$axis == "herb"),
              aes(label = ifelse(grepl("spacer_", as.character(after_stat(.data$stratum))), "",
                                 as.character(after_stat(.data$stratum)))),
              hjust = 1, nudge_x = -0.03,
              size = font_size, family = plot_font, fontface = font_face) +

    # molecule
    geom_text(stat = ggalluvial::StatStratum,
              data = function(x) dplyr::filter(x, .data$axis == "molecule"),
              aes(label = ifelse(grepl("spacer_", as.character(after_stat(.data$stratum))), "",
                                 as.character(after_stat(.data$stratum)))),
              hjust = 0.5, nudge_x = 0,
              size = font_size, family = plot_font, fontface = font_face) +

    # target
    geom_text(stat = ggalluvial::StatStratum,
              data = function(x) dplyr::filter(x, .data$axis == "target"),
              aes(label = yulab.utils::str_wrap(
                ifelse(grepl("spacer_", as.character(after_stat(.data$stratum))), "",
                       as.character(after_stat(.data$stratum))),50)),
              hjust = 0, nudge_x = 0.03,
              size = font_size, family = plot_font, fontface = target_fontface) +
    scale_x_discrete(limits = axis_order, expand = c(0, 0)) +
    scale_fill_identity() +
    guides(fill = "none") +
    .theme_tcm_void(base_size = 7, base_family = plot_font) +
    theme(plot.margin = margin(t = 4, r = 28, b = 4, l = 28, unit = "pt")) +
    coord_cartesian(clip = "off")

  return(sankeyPlot)

}

#' @rdname tcm_sankey
#' @export
TCM_sankey <- function(...) {
  warning("TCM_sankey is deprecated. Please use tcm_sankey instead.")
  tcm_sankey(...)
}









