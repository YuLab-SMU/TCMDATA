#' TCM Sankey Plot for herb-molecule-target
#'
#'
#' @param data A data frame containing at least three columns: \code{herb}, \code{molecule}, and \code{target}.
#' @param axis_order Character vector specifying the order of axes in the Sankey diagram. Default is \code{c("herb", "molecule", "target")}.
#' @param herb_cols Character vector defining the base color palette for the herb layer.
#' @param mol_cols Character vector defining the base color palette for the molecule.
#' @param target_cols Character vector defining the base color palette for the target.
#' @param plot_font Character string specifying the font family used for text. Default is \code{"sans"}.
#' @param font_face Character string specifying the font face. Default is \code{"bold"}.
#' @param font_size Numeric value controlling the size of node labels. Default is \code{3.5}.
#' @param width Numeric value controlling the width of both nodes and flows. Default is \code{0.05}.
#' @param alpha Numeric value controlling the transparency of the flows. Default is \code{0.3}.
#' @param knot.pos Numeric value (between 0 and 1) determining the curvature position of flow lines. Default is \code{0.3}.
#'
#' @import ggplot2
#' @importFrom ggalluvial geom_stratum geom_flow to_lodes_form
#' @importFrom dplyr count mutate group_by ungroup arrange select filter lead case_when all_of desc
#' @importFrom RColorBrewer brewer.pal
#' @importFrom yulab.utils str_wrap
#' @importFrom stats setNames
#' @importFrom rlang .data
#'
#' @return A \code{ggplot} object representing the herb–compound–target Sankey diagram.
#' @export

TCM_sankey <- function(
    data,
    axis_order = c("herb", "molecule", "target"),
    herb_cols = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                   "#8c564b", "#e377c2", "#7f7f7f"),
    mol_cols = c("#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a",
                   "#ff9896", "#c5b0d5", "#E41A1C", "#377EB8"),
    target_cols = c("#c49c94", "#f7b6d2", "#dbdb8d", "#c7e9c0", "#f4cae4",
                   "#e6f598", "#ffeda0", "#BB5234", "#BB7813", "#FF6158"),
    plot_font = "sans",
    font_face = "bold",
    font_size = 3.6,
    width = 0.05,
    alpha = 0.3,
    knot.pos = 0.3)
{
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
  herb_colors <- make_colors(items = dfForLodes$herb, colors = herb_cols, insert = brewer.pal(12, "Paired"))
  mol_colors <- make_colors(items = dfForLodes$molecule, colors = mol_cols, insert = brewer.pal(8, "Set1"))
  target_cols <- make_colors(items = dfForLodes$target,colors = target_cols, insert = brewer.pal(8, "Set2"))
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
    geom_text(stat = "stratum",
              data = function(x) dplyr::filter(x, .data$axis == "herb"),
              aes(label = ifelse(grepl("spacer_", as.character(after_stat(.data$stratum))), "",
                                 as.character(after_stat(.data$stratum)))),
              hjust = 1, nudge_x = -0.03,
              size = font_size, family = plot_font, fontface = font_face) +

    # molecule
    geom_text(stat = "stratum",
              data = function(x) dplyr::filter(x, .data$axis == "molecule"),
              aes(label = ifelse(grepl("spacer_", as.character(after_stat(.data$stratum))), "",
                                 as.character(after_stat(.data$stratum)))),
              hjust = 0.5, nudge_x = 0,
              size = font_size, family = plot_font, fontface = font_face) +

    # target
    geom_text(stat = "stratum",
              data = function(x) dplyr::filter(x, .data$axis == "target"),
              aes(label = yulab.utils::str_wrap(
                ifelse(grepl("spacer_", as.character(after_stat(.data$stratum))), "",
                       as.character(after_stat(.data$stratum))),50)),
              hjust = 0, nudge_x = 0.03,
              size = font_size, family = plot_font, fontface = font_face) +
    scale_x_discrete(limits = c("herb", "molecule", "target"), expand = c(0, 0)) +
    scale_fill_identity() +
    guides(fill = "none") +
    theme_void() +
    theme(plot.margin = margin(t = 0.8, r = 4, b = 1, l = 4, unit = "cm")) +
    coord_cartesian(clip = "off")

  return(sankeyPlot)

}










