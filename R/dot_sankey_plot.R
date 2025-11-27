#' Extract enrichment table for Sankey-Dot plotting
#'
#' @param x enrichResult object or data.frame containing enrichment results.
#' @param n Integer. The number of pathways used for visualization.
#' @importFrom dplyr mutate select arrange slice desc
#' @importFrom tidyr drop_na
#' @importFrom tibble as_tibble
#' @importFrom rlang .data

#' @return a tibble with necessary columns for sankey-dot plot.
#' @keywords internal
getenrichres <- function(x, n = 10) {

  if (inherits(x, "enrichResult") || isS4(x)) {
    df <- x@result
  } else if (is.data.frame(x)) {
    df <- x
  } else {
    stop("x must be an enrichResult (S4) or a data.frame.")
  }
  df <- tidyr::drop_na(df)

  req_cols <- c("geneID", "Description", "Count", "GeneRatio")
  if (!all(req_cols %in% names(df))) {
    stop("Input must contain columns: geneID, Description, Count, GeneRatio.")
  }

  if (!"RichFactor" %in% names(df)) df$RichFactor <- NA_real_
  if (!"FoldEnrichment" %in% names(df)) df$FoldEnrichment <- NA_real_

  if ("p.adjust" %in% names(df)) {
    sig_col <- "p.adjust"
  }
  else {
    stop("No p.adjust column found in input.")
  }

  df <- df |>
    dplyr::mutate(
      Count = as.numeric(.data$Count),
      GeneRatio = ifelse(
        grepl("/", .data$GeneRatio),
        sapply(strsplit(.data$GeneRatio, "/"), function(x) {
          x <- as.numeric(x); x[1] / x[2]
        }), as.numeric(.data$GeneRatio)),
      sig = .data[[sig_col]]) |>
    dplyr::arrange(.data$sig) |>
    dplyr::slice(1:n) |>
    dplyr::select(id = .data$geneID, desc = .data$Description, count = .data$Count,
                  .data$sig, GeneRatio = .data$GeneRatio, RichFactor = .data$RichFactor, FoldEnrichment = .data$FoldEnrichment)

  df <- tibble::as_tibble(df)

  return(df)
}



#' Insert spacer nodes between strata for better visual separation in Sankey plots
#'
#' @param dat A data frame or tibble.
#' @return A tibble with additional "spacer" rows and updated factor levels
#'   in the `stratum` column.
#' @importFrom dplyr filter pull bind_rows
#' @importFrom rlang .data
#' @importFrom purrr map_df
#' @importFrom tibble as_tibble tibble
#' @keywords internal

insert_spacer_nodes <- function(dat) {
  d <- tibble::as_tibble(dat)
  new_lvl <- list()

  for (ax in unique(d$axis)) {
    old <- d %>%
      dplyr::filter(.data$axis == ax) %>%
      dplyr::pull(.data$stratum) %>%
      unique()

    nm   <- paste0("spacer_", ax, "_", seq_along(old)[-length(old)])
    ord  <- character(0)
    ord[2 * (seq_along(old) - 1) + 1] <- old
    if (length(nm)) ord[2 * seq_along(nm)] <- nm

    ref_y <- d %>%
      dplyr::filter(.data$axis == ax) %>%
      dplyr::pull(.data$y_pos)

    ref_y <- ref_y[1]

    if (length(nm)) {
      spacer_rows <- purrr::map_df(nm, function(spacer_name) {
        tibble::tibble(
          Freq     = 1,
          alluvium = spacer_name,
          axis     = ax,
          stratum  = spacer_name,
          y_pos    = ref_y
        )
      })
      d <- dplyr::bind_rows(d, spacer_rows)
    }
    new_lvl[[ax]] <- ord
  }

  d$stratum <- factor(d$stratum, levels = unlist(new_lvl))
  return(d)
}



#' set colors for sankey and dot plot
#' @param items Character vector of unique item names to assign colors to.
#' @param colors Character vector of base colors to use for mapping.
#' @param insert Optional character vector of additional colors to insert into the palette.
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom stats setNames
#' @keywords internal

make_colors <- function(items, colors, insert = NULL) {
  items <- unique(items)
  n   <- length(items)
  base <- if (!is.null(insert)) c(colors, insert) else colors
  if (n > length(base)) base <- colorRampPalette(base)(n)
  stats::setNames(base[seq_len(n)], items)
}



#' Dot–Sankey plot for enrichment results
#'
#' @param enrich_obj enrichResult object or data.frame containing enrichment results.
#' @param n Integer. The number of pathways to visualize (top `n` ranked by significance).
#' @param axis_order Character vector of length 2, specifying the order of axes in
#'   the Sankey diagram. Default is `c("Pathway", "Gene")`.
#' @param id_colors Character vector of base colors for gene nodes.
#' @param desc_colors Character vector of base colors for pathway nodes.
#' @param id_y_pos Numeric constant controlling the y-position of “Gene” strata in the Sankey diagram. Default is `1.1`.
#' @param desc_y_pos Numeric constant controlling the y-position of “Pathway” strata in the Sankey diagram. Default is `1.0`.
#' @param pathway_wrap Integer. The maximum line width for pathway label wrapping. Default is `50`.
#' @param sankey_text_size Numeric. Font size for text labels in the Sankey diagram. Default is `4`.
#' @param bubble_size_range Numeric vector of length 2. Range of point sizes in the dot plot. Default is `c(3, 8)`.
#' @param dot_palette Character. Name of the RColorBrewer palette used for color gradients in the dot plot. Default is `"RdBu"`.
#' @param dot_x_var Character. Variable used for the x-axis in the dot plot.
#' @param bubble_p_label Character. The column name for significance values used to color dotsd. Default is `p.adjust`.
#' @param sankey_width Numeric. Relative width of the Sankey panel in the combined plot. Default is `2`.
#' @param dot_width Numeric. Relative width of the dot plot panel in the combined plot. Default is `1`.
#' @param font_family Character. Font family for all text elements. Default is `"Arial"`.
#' @param font_face Character. Font face for all text. Default is `"bold"`.
#' @param sankey_lab Character. Label for the x-axis of the Sankey diagram. Default is `"Gene-Pathway"`.
#' @param seed Integer. Random seed for reproducibility of layout. Default is `2025`.
#' @param ... Additional arguments passed to internal helper functions.
#'
#' @import ggplot2
#' @importFrom ggalluvial geom_stratum geom_flow StatStratum to_lodes_form
#' @importFrom yulab.utils str_wrap
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @importFrom dplyr filter mutate count arrange select distinct left_join case_when desc
#' @importFrom tidyr separate_rows
#' @importFrom ggfun get_legend
#' @importFrom aplot plot_list
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales pretty_breaks
#' @return a ggplot object containing dot plot and sankey plot.
#' @export

ggdot_sankey <- function(
    enrich_obj,
    n = 10,
    axis_order = c("Pathway","Gene"),
    id_colors = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
                  "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F"),
    desc_colors = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F","#8491B4",
                     "#91D1C2", "#C8DC96", "#7E6148", "#B09C85","#6A5ACD", "#A0522D"),
    id_y_pos = 1.1,
    desc_y_pos = 1.0,
    pathway_wrap = 50,
    sankey_text_size = 4,
    bubble_size_range = c(3, 8),
    dot_palette = "RdBu",
    dot_x_var = c("GeneRatio", "RichFactor", "FoldEnrichment"),
    bubble_p_label = "p.adjust",
    sankey_width = 2,
    dot_width = 1,
    font_family = "sans",
    font_face = "bold",
    sankey_lab = "Gene-Pathway",
    seed = 2025,
    ...){

  dot_x_var <- match.arg(dot_x_var)

  bubble_x_label <- dot_x_var

  set.seed(seed)

  ## 1. get enrichment result
  et <- getenrichres(enrich_obj, n = n, ...)

  dfForLodes <- et %>%
    tidyr::separate_rows("id", convert = TRUE, sep = "/") %>%
    dplyr::count(Gene = .data$id , Pathway = .data$desc, order = .data$count, name = "Freq") %>%
    dplyr::arrange(dplyr::desc(.data$order))

  ## 2. prepare sankey data
  sankeyData <- ggalluvial::to_lodes_form(dfForLodes,
                                          key = "axis",
                                          axes = axis_order) %>%
    dplyr::mutate(y_pos = dplyr::case_when(
      .data$axis == "Gene" ~ id_y_pos,
      .data$axis == "Pathway" ~ desc_y_pos),
      alluvium = as.character(.data$alluvium),
      stratum = as.character(.data$stratum))

  sankeyData <- insert_spacer_nodes(sankeyData)

  ## 3. make colors
  gene_colors <- make_colors(items  = dfForLodes$Gene,colors = id_colors,insert = RColorBrewer::brewer.pal(12, "Paired"))
  pathway_colors <- make_colors(items  = dfForLodes$Pathway,colors = desc_colors,insert = RColorBrewer::brewer.pal(8, "Set1"))
  spacer_strata <- grep("spacer_", unique(sankeyData$stratum), value = TRUE)
  spacer_colors <- stats::setNames(rep("transparent", length(spacer_strata)), spacer_strata)
  nodeColors <- c(gene_colors, pathway_colors,spacer_colors)

  sankeyData <- sankeyData %>%
    mutate(axis = factor(.data$axis, levels = axis_order)) %>%
    mutate(node_color = nodeColors[as.character(.data$stratum)]) %>%   # 节点色
    mutate(flow_color = .data$node_color)

  ## 4. sankey plot
  base_theme <- theme(text = element_text(family = font_family, face = font_face),
                      plot.margin = unit(c(0, 0, 0, 0), "cm"))

  sankeyPlot <- ggplot(
    data = sankeyData,
    aes(x = .data$axis, stratum = .data$stratum, alluvium = .data$alluvium, y = .data$y_pos)
  ) +
    geom_stratum(aes(fill = .data$node_color), color = NA, width = 0.05) +
    geom_flow(aes(fill = .data$flow_color), alpha = 0.3, width = 0.05,
              knot.pos = 0.3, color = "transparent") +
    geom_text(
      stat = ggalluvial::StatStratum,
      data = function(x) dplyr::filter(x, .data$axis == "Gene"),
      aes(label = ifelse(
        grepl("^spacer_", as.character(after_stat(.data$stratum))),
        "", as.character(after_stat(.data$stratum))
      )),
      hjust = 0, nudge_x = 0.03,
      size = sankey_text_size, family = font_family, fontface = font_face) +
    geom_text(
      stat = ggalluvial::StatStratum,
      data = function(x) dplyr::filter(x, .data$axis == "Pathway"),
      aes(label = yulab.utils::str_wrap(
        ifelse(
          grepl("^spacer_", as.character(after_stat(.data$stratum))),
          "", as.character(after_stat(.data$stratum))), width = pathway_wrap)),
      hjust = 0, nudge_x = 0.03, size = sankey_text_size, family = font_family, fontface = font_face) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_identity() +
    guides(fill = "none") +
    theme_void() +
    base_theme +
    labs(x = sankey_lab) +
    theme(axis.title.x = element_text(margin = margin(t = 6), size = 16))

  ## 5. dot plot
  # extract plot data
  sankeyPlotData <- ggplot2::ggplot_build(sankeyPlot)

  leftNodes <- sankeyPlotData$data[[1]] %>%
    dplyr::filter(.data$x == min(.data$x)) %>%
    dplyr::mutate(node_name = as.character(.data$stratum),
           node_ymin = .data$ymin,
           node_ymax = .data$ymax,
           node_center_y = (.data$ymin + .data$ymax) / 2) %>%
    dplyr::filter(!grepl("spacer_", .data$node_name)) %>%
    dplyr::select(.data$node_name, .data$node_center_y, .data$ymin, .data$ymax)

  dotData <- et %>%
    distinct(.data$desc, .keep_all = TRUE) %>%
    left_join(leftNodes, by = c("desc" = "node_name"))

  dotPlot <- ggplot(
    data = dotData,
    aes(x = .data[[dot_x_var]], y = .data$node_center_y, color = -log10(.data$sig))
  ) +
    geom_point(aes(size = .data$count), stroke = 0.5) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0.003)) +
    scale_color_distiller(
      palette = dot_palette,
      direction = -1,
      name = paste0("-log10(", bubble_p_label, ")")
    ) +
    scale_size_continuous(
      range = bubble_size_range,
      name = "Count",
      breaks = scales::pretty_breaks(n = 4)
    ) +
    guides(
      color = guide_colorbar(order = 1, barwidth = 0.8, barheight = 4),
      size = guide_legend(
        order = 2,
        keywidth = 0.8, keyheight = 0.8,
        override.aes = list(color = "black", alpha = 1))) +
    labs(
      size = "Count",
      color = paste0("-log10(", bubble_p_label, ")"),
      x = bubble_x_label,
      y = NULL
    ) +
    theme_bw() + base_theme +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(margin = margin(t = 4), size = 12),
      axis.title.x = element_text(margin = margin(t = 6), size = 16))

  ## combine dot plot and sankey plot
  yRange <- sankeyPlotData$layout$panel_params[[1]]$y.range

  sankeyPlot <- sankeyPlot + coord_cartesian(clip = "off", ylim = yRange)+
    theme(plot.margin = margin(0, 2.5, 0, 0.05, "cm"))

  dotPlot <- dotPlot +
    annotate("rect",
             xmin = -Inf, xmax = Inf,
             ymin = -Inf, ymax = max(leftNodes$ymax),
             fill = NA, color = "black", linewidth = 0.3) +
    coord_cartesian(ylim = yRange)


  leg <- ggfun::get_legend(dotPlot)
  combinedPlot <- plot_list(dotPlot + theme(legend.position = "none"), sankeyPlot,
                            ncol = 2, widths = c(dot_width, sankey_width))

  combinedPlot <- aplot::plot_list(leg, combinedPlot, ncol = 2, widths = c(1, 5))

  return(combinedPlot)
}
