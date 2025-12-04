#' Prepare input data for enrich_circle_plot
#'
#' @param x enrichResult object or data.frame from clusterProfiler
#' with columns: ID, ONTOLOGY, BgRatio, p.adjust, geneID, Count, GeneRatio and RichFactor.
#' @param up_genes  Character vector of up-regulated genes (same ID space as \code{geneID}).
#' @param down_genes Character vector of down-regulated genes (same ID space as \code{geneID}).
#' @param top Numeric. Number of pathways to keep per category (by ascending \code{p.adjust}). Default 10.
#'
#' @importFrom dplyr mutate select group_by summarise arrange slice_head ungroup desc all_of
#' @importFrom tidyr drop_na separate_rows
#' @importFrom rlang .data
#' @importFrom methods is
#'
#' @export
getGores <- function(x,
                     up_genes = NULL,
                     down_genes = NULL,
                     top = 10) {

  if (inherits(x, "enrichResult") || (isS4(x) && methods::is(x, "enrichResult"))) {
    df <- x@result
  } else if (is.data.frame(x)) {
    df <- x
  } else {
    stop("x must be an enrichResult (S4) or a data.frame.")
  }

  req_cols <- c("ID","ONTOLOGY","BgRatio","p.adjust","geneID")
  miss <- setdiff(req_cols, names(df))
  if (length(miss) > 0) stop("Missing required columns: ", paste(miss, collapse = ", "))

  # select necessary columns
  df <- df %>%
    tidyr::drop_na(dplyr::all_of(c("ID", "ONTOLOGY", "BgRatio", "p.adjust", "geneID"))) %>%
    dplyr::mutate(
      termnumber = as.integer(sub("/.*", "", .data$BgRatio)),
      totalnumber = as.integer(sub(".*/", "", .data$BgRatio)),
      goterm = .data$ID,
      category = .data$ONTOLOGY,
      rich_factor = .data$RichFactor)

  # combine with up-down genes information
  long_df <- df %>%
    dplyr::select(dplyr::all_of(c("goterm", "category", "totalnumber", "termnumber", "p.adjust", "rich_factor", "geneID"))) %>%
    tidyr::separate_rows("geneID", sep = "/") %>%
    dplyr::mutate(
      is_up = if (is.null(up_genes)) NA_integer_ else as.integer(.data$geneID %in% up_genes),
      is_down = if (is.null(down_genes)) NA_integer_ else as.integer(.data$geneID %in% down_genes))

  input_full <- long_df %>%
    dplyr::group_by(.data$goterm, .data$category, .data$totalnumber,
                    .data$termnumber, .data$p.adjust, .data$rich_factor) %>%
    dplyr::summarise(
      up_regulated = if (all(is.na(.data$is_up))) NA_integer_ else sum(.data$is_up, na.rm = TRUE),
      down_regulated = if (all(is.na(.data$is_down))) NA_integer_ else sum(.data$is_down, na.rm = TRUE),
      .groups = "drop") %>%
    dplyr::arrange(.data$category, .data$p.adjust)

  # extract the top n
  input_top <- input_full %>%
    dplyr::group_by(.data$category) %>%
    dplyr::arrange(.data$p.adjust, .by_group = TRUE) %>%
    dplyr::slice_head(n = top) %>%
    dplyr::ungroup()

  # add additional columns
  res <- input_top %>%
    dplyr::mutate(ID = .data$goterm,
                  up_counts = .data$up_regulated,
                  down_counts = .data$down_regulated,
                  gene_num.min = 0,
                  "-log10(p.adjust)" = -log10(.data$p.adjust),
                  gene_num.max = log10(.data$totalnumber),
                  BgRatio1 = .data$termnumber,
                  BgRatio2 = .data$totalnumber,
                  RichFactor = .data$rich_factor) %>%
    dplyr::arrange(.data$category, desc(.data$RichFactor)) %>%
    as.data.frame()

  res$ID <- factor(res$ID, levels = res$ID)
  rownames(res) <- as.character(res$ID)

  return(res)
}


#' GO circle plot
#'
#' Visualize GO enrichment results as a circular plot highlighting p.adjust,
#' up/down-regulated genes, and RichFactor.
#'
#' @param x enrichResult object or data.frame from clusterProfiler.
#' @param top Number of top terms to show per category.
#' @param up_genes Vector of up-regulated genes (optional).
#' @param down_genes Vector of down-regulated genes (optional).
#' @param cat_col Named colors for GO categories (BP/CC/MF).
#' @param up_col Colors for up-regulated gene bars.
#' @param down_col Colors for down-regulated gene bars.
#' @param padjust_col Gradient colors for -log10(p.adjust).
#' @param bg.col Background color of outer ring.
#' @param max_width Maximum width of up/down bars.
#' @param fontsize Numeric. Font size for labels and legends.
#' @param fontface Character. Font style, e.g. "plain", "bold", "italic".
#' @param fontfamily Character. Font family, e.g. "sans", "serif", "mono".
#' @param ... Additional parameters for flexibility.
#'
#' @importFrom circlize circos.clear circos.par circos.genomicInitialize circos.track
#' @importFrom circlize circos.genomicTrack circos.genomicTrackPlotRegion circos.genomicRect
#' @importFrom circlize get.cell.meta.data circos.axis circos.text colorRamp2
#' @importFrom ComplexHeatmap Legend packLegend draw
#' @importFrom grid gpar unit
#' @importFrom graphics par
#' @importFrom grDevices colorRampPalette
#' @importFrom dplyr select mutate arrange group_by summarise ungroup all_of
#'
#' @return Invisibly returns TRUE after drawing the plot.
#' @export

gocircle_plot <- function(
    x,
    top,
    up_genes = NULL,
    down_genes = NULL,
    cat_col = c(BP = "#E64B35", CC = "#4DBBD5", MF = "#00A087"),
    up_col = "#AE2A8A",
    down_col = "#7B87BF",
    padjust_col = c("#FF906F", "#861D30"),
    bg.col = "gray95",
    max_width = 2.5,
    fontsize = 8,
    fontface = "bold",
    fontfamily = "sans",
    ...){

  x1 <- getGores(x, up_genes = up_genes, down_genes = down_genes, top = top)

  circlize::circos.clear()
  on.exit(circlize::circos.clear(), add = TRUE)
  circlize::circos.par(gap.degree = 1.2, start.degree = 90)

  plot_data <- x1[c("ID", "gene_num.min", "BgRatio2")]
  plot_data$BgRatio2 <- log10(plot_data$BgRatio2)
  circlize::circos.genomicInitialize(plot_data, plotType = NULL, major.by = 1)

  ## first round
  sector_colors <- unname(cat_col[as.character(x1$category)])
  circlize::circos.track(
    ylim = c(0, 10), track.height = 0.1, bg.border = NA, bg.col = sector_colors,
    panel.fun = function(x, y) {
      ylim <- circlize::get.cell.meta.data("ycenter")
      xlim <- circlize::get.cell.meta.data("xcenter")
      sector.name <- circlize::get.cell.meta.data("sector.index")

      circlize::circos.axis(
        h = "top", labels.cex = 0.4, major.at = c(0, 1, 2, 3, 4),
        labels = c("0", "10", "100", "1000", "10000"), minor.ticks = 0,
        labels.facing = "clockwise", labels.niceFacing = FALSE)

      circlize::circos.text(xlim, ylim, sector.name, cex = 0.5, col = "black", font = 2, niceFacing = FALSE)
    })

  ## second round
  plot_data <- x1[c("ID", "gene_num.min", "BgRatio1", "-log10(p.adjust)")]
  plot_data$BgRatio1 <- log10(plot_data$BgRatio1)
  p_max <- round(max(x1[["-log10(p.adjust)"]])) + 1
  p_col_fun <- circlize::colorRamp2(breaks = 0:p_max,
                                    col = colorRampPalette(padjust_col)(p_max + 1))

  circlize::circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.1, bg.border = NA, stack = TRUE,
    panel.fun = function(region, value, ...) {
      sector.name <- circlize::get.cell.meta.data("sector.index")
      pvalue_color <- x1[sector.name, "-log10(p.adjust)"]
      circlize::circos.genomicRect(region, value, col = p_col_fun(pvalue_color), border = NA, ...)
      ylim <- circlize::get.cell.meta.data("ycenter")
      xlim <- (region[1, 1] + region[1, 2]) / 2
      term_number <- x1[sector.name, "BgRatio1"]
      circlize::circos.text(xlim, ylim, term_number, cex = 0.6, col = "black", font = 2, niceFacing = FALSE)
    }
  )

  ## third round
  total_genes <- x1$up_counts + x1$down_counts
  x1$up_width <- ifelse(total_genes > 0, x1$up_counts / total_genes * max_width, 0)
  x1$down_width <- ifelse(total_genes > 0, x1$down_counts / total_genes * max_width, 0)

  plot_data_up <- data.frame(
    id = x1$ID,
    start = 0,
    end = x1$up_width,
    type = 1)

  plot_data_down <- data.frame(
    id = x1$ID,
    start = x1$up_width,
    end = x1$up_width + x1$down_width,
    type = 2)

  plot_data <- rbind(plot_data_up, plot_data_down)

  label_data <- data.frame(
    up_width = x1$up_width,
    total_width = x1$up_width + x1$down_width,
    up_counts = x1$up_counts,
    down_counts = x1$down_counts,
    row.names = x1$ID)

  circlize::circos.genomicTrackPlotRegion(
    plot_data, track.height = 0.1, bg.border = NA, stack = TRUE,
    panel.fun = function(region, value, ...) {
      rect_col_fun <- circlize::colorRamp2(breaks = c(1, 2), col = c(up_col, down_col))
      circlize::circos.genomicRect(region, value, col = rect_col_fun(value[[1]]), border = NA, ...)
      ylim <- circlize::get.cell.meta.data("cell.bottom.radius") - 0.35
      sector.index <- circlize::get.cell.meta.data("sector.index")

      if (label_data[sector.index, "up_counts"] > 0) {
        xlim <- label_data[sector.index, "up_width"] / 2
        circlize::circos.text(xlim, ylim, label_data[sector.index, "up_counts"], cex = 0.6, col = up_col, font = graphics::par("font"), niceFacing = FALSE)
      }

      if (label_data[sector.index, "down_counts"] > 0) {
        xlim <- (label_data[sector.index, "total_width"] + label_data[sector.index, "up_width"]) / 2
        circlize::circos.text(xlim, ylim, label_data[sector.index, "down_counts"], cex = 0.6, col = down_col, font = graphics::par("font"), niceFacing = FALSE)
      }
    }
  )

  ## fourth round
  plot_data <- x1 %>% dplyr::select(dplyr::all_of(c("ID", "gene_num.min", "gene_num.max", "RichFactor")))
  plot_data$RichFactor <- pmin(plot_data$RichFactor, 1.0)
  circlize::circos.genomicTrack(
    plot_data, ylim = c(0, 1.0), track.height = 0.25, bg.col = bg.col, bg.border = NA,
    panel.fun = function(region, value, ...) {
      sector.name <- circlize::get.cell.meta.data("sector.index")
      circlize::circos.genomicRect(
        region, value,
        col = cat_col[as.character(x1[sector.name, "category"])],
        border = NA, ytop.column = 1, ybottom = 0, ...)
    })

  ## legend
  text_gp <- grid::gpar(fontsize = fontsize, fontface = fontface, fontfamily = fontfamily)

  category_legend <- ComplexHeatmap::Legend(
    labels = names(cat_col), type = "points", pch = NA,
    background = unname(cat_col),
    labels_gp = text_gp)

  updown_legend <- ComplexHeatmap::Legend(
    labels = c("Up", "Down"), type = "points", pch = NA,
    background = c(up_col, down_col),
    labels_gp = text_gp)

  padj_legend <- ComplexHeatmap::Legend(
    col_fun = circlize::colorRamp2(round(seq(0, p_max, length.out = 6), 0), colorRampPalette(padjust_col)(6)),
    title = "-log10(p.adjust)",
    labels_gp = text_gp,
    title_gp = text_gp)

  lgd_list_vertical <- ComplexHeatmap::packLegend(category_legend, updown_legend)
  ComplexHeatmap::draw(x = grid::unit(0.43, "npc"), y = grid::unit(0.5, "npc"), lgd_list_vertical, just = c("left"))
  ComplexHeatmap::draw(padj_legend, x = grid::unit(0.51, "npc"), y = grid::unit(0.5, "npc"), just = c("left"))

  circlize::circos.clear()

  invisible(TRUE)
}
