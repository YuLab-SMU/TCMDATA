## Single-gene diagnostic analysis
## Functions for evaluating individual gene diagnostic ability after ML-based consensus feature selection.

## Internal formatting helpers
#' Format a p-value as a compact string
#' @keywords internal
#' @noRd
.fmt_p <- function(pval) {
  ifelse(pval < 0.001, sprintf("%.2e", pval),
         sprintf("%.4f", pval))
}

#' Convert a p-value to significance stars
#' @keywords internal
#' @noRd
.signif_stars <- function(pval) {
  ifelse(pval < 0.0001, "****",
  ifelse(pval < 0.001, "***",
  ifelse(pval < 0.01, "**",
  ifelse(pval < 0.05, "*", "ns"))))
}

## Input resolution helper
## make input format consistent
#' @keywords internal
#' @noRd
.resolve_expr_group <- function(genes,
                                ml_data = NULL,
                                expr_mat = NULL,
                                group = NULL,
                                use = c("auto", "train", "test")) {
  use <- match.arg(use)

  if (!is.null(ml_data)) {
    stopifnot(inherits(ml_data, "tcm_ml_data"))
    has_test <- !is.null(ml_data$test_x) && !is.null(ml_data$test_y)

    if (use == "auto")  use <- if (has_test) "test" else "train"
    if (use == "test" && !has_test)
      stop("ml_data has no test set (Mode A). Use use = 'train'.", call. = FALSE)

    if (use == "test") {
      expr_df <- ml_data$test_x
      grp <- ml_data$test_y
    } else {
      expr_df <- ml_data$train_x
      grp <- ml_data$train_y
    }
    lvls <- ml_data$levels

    ## Subset to requested genes
    safe_genes <- make.names(genes, unique = TRUE)
    avail <- intersect(safe_genes, colnames(expr_df))
    if (length(avail) == 0L)
      stop("None of the supplied genes found in the expression data.", call. = FALSE)
    name_map <- stats::setNames(genes[match(avail, safe_genes)], avail)
    expr_df <- expr_df[, avail, drop = FALSE]

  } else {
    ## from raw expr_mat + group
    if (is.null(expr_mat) || is.null(group))
      stop("Provide either `ml_data` or both `expr_mat` and `group`.",
           call. = FALSE)
    expr_mat <- as.matrix(expr_mat)
    group <- factor(group)
    stopifnot(ncol(expr_mat) == length(group), nlevels(group) == 2L)

    avail <- intersect(genes, rownames(expr_mat))
    if (length(avail) == 0L)
      stop("None of the supplied genes found in rownames(expr_mat).",
           call. = FALSE)

    expr_df <- as.data.frame(t(expr_mat[avail, , drop = FALSE]))
    colnames(expr_df) <- make.names(avail, unique = TRUE)
    name_map <- stats::setNames(avail, colnames(expr_df))
    grp  <- group
    lvls <- levels(group)
  }

  res <- list(expr = expr_df, group = grp, levels = lvls, name_map = name_map)
  return(res)
}


## get gene auc values

#' Compute single-gene AUC for diagnostic evaluation
#'
#' For each gene the function treats its expression as a univariate predictor,
#' computes a ROC curve via \pkg{pROC}, and returns a summary table with AUC,
#' 95 \% DeLong confidence interval, optimal cutoff (Youden's J), sensitivity,
#' specificity and a one-sided test against AUC = 0.5.
#'
#' @param genes Character vector of gene symbols to evaluate.
#' @param ml_data A \code{tcm_ml_data} object from \code{\link{prepare_ml_data}}
#'   (preferred input).  When provided, \code{expr_mat} / \code{group} are
#'   ignored.
#' @param expr_mat Numeric matrix (genes \eqn{\times} samples) with
#'   \code{rownames} = gene symbols.  Used only when \code{ml_data} is
#'   \code{NULL}.
#' @param group Factor or character vector of sample labels (length =
#'   \code{ncol(expr_mat)}, two levels).
#' @param use Which data split to evaluate when using \code{ml_data}:
#'   \code{"auto"} (default; test set if available, else train),
#'   \code{"train"}, or \code{"test"}.
#' @param ci Logical; compute 95 \% DeLong confidence interval?
#'   Default \code{TRUE}.
#' @param p_adjust_method Method for multiple-testing correction passed to
#'   \code{\link[stats]{p.adjust}}.  Default \code{"BH"} (Benjamini--Hochberg
#'   FDR).  Other common choices: \code{"bonferroni"}, \code{"holm"},
#'   \code{"none"}.  When only a single gene is supplied no correction is
#'   needed and \code{p_adj} equals \code{p_value} regardless of method.
#'
#' @return A \code{data.frame} sorted by descending AUC with columns:
#'   \code{gene}, \code{auc}, \code{ci_lower}, \code{ci_upper},
#'   \code{direction}, \code{optimal_cutoff}, \code{sensitivity},
#'   \code{specificity}, \code{youden_J}, \code{p_value}, \code{p_adj}
#'   (BH-corrected).
#'
#' @examples
#' \dontrun{
#'   consensus <- get_ml_consensus(ml_list, min_methods = 2)
#'   auc_tbl <- get_gene_auc(consensus, ml_data = ml_data)
#'   auc_tbl
#' }
#' @importFrom stats p.adjust qnorm pnorm
#' @export
get_gene_auc <- function(genes,
                         ml_data = NULL,
                         expr_mat = NULL,
                         group = NULL,
                         use = c("auto", "train", "test"),
                         ci = TRUE,
                         p_adjust_method = "BH") {

  p_adjust_method <- match.arg(p_adjust_method,
                               choices = stats::p.adjust.methods)

  if (!requireNamespace("pROC", quietly = TRUE))
    stop("Package 'pROC' is required. Install with install.packages('pROC').",
         call. = FALSE)

  dat <- .resolve_expr_group(genes, ml_data = ml_data, expr_mat = expr_mat,
                             group = group, use = use)
  expr_df <- dat$expr
  grp <- dat$group
  lvls <- dat$levels
  name_map <- dat$name_map

  rows <- lapply(colnames(expr_df), function(col) {
    x <- expr_df[[col]]
    roc_obj <- pROC::roc(response = grp, predictor = x,
                         levels = lvls, quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))

    ci_lo <- ci_hi <- NA_real_
    if (isTRUE(ci)) {
      ci_obj <- pROC::ci.auc(roc_obj, method = "delong")
      ci_lo  <- ci_obj[1L]
      ci_hi  <- ci_obj[3L]
    }

    ## Optimal cutoff by Youden's J
    best <- pROC::coords(roc_obj, x = "best", best.method = "youden",
                         ret = c("threshold", "sensitivity", "specificity"),
                         transpose = FALSE)
    ## coords may return multiple rows — pick first
    if (nrow(best) > 1L) best <- best[1L, , drop = FALSE]

    ## Test: AUC != 0.5 (two-sided, via DeLong SE from CI)
    p_val <- tryCatch({
      ci_obj_p <- if (isTRUE(ci)) ci_obj else pROC::ci.auc(roc_obj, method = "delong")
      se_auc <- (ci_obj_p[3L] - ci_obj_p[1L]) / (2 * stats::qnorm(0.975))
      z <- (auc_val - 0.5) / se_auc
      2 * stats::pnorm(-abs(z))
    }, error = function(e) NA_real_)

    data.frame(
      gene = name_map[[col]],
      auc = auc_val,
      ci_lower = ci_lo,
      ci_upper = ci_hi,
      direction = roc_obj$direction,
      optimal_cutoff = best$threshold,
      sensitivity = best$sensitivity,
      specificity = best$specificity,
      youden_J = best$sensitivity + best$specificity - 1,
      p_value = p_val,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  out <- out[order(-out$auc), ]
  out$p_adj <- stats::p.adjust(out$p_value, method = p_adjust_method)
  rownames(out) <- NULL
  return(out)
}


## plot roc curves 

#' ROC curves for individual diagnostic genes
#'
#' Overlays one ROC curve per gene on a single panel (\code{combine = TRUE},
#' default), or draws one ROC per facet (\code{combine = FALSE}).
#' This is the natural follow-up after \code{\link{get_ml_consensus}}: the
#' consensus genes are fed back to evaluate their individual diagnostic power.
#'
#' @inheritParams get_gene_auc
#' @param combine Logical; if \code{TRUE} (default) all curves are overlaid on one panel.  If \code{FALSE} each gene gets its own facet.
#' @param palette Character vector of colours (one per gene), or \code{NULL} for the built-in 12-colour publication palette.
#' @param show_ci Logical; Whether to add a translucent 95 \% CI band around each curve. Default \code{FALSE}.  Requires \code{pROC::ci.se()}.
#' @param ncol Integer; number of columns when \code{combine = FALSE}. Default \code{NULL} (auto).
#' @param label_size Legend text size. Default 9.
#'
#' @return A \code{ggplot} object. \code{attr(, "auc_table")} contains the same data.frame produced by \code{\link{get_gene_auc}}.
#'
#' @examples
#' \dontrun{
#'   consensus <- get_ml_consensus(ml_list, min_methods = 2)
#'   plot_gene_roc(consensus, ml_data = ml_data)
#'   plot_gene_roc(consensus, ml_data = ml_data, combine = FALSE)
#' }
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_abline geom_ribbon
#'   scale_colour_manual labs theme_bw theme element_text element_line
#'   element_rect facet_wrap
#' @importFrom rlang .data
plot_gene_roc <- function(genes,
                          ml_data = NULL,
                          expr_mat = NULL,
                          group = NULL,
                          use = c("auto", "train", "test"),
                          combine = TRUE,
                          palette = NULL,
                          show_ci = FALSE,
                          ncol = NULL,
                          label_size = 9,
                          p_adjust_method = "BH") {

  if (!requireNamespace("pROC", quietly = TRUE))
    stop("Package 'pROC' is required.", call. = FALSE)

  dat <- .resolve_expr_group(genes, ml_data = ml_data, expr_mat = expr_mat,
                             group = group, use = use)
  expr_df <- dat$expr
  grp <- dat$group
  lvls <- dat$levels
  name_map <- dat$name_map

  ## compute ROC for each gene
  frames <- list()
  ci_frames <- list()
  auc_vals <- numeric(0)

  for (col in colnames(expr_df)) {
    gname <- name_map[[col]]
    x <- expr_df[[col]]
    roc_obj <- pROC::roc(response = grp, predictor = x,
                         levels = lvls, quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    auc_vals[gname] <- auc_val

    label <- sprintf("%s (AUC = %.3f)", gname, auc_val)
    frames[[gname]] <- data.frame(
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities,
      gene = label,
      stringsAsFactors = FALSE
    )

    ## Optional CI band
    if (isTRUE(show_ci)) {
      ci_se <- tryCatch(
        pROC::ci.se(roc_obj, specificities = seq(0, 1, 0.025)),
        error = function(e) NULL
      )
      if (!is.null(ci_se)) {
        ci_mat <- as.data.frame(ci_se)
        ci_frames[[gname]] <- data.frame(
          fpr = 1 - as.numeric(rownames(ci_mat)),
          tpr_lo = ci_mat[, 1L],
          tpr_hi = ci_mat[, 3L],
          gene = label,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(frames) == 0L) {
    message("No ROC data available")
    return(invisible(NULL))
  }

  ## Sort by AUC descending for legend order
  ord <- names(sort(auc_vals, decreasing = TRUE))
  roc_df <- do.call(rbind, frames[ord])
  roc_df$gene <- factor(roc_df$gene, levels = unique(roc_df$gene))

  ## Colour palette ─ 12-colour publication palette
  n <- length(ord)
  if (is.null(palette)) {
    palette <- c("#e74c3c", "#3498db", "#2ecc71", "#9b59b6",
                 "#e67e22", "#1abc9c", "#34495e", "#e84393",
                 "#00b894", "#6c5ce7", "#fdcb6e", "#636e72")
    palette <- rep_len(palette, n)
  } else {
    palette <- rep_len(palette, n)
  }
  names(palette) <- levels(roc_df$gene)

  p <- ggplot(roc_df, aes(
    x = .data$fpr, y = .data$tpr, colour = .data$gene)) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey70", linewidth = 0.4)

  ## CI band (ribbon)
  if (isTRUE(show_ci) && length(ci_frames) > 0L) {
    ci_df <- do.call(rbind, ci_frames[ord])
    ci_df$gene <- factor(ci_df$gene, levels = levels(roc_df$gene))
    p <- p + geom_ribbon(
      data = ci_df,
      aes(x = .data$fpr, ymin = .data$tpr_lo, ymax = .data$tpr_hi,
          fill = .data$gene),
      alpha = 0.12, colour = NA, inherit.aes = FALSE
    ) + ggplot2::scale_fill_manual(values = palette, guide = "none")
  }

  p <- p +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = palette) +
    labs(x = "1 \u2212 Specificity", y = "Sensitivity",
         colour = NULL,
         subtitle = "Single-gene diagnostic ROC") +
    theme_bw(base_size = 12)

  if (isTRUE(combine)) {
    p <- p + theme(
      legend.position = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = "white", colour = "grey80",
                                       linewidth = 0.3),
      legend.text  = element_text(size = label_size),
      panel.grid.minor = element_line(linewidth = 0.2),
      plot.subtitle = element_text(size = 10, colour = "grey40")
    )
  } else {
    ## Faceted mode: one panel per gene 
    p <- p +
      facet_wrap(~ gene, ncol = ncol) +
      theme(
        legend.position  = "none",
        strip.text       = element_text(face = "italic", size = 10),
        strip.background = ggplot2::element_rect(fill = "grey96",
                                                 colour = "grey80"),
        panel.grid.minor = element_line(linewidth = 0.2),
        plot.subtitle    = element_text(size = 10, colour = "grey40")
      )
  }

  ## Attach AUC table as attribute
  auc_tbl <- get_gene_auc(genes, ml_data = ml_data, expr_mat = expr_mat,
                           group = group, use = use, ci = TRUE,
                           p_adjust_method = p_adjust_method)
  attr(p, "auc_table") <- auc_tbl

  return(p)
}


## plot_gene_boxplot

#' Two-group expression boxplots for diagnostic genes
#'
#' Draws a faceted boxplot comparing expression levels between two groups for
#' each gene, with optional jittered points, violin overlay and statistical
#' test annotation.  Designed for publication-quality figures with zero extra
#' dependencies beyond \pkg{ggplot2} and base R \code{stats}.
#'
#' @inheritParams get_gene_auc
#' @param test_method Statistical test: \code{"wilcox"} (Wilcoxon rank-sum,
#'   default) or \code{"t.test"}.
#' @param palette Character(2); fill colours for the two groups.
#'   Default \code{c("#E74C3C", "#3498DB")} (red / blue).
#' @param show_points Logical; overlay jittered data points? Default \code{TRUE}.
#' @param point_alpha Alpha transparency for jittered points. Default 0.35.
#' @param violin Logical; add a violin layer behind the boxes?
#'   Default \code{FALSE}.
#' @param ncol Integer; number of columns in \code{facet_wrap}.
#'   Default \code{NULL} (auto-computed).
#' @param scales Passed to \code{\link[ggplot2]{facet_wrap}}.
#'   Default \code{"free_y"}.
#' @param p_label How to annotate p-values: \code{"p.format"} (numeric,
#'   default), \code{"p.signif"} (stars: ns / * / ** / *** / ****), or
#'   \code{"p.adj"} (adjusted, numeric).
#' @param p_adjust_method Method for multiple-testing correction passed to
#'   \code{\link[stats]{p.adjust}}.  Default \code{"BH"}.  Applies both to the
#'   \code{p_adj} column in \code{attr(, "test_table")} and to the annotation
#'   when \code{p_label = "p.adj"}.  When only a single gene is tested,
#'   \code{p_adj} equals \code{p_value} regardless of method.
#' @param base_size Base font size for \code{\link[ggplot2]{theme_bw}}.
#'   Default 12.
#'
#' @return A \code{ggplot} object.
#'   \code{attr(, "test_table")} contains a \code{data.frame}:
#'   \code{gene}, \code{statistic}, \code{p_value}, \code{p_adj},
#'   \code{method}, \code{group1_mean}, \code{group2_mean}, \code{log2FC}.
#'
#' @examples
#' \dontrun{
#'   consensus <- get_ml_consensus(ml_list, min_methods = 2)
#'   plot_gene_boxplot(consensus, ml_data = ml_data)
#'   plot_gene_boxplot(consensus, expr_mat = expr_mat, group = group,
#'                     violin = TRUE, p_label = "p.signif")
#' }
#' @export
#' @importFrom ggplot2 ggplot aes geom_boxplot geom_jitter geom_violin
#'   facet_wrap scale_fill_manual labs theme_bw theme element_text annotate
#' @importFrom stats wilcox.test t.test p.adjust
#' @importFrom rlang .data
plot_gene_boxplot <- function(genes,
                              ml_data = NULL,
                              expr_mat = NULL,
                              group = NULL,
                              use = c("auto", "train", "test"),
                              test_method = c("wilcox", "t.test"),
                              palette = c("#E74C3C", "#3498DB"),
                              show_points = TRUE,
                              point_alpha = 0.35,
                              violin = FALSE,
                              ncol = NULL,
                              scales = "free_y",
                              p_label = c("p.format", "p.signif", "p.adj"),
                              p_adjust_method = "BH",
                              base_size = 12) {

  test_method <- match.arg(test_method)
  p_label <- match.arg(p_label)
  p_adjust_method <- match.arg(p_adjust_method,
                               choices = stats::p.adjust.methods)
  palette <- rep_len(palette, 2L)

  dat <- .resolve_expr_group(genes, ml_data = ml_data, expr_mat = expr_mat,
                             group = group, use = use)
  expr_df <- dat$expr
  grp <- dat$group
  lvls <- dat$levels
  name_map <- dat$name_map

  long_list <- lapply(colnames(expr_df), function(col) {
    data.frame(
      gene = name_map[[col]],
      group = grp,
      expression = expr_df[[col]],
      stringsAsFactors = FALSE
    )
  })
  long_df <- do.call(rbind, long_list)

  ## Preserve gene order from input
  gene_order <- unique(name_map[colnames(expr_df)])
  long_df$gene  <- factor(long_df$gene, levels = gene_order)
  long_df$group <- factor(long_df$group, levels = lvls)

  test_rows <- lapply(gene_order, function(g) {
    sub <- long_df[long_df$gene == g, ]
    vals1 <- sub$expression[sub$group == lvls[1]]
    vals2 <- sub$expression[sub$group == lvls[2]]

    res <- if (test_method == "wilcox") {
      stats::wilcox.test(vals1, vals2)
    } else {
      stats::t.test(vals1, vals2)
    }

    m1 <- mean(vals1, na.rm = TRUE)
    m2 <- mean(vals2, na.rm = TRUE)
    ## log2FC: positive = higher in group1 (positive class)
    lfc <- if (m2 > 0) log2(m1 / m2) else NA_real_

    data.frame(
      gene = g,
      statistic = unname(res$statistic),
      p_value = res$p.value,
      method = test_method,
      group1_mean = m1,
      group2_mean = m2,
      log2FC = lfc,
      stringsAsFactors = FALSE
    )
  })
  test_tbl <- do.call(rbind, test_rows)
  test_tbl$p_adj <- stats::p.adjust(test_tbl$p_value, method = p_adjust_method)

  test_tbl$label <- switch(p_label,
    p.format = paste0("p = ", .fmt_p(test_tbl$p_value)),
    p.signif = .signif_stars(test_tbl$p_value),
    p.adj    = paste0("adj.p = ", .fmt_p(test_tbl$p_adj))
  )
  test_tbl$gene <- factor(test_tbl$gene, levels = gene_order)

  y_pos <- stats::aggregate(expression ~ gene, data = long_df,
                            FUN = max, na.rm = TRUE)
  names(y_pos)[2L] <- "ymax"
  annot <- merge(test_tbl[, c("gene", "label")], y_pos, by = "gene")
  annot$y <- annot$ymax * 1.08  # bracket height
  annot$y_lab <- annot$ymax * 1.14  # text height

  p <- ggplot(long_df, aes(
    x = .data$group, y = .data$expression, fill = .data$group))

  if (isTRUE(violin)) {
    p <- p + geom_violin(alpha = 0.15, colour = NA, width = 0.9,
                         show.legend = FALSE)
  }

  ## Boxplot
  p <- p + geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85,
                        colour = "grey30", linewidth = 0.35)

  ## Jittered points
  if (isTRUE(show_points)) {
    p <- p + geom_jitter(width = 0.15, size = 1.2, alpha = point_alpha,
                         show.legend = FALSE)
  }

  ## Bracket + p-value annotation (manual segments + text)
  ## Using geom_segment + geom_text per facet via a small data.frame
  bracket_df <- data.frame(
    gene = annot$gene,
    x = 1, xend = 2,
    y = annot$y,
    label = annot$label,
    y_lab = annot$y_lab,
    stringsAsFactors = FALSE
  )
  bracket_df$gene <- factor(bracket_df$gene, levels = gene_order)

  p <- p +
    ## horizontal bracket line
    ggplot2::geom_segment(
      data = bracket_df,
      aes(x = .data$x, xend = .data$xend, y = .data$y, yend = .data$y),
      inherit.aes = FALSE, colour = "grey30", linewidth = 0.4
    ) +
    ## left tick
    ggplot2::geom_segment(
      data = bracket_df,
      aes(x = .data$x, xend = .data$x,
          y = .data$y, yend = .data$y - (.data$y_lab - .data$y) * 0.3),
      inherit.aes = FALSE, colour = "grey30", linewidth = 0.4
    ) +
    ## right tick
    ggplot2::geom_segment(
      data = bracket_df,
      aes(x = .data$xend, xend = .data$xend,
          y = .data$y, yend = .data$y - (.data$y_lab - .data$y) * 0.3),
      inherit.aes = FALSE, colour = "grey30", linewidth = 0.4
    ) +
    ## p-value text
    ggplot2::geom_text(
      data = bracket_df,
      aes(x = 1.5, y = .data$y_lab, label = .data$label),
      inherit.aes = FALSE, size = 3.5, colour = "grey20"
    )

  p <- p +
    facet_wrap(~ gene, scales = scales, ncol = ncol) +
    scale_fill_manual(values = stats::setNames(palette, lvls)) +
    labs(x = NULL, y = "Expression", fill = "Group") +
    theme_bw(base_size = base_size) +
    theme(
      strip.text       = element_text(face = "italic", size = base_size),
      strip.background = ggplot2::element_rect(fill = "grey96", colour = "grey80"),
      legend.position  = "bottom",
      panel.grid.minor = element_line(linewidth = 0.2)
    )

  ## Reorder test_table columns
  test_tbl <- test_tbl[, c("gene", "statistic", "p_value", "p_adj",
                            "method", "group1_mean", "group2_mean",
                            "log2FC", "label")]
  attr(p, "test_table") <- test_tbl

  return(p)
}
