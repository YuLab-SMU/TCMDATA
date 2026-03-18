## Elastic Net / LASSO / Ridge plots

#' Plot Elastic Net / LASSO / Ridge CV curve
#'
#' Wrapper around the default \code{plot.cv.glmnet()} method from the
#' \pkg{glmnet} package. Displays the CV metric (AUC / deviance) as a
#' function of \eqn{\log(\lambda)}, with ± 1 SE error bars and vertical
#' dashed lines at \code{lambda.min} and \code{lambda.1se}.
#'
#' @param result A \code{tcm_ml} object from [ml_enet()], [ml_lasso()] or [ml_ridge()].
#' @param ... Additional arguments passed to \code{plot.cv.glmnet()}.
#' @return Invisible \code{NULL}. Called for its side effect (a base R plot).
#' @export
plot_enet_cv <- function(result, ...) {
  if (!inherits(result, "tcm_ml") || is.null(result$cv_fit))
    stop("result must be a tcm_ml object from ml_enet() / ml_lasso() / ml_ridge()")
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("Package 'glmnet' is required")
  graphics::plot(result$cv_fit, ...)
  invisible(NULL)
}

#' Plot Elastic Net / LASSO / Ridge coefficient path
#'
#' Wrapper around the default \code{plot.glmnet()} method.
#' Each curve traces one feature's coefficient across the \eqn{\lambda} grid.
#' When \code{top_n} is set, only the features with the largest absolute
#' coefficients at the selected \eqn{\lambda} are shown.
#'
#' @param result A \code{tcm_ml} object from [ml_enet()], [ml_lasso()] or [ml_ridge()].
#' @param top_n Integer or \code{NULL}. If non-\code{NULL} and fewer features
#'   than the total, only the \code{top_n} features with largest absolute
#'   coefficient at the selected \eqn{\lambda} are drawn. Default 20.
#' @param xvar Variable on x-axis: \code{"lambda"} (default), \code{"norm"}, or \code{"dev"}.
#' @param ... Additional arguments passed to \code{plot.glmnet()}.
#' @return Invisible \code{NULL}. Called for its side effect (a base R plot).
#' @export
#' @importFrom graphics plot abline
plot_enet_path <- function(result, top_n = 20, xvar = "lambda", ...) {
  if (!inherits(result, "tcm_ml") || is.null(result$full_fit))
    stop("result must be a tcm_ml object from ml_enet() / ml_lasso() / ml_ridge()")
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("Package 'glmnet' is required")

  fit <- result$full_fit

  ## Subset to top_n features by |coef| at selected lambda

  if (!is.null(top_n) && nrow(fit$beta) > top_n) {
    idx <- which.min(abs(fit$lambda - result$lambda_used))
    coef_at_sel <- abs(as.matrix(fit$beta)[, idx])
    keep <- names(sort(coef_at_sel, decreasing = TRUE))[seq_len(top_n)]
    fit$beta <- fit$beta[keep, , drop = FALSE]
  }

  graphics::plot(fit, xvar = xvar, ...)
  if (xvar == "lambda") {
    graphics::abline(v = log(result$lambda_used),
                     lty = 2, col = "red")
  }
  invisible(NULL)
}

#' Bar chart of LASSO / Elastic Net / Ridge coefficients
#'
#' Horizontal bar chart coloured by sign (positive / negative).
#'
#' @param result A \code{tcm_ml} object from [ml_enet()], [ml_lasso()] or
#'   [ml_ridge()].
#' @param top_n Number of top features to show. Default 20.
#' @return A \code{ggplot} object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_col coord_flip scale_fill_manual
#'   labs theme_bw
#' @importFrom rlang .data
plot_enet_coefs <- function(result, top_n = 20) {
  if (!inherits(result, "tcm_ml") || is.null(result$coefficients))
    stop("result must be a tcm_ml object from ml_enet() / ml_lasso() / ml_ridge()")

  df <- result$coefficients
  df <- df[df$coefficient != 0, , drop = FALSE]
  df <- df[order(-df$abs_coef), ]
  if (nrow(df) > top_n) df <- df[seq_len(top_n), ]
  df$gene <- factor(df$gene, levels = rev(df$gene))
  df$sign <- ifelse(df$coefficient > 0, "Positive", "Negative")

  p <- ggplot(df, aes(
    x = .data$gene, y = .data$coefficient, fill = .data$sign
  )) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c(Positive = "#e74c3c", Negative = "#3498db")) +
    labs(x = NULL, y = "Coefficient", fill = "Direction") +
    theme_bw()
  return(p)
}


## RF plots

#' Plot Boruta feature selection result
#'
#' @param result A \code{tcm_ml} object from [ml_rf()].
#' @param top_n Integer; show at most this many features.
#' @return A \code{ggplot} object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_col geom_vline scale_fill_manual
#'   labs theme_bw
#' @importFrom rlang .data
plot_rf_boruta <- function(result, top_n = 20) {
  if (!inherits(result, "tcm_ml") || is.null(result$boruta))
    stop("result must be a tcm_ml object from ml_rf()")
  if (!requireNamespace("Boruta", quietly = TRUE))
    stop("Package 'Boruta' is required")

  imp <- Boruta::attStats(result$boruta)
  imp$feature <- rownames(imp)
  imp <- imp[order(imp$medianImp, decreasing = TRUE), ]

  if (!is.null(top_n) && nrow(imp) > top_n)
    imp <- imp[seq_len(top_n), ]

  imp$feature <- factor(imp$feature, levels = rev(imp$feature))

  shadow_ref <- stats::median(
    result$boruta$ImpHistory[, "shadowMax"], na.rm = TRUE
  )

  p <- ggplot(imp, aes(
    x = .data$medianImp, y = .data$feature, fill = .data$decision
  )) +
    geom_col() +
    geom_vline(
      xintercept = shadow_ref, linetype = "dashed", colour = "grey40"
    ) +
    scale_fill_manual(
      values = c(Confirmed = "#2ecc71", Tentative = "#f1c40f",
                 Rejected = "#e74c3c")
    ) +
    labs(x = "Importance (median)", y = NULL, fill = "Decision") +
    theme_bw()
  return(p)
}

#' Plot RF variable importance
#'
#' @param result A \code{tcm_ml} object from [ml_rf()].
#' @param top_n Integer; number of top features to show.
#' @return A \code{ggplot} object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_col coord_flip labs theme_bw
#' @importFrom rlang .data
plot_rf_importance <- function(result, top_n = 20) {
  if (!inherits(result, "tcm_ml") || is.null(result$importance))
    stop("result must be a tcm_ml object from ml_rf()")

  df <- result$importance
  df <- df[order(df$MeanDecreaseGini, decreasing = TRUE), ]
  if (nrow(df) > top_n) df <- df[seq_len(top_n), ]
  df$gene <- factor(df$gene, levels = rev(df$gene))

  p <- ggplot(df, aes(
    x = .data$gene, y = .data$MeanDecreaseGini
  )) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(x = NULL, y = "Mean Decrease Gini") +
    theme_bw()
  return(p)
}


## SVM-RFE plots

#' Plot SVM-RFE accuracy profile
#' @param result A \code{tcm_ml} object from [ml_svm_rfe()].
#' @return A \code{ggplot} object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_vline labs theme_bw
#' @importFrom rlang .data
plot_svm_rfe_curve <- function(result) {
  if (!inherits(result, "tcm_ml") || is.null(result$profile))
    stop("result must be a tcm_ml object from ml_svm_rfe()")

  df <- result$profile
  best_n <- length(result$genes)

  metric <- if ("ROC" %in% names(df)) "ROC" else "Accuracy"
  df$.metric <- df[[metric]]

  p <- ggplot(df, aes(
    x = .data$Variables, y = .data$.metric
  )) +
    geom_line(colour = "steelblue") +
    geom_point(colour = "steelblue") +
    geom_vline(
      xintercept = best_n, linetype = "dashed", colour = "red"
    ) +
    labs(x = "Number of features", y = paste0(metric, " (CV)")) +
    theme_bw()
  return(p)
}


## XGBoost plots

#' Plot XGBoost feature importance (Gain)
#'
#' @param result A \code{tcm_ml} object from [ml_xgboost()].
#' @param top_n Integer; number of top features to show.
#' @return A \code{ggplot} object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_col coord_flip labs theme_bw
#' @importFrom rlang .data
plot_xgb_importance <- function(result, top_n = 20) {
  if (!inherits(result, "tcm_ml") || is.null(result$importance))
    stop("result must be a tcm_ml object from ml_xgboost()")

  df <- result$importance
  df <- df[order(-df$importance), ]
  if (nrow(df) > top_n) df <- df[seq_len(top_n), ]
  df$gene <- factor(df$gene, levels = rev(df$gene))

  p <- ggplot(df, aes(
    x = .data$gene, y = .data$importance
  )) +
    geom_col(fill = "#e67e22") +
    coord_flip() +
    labs(x = NULL, y = "Gain") +
    theme_bw()
  return(p)
}

## Multi-method plots

#' ROC curves for ML methods
#'
#' Supports all modes.
#' \itemize{
#'   \item Mode B / C: ROC is computed on the held-out test set.
#'   \item Mode A: uses out-of-fold (LASSO / Enet / XGBoost / SVM-RFE) or OOB (RF)
#'         predictions.
#'         SVM-RFE falls back to resubstitution when OOF is unavailable.
#' }
#'
#' @param ml_list A \code{tcm_ml_list} from [run_ml_screening()].
#' @return A \code{ggplot} object, or \code{NULL} if no data available.
#' @export
#' @importFrom ggplot2 ggplot aes geom_line geom_abline scale_colour_manual
#'   labs theme_bw theme element_text element_line element_rect
#' @importFrom rlang .data
plot_ml_roc <- function(ml_list) {
  if (!inherits(ml_list, "tcm_ml_list"))
    stop("ml_list must be a tcm_ml_list object")
  if (!requireNamespace("pROC", quietly = TRUE))
    stop("Package 'pROC' is required")

  ml_data <- attr(ml_list, "ml_data")
  if (is.null(ml_data)) ml_data <- ml_list[[1]]$ml_data

  has_test <- !is.null(ml_data$test_x)

  frames <- list()

  for (nm in names(ml_list)) {
    res <- ml_list[[nm]]
    genes <- res$genes
    if (length(genes) == 0) next

    pred_prob <- tryCatch({
      if (has_test) {
        ## Mode B / C: predict on test set 
        if (nm %in% c("lasso", "enet", "ridge")) {
          sub <- ml_data$test_x
          as.numeric(stats::predict(
            res$cv_fit, newx = as.matrix(sub),
            s = res$lambda_used, type = "response"
          ))
        } else if (nm == "rf") {
          sub <- ml_data$test_x[, genes, drop = FALSE]
          as.numeric(stats::predict(
            res$rf_fit, newdata = sub, type = "prob"
          )[, ml_data$levels[1]])
        } else if (nm == "svm_rfe") {
          sub <- ml_data$test_x[, genes, drop = FALSE]
          as.numeric(stats::predict(
            res$rfe_result, newdata = sub
          )[, ml_data$levels[1]])
        } else if (nm == "xgboost") {
          sub <- ml_data$test_x
          as.numeric(stats::predict(res$xgb_fit, newdata = as.matrix(sub)))
        } else { NULL }

      } else {
        ## Mode A: use OOF / OOB predictions 
        if (nm %in% c("lasso", "enet", "ridge", "xgboost", "svm_rfe")) {
          ## oof_prob stored as P(levels[1])
          if (!is.null(res$oof_prob)) res$oof_prob else NULL
        } else if (nm == "rf") {
          ## OOB vote probabilities
          if (!is.null(res$rf_fit$votes)) {
            as.numeric(res$rf_fit$votes[, ml_data$levels[1]])
          } else { NULL }
        } else { NULL }
      }
    }, error = function(e) NULL)

    if (is.null(pred_prob)) next

    y_true <- if (has_test) ml_data$test_y else ml_data$train_y
    roc_obj <- pROC::roc(y_true, pred_prob, levels = ml_data$levels,
                         quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))

    coords <- data.frame(
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities,
      method = sprintf("%s (AUC = %.3f)", toupper(nm), auc_val),
      stringsAsFactors = FALSE
    )
    frames[[nm]] <- coords
  }

  if (length(frames) == 0) {
    message("No ROC data available")
    return(invisible(NULL))
  }

  roc_df <- do.call(rbind, frames)

  ## Build subtitle
  if (has_test) {
    subtitle <- "Test-set ROC"
  } else {
    tags <- vapply(names(frames), function(nm) {
      if (nm == "rf") "OOB" else "OOF"
    }, character(1))
    subtitle <- paste0("Mode A ROC (", paste(unique(tags), collapse = " / "), ")")
  }

  ## Colour palette
  n_method <- length(frames)
  pal <- c("#e74c3c", "#3498db", "#2ecc71", "#9b59b6",
           "#e67e22", "#1abc9c")[seq_len(n_method)]
  names(pal) <- unique(roc_df$method)

  p <- ggplot(roc_df, aes(
    x = .data$fpr, y = .data$tpr, colour = .data$method
  )) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dashed", colour = "grey70", linewidth = 0.4) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = pal) +
    labs(x = "1 \u2212 Specificity", y = "Sensitivity",
         colour = NULL, subtitle = subtitle) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = "white", colour = "grey80",
                                       linewidth = 0.3),
      legend.text = element_text(size = 9),
      panel.grid.minor = element_line(linewidth = 0.2),
      plot.subtitle = element_text(size = 10, colour = "grey40")
    )

  return(p)
}


#' Extract gene sets from ML results
#'
#' Flexibly accepts the output of [run_ml_screening()], individually fitted
#' \code{tcm_ml} objects, or plain character vectors, and returns a named list
#' of character vectors ready for downstream plotting.
#'
#' @param ... One of the following:
#' \describe{
#'   \item{A single \code{tcm_ml_list}}{from [run_ml_screening()].}
#'   \item{Multiple \code{tcm_ml} objects}{e.g.
#'     \code{get_ml_gene_sets(res_lasso, res_rf, res_xgb)}.}
#'   \item{A single named list}{of character vectors or \code{tcm_ml} objects.}
#'   \item{Multiple character vectors}{gene sets directly.}
#' }
#' @param set_names Optional character vector of names for the sets.
#'   Ignored when names can be inferred from the input.
#'
#' @return A named list of character vectors (one per method / set).
#'   Can be passed directly to:
#'   \itemize{
#'     \item \code{do.call(getvenndata, c(gene_sets, list(set_names = names(gene_sets))))}
#'       then \code{ggvenn_plot()} (\eqn{\le 4} sets).
#'     \item \code{aplotExtra::upset_plot(gene_sets)} (\eqn{> 4} sets).
#'   }
#'
#' @examples
#' \dontrun{
#'   ## From run_ml_screening()
#'   gene_sets <- get_ml_gene_sets(ml_list)
#'
#'   ## From individual models
#'   gene_sets <- get_ml_gene_sets(res_lasso, res_rf, res_xgb)
#'
#'   ## Venn diagram (2-4 sets)
#'   venn_df <- do.call(getvenndata,
#'                      c(gene_sets, list(set_names = names(gene_sets))))
#'   ggvenn_plot(venn_df)
#'
#'   ## UpSet plot (usually more than 4 sets)
#'   aplotExtra::upset_plot(gene_sets)
#' }
#' @export
get_ml_gene_sets <- function(..., set_names = NULL) {

  args <- list(...)

  if (length(args) == 1L && is.list(args[[1]]) &&
      !inherits(args[[1]], "tcm_ml")) {
    ## Single list-like input: tcm_ml_list or plain named list
    input <- args[[1]]
  } else {
    ## Multiple positional arguments
    input <- args
  }

  ## Convert each element to a character vector of genes
  gene_sets <- lapply(input, function(x) {
    if (inherits(x, "tcm_ml")) return(x$genes)
    if (is.character(x)) return(x)
    stop("Each element must be a tcm_ml object or a character vector of genes.",
         call. = FALSE)
  })

  ## Assign names
  if (!is.null(set_names)) {
    if (length(set_names) != length(gene_sets))
      stop("length(set_names) must equal the number of gene sets.", call. = FALSE)
    names(gene_sets) <- set_names
  } else if (is.null(names(gene_sets)) || any(names(gene_sets) == "")) {
    ## Try to infer from tcm_ml$method
    inferred <- vapply(input, function(x) {
      if (inherits(x, "tcm_ml")) toupper(x$method) else ""
    }, character(1))
    if (all(nchar(inferred) > 0)) {
      names(gene_sets) <- inferred
    } else {
      names(gene_sets) <- paste0("Set", seq_along(gene_sets))
    }
  }

  ## Drop empty sets
  gene_sets <- gene_sets[lengths(gene_sets) > 0]
  if (length(gene_sets) < 2)
    stop("Need at least 2 non-empty gene sets.", call. = FALSE)

  return(gene_sets)
}


#' Venn diagram of selected genes
#'
#' Convenience wrapper: calls [get_ml_gene_sets()] to extract gene sets,
#' then draws a Venn diagram via [getvenndata()] + [ggvenn_plot()].
#' Accepts the same flexible inputs as [get_ml_gene_sets()].
#' For > 4 sets consider using \code{aplotExtra::upset_plot(get_ml_gene_sets(...))}.
#'
#' @inheritParams get_ml_gene_sets
#' @return A \code{ggplot} object.
#' @export
plot_ml_venn <- function(..., set_names = NULL) {

  gene_sets <- get_ml_gene_sets(..., set_names = set_names)

  if (length(gene_sets) > 4) {
    message("ggvenn_plot supports at most 4 sets; using the first 4.\n",
            "Tip: for all sets use aplotExtra::upset_plot(get_ml_gene_sets(...))")
    gene_sets <- gene_sets[seq_len(4)]
  }

  venn_df <- do.call(getvenndata,
                     c(gene_sets, list(set_names = names(gene_sets))))
  ggvenn_plot(venn_df)
}
