#' Prepare expression data for ML feature selection
#'
#' Supports three modes:
#' - **Mode A** (default): full-data CV, no hold-out (`split = FALSE`).
#' - **Mode B**: internal train/test split (`split = TRUE`).
#' - **Mode C**: external validation (`test_expr` + `test_group`).
#'
#' @param expr_mat Numeric matrix (genes x samples). Row names = gene symbols.
#' @param group Factor/character of length `ncol(expr_mat)`, two levels.
#' @param positive_class Which level to treat as the positive class (`levels[[1]]`).
#'   Default `NULL`: if `group` is already an ordered factor, its level order is
#'   preserved; otherwise levels are sorted alphabetically.
#' @param genes Optional character vector of candidate genes to keep. Default is NULL.
#' @param split Logical. `TRUE` = Mode B. Default is `FALSE`.
#' @param train_ratio Fraction for training when `split = TRUE`. Default is 0.7.
#' @param train_idx Optional integer indices overriding `train_ratio`.
#' @param test_expr External validation matrix (genes x samples) for Mode C.
#' @param test_group Labels for `test_expr` (required if `test_expr` given).
#' @param seed Random seed. Default is 2025.
#'
#' @return A `tcm_ml_data` list: `train_x`, `train_y`, `test_x`, `test_y`, `gene_names`, `levels`, `full_cv`.
#' @importFrom stats var
#' @examples
#' \dontrun{
#'   ## Mode A: full CV
#'   ml_data <- prepare_ml_data(expr_mat, group, positive_class = "Disease")
#'   ## Mode B: internal split
#'   ml_data <- prepare_ml_data(expr_mat, group, split = TRUE, train_ratio = 0.7)
#' }
#' @export
prepare_ml_data <- function(expr_mat,
                            group,
                            positive_class = NULL,
                            genes = NULL,
                            split = FALSE,
                            train_ratio = 0.7,
                            train_idx = NULL,
                            test_expr = NULL,
                            test_group = NULL,
                            seed = 2025) {

  if (isTRUE(split)) {
     .check_ml_deps("caret")
   }

  if (is.data.frame(group) || is.matrix(group)) {
    if (ncol(group) == 1) {
      group <- group[[1]]
    } else if ("group" %in% colnames(group)) {
      group <- group$group
    } else {
      stop("`group` must be a vector or factor, not a data.frame or matrix.")
    }
  }

  expr_mat <- as.matrix(expr_mat)
  group <- factor(group)
  stopifnot(is.numeric(expr_mat),
            ncol(expr_mat) == length(group),
            nlevels(group) == 2L)

  ## Re-order levels so that positive_class is levels[1]
  if (!is.null(positive_class)) {
    if (!positive_class %in% levels(group))
      stop(sprintf("positive_class '%s' not found in group levels: %s",
                   positive_class, paste(levels(group), collapse = ", ")))
    other <- setdiff(levels(group), positive_class)
    group <- factor(group, levels = c(positive_class, other))
  }
  levels(group) <- make.names(levels(group))

  if (!is.null(genes)) {
    genes <- intersect(genes, rownames(expr_mat))
    if (length(genes) == 0L)
      stop("None of the supplied genes found in rownames(expr_mat).")
    expr_mat <- expr_mat[genes, , drop = FALSE]
  }

  rv <- apply(expr_mat, 1, stats::var, na.rm = TRUE)
  expr_mat <- expr_mat[rv > 0, , drop = FALSE]
  dat <- as.data.frame(t(expr_mat))
  colnames(dat) <- make.names(colnames(dat), unique = TRUE)

  ## Mode C: external validation
  if (!is.null(test_expr)) {
    if (is.null(test_group))
      stop("test_group is required when test_expr is provided.")
    test_expr <- as.matrix(test_expr)
    test_group <- factor(test_group, levels = levels(group))
    common <- intersect(rownames(expr_mat), rownames(test_expr))
    if (length(common) == 0L)
      stop("No shared genes between train and test data.")
    dat <- as.data.frame(t(expr_mat[common, , drop = FALSE]))
    test_dat <- as.data.frame(t(test_expr[common, , drop = FALSE]))
    colnames(dat) <- make.names(colnames(dat), unique = TRUE)
    colnames(test_dat) <- make.names(colnames(test_dat), unique = TRUE)

    out <- list(train_x = dat, train_y = group,
                test_x = test_dat, test_y = test_group,
                gene_names = colnames(dat), levels = levels(group),
                full_cv = FALSE)
    mode_msg <- sprintf("Mode C: Train %d | Test %d | %d genes",
                        nrow(dat), nrow(test_dat), ncol(dat))

  ## Mode B: internal split
  } else if (isTRUE(split)) {
    if (!is.null(train_idx)) {
      idx <- train_idx
    } else {
      set.seed(seed)
      idx <- caret::createDataPartition(group, p = train_ratio, list = FALSE)[, 1]
    }
    out <- list(train_x = dat[idx, , drop = FALSE],
                train_y = group[idx],
                test_x  = dat[-idx, , drop = FALSE],
                test_y  = group[-idx],
                gene_names = colnames(dat), levels = levels(group),
                full_cv = FALSE)
    mode_msg <- sprintf("Mode B: Train %d | Test %d | %d genes",
                        nrow(out$train_x), nrow(out$test_x), ncol(dat))

  ## Mode A: full-data CV (default)
  } else {
    out <- list(train_x = dat, train_y = group,
                test_x = NULL, test_y = NULL,
                gene_names = colnames(dat), levels = levels(group),
                full_cv = TRUE)
    mode_msg <- sprintf("Mode A: %d samples | %d genes",
                        nrow(dat), ncol(dat))
  }

  class(out) <- "tcm_ml_data"
  message(mode_msg)
  return(out)
}


## Internal: check required packages for ML modules
#' @keywords internal
#' @noRd
.check_ml_deps <- function(pkgs) {
  missing <- vapply(pkgs,
                    function(p) !requireNamespace(p, quietly = TRUE),
                    logical(1))
  if (any(missing))
    stop(sprintf("Package(s) required: %s\n  install.packages(c(%s))",
                 paste(pkgs[missing], collapse = ", "),
                 paste(shQuote(pkgs[missing]), collapse = ", ")),
         call. = FALSE)
}


## Internal: caret trainControl for AUC-based repeated CV
#' @param method CV method passed to [caret::trainControl()]. Default \code{"repeatedcv"}.
#' @param number Number of CV folds. Default \code{5}.
#' @param repeats Number of CV repeats. Default \code{5}.
#' @param seed Ignored (kept for signature compatibility; callers should call
#'   [set.seed()] immediately before the \code{caret::train()} / \code{caret::rfe()} call).
#' @param class_probs Logical. Whether to compute class probabilities. Default \code{TRUE}.
#' @param summary_func Summary function for [caret::trainControl()]. \code{NULL} (default)
#'   uses \code{caret::twoClassSummary} when \code{class_probs = TRUE}, otherwise
#'   \code{caret::defaultSummary}.
#' @keywords internal
#' @noRd
.make_train_ctrl <- function(method = "repeatedcv",
                             number = 5,
                             repeats = 5,
                             seed = 2025,
                             class_probs = TRUE,
                             summary_func = NULL) {
  if (is.null(summary_func)) {
    summary_func <- if (class_probs) caret::twoClassSummary else caret::defaultSummary
  }
  
  args <- list(
    method = method,
    number = number,
    classProbs = class_probs,
    summaryFunction = summary_func,
    savePredictions = "final",
    allowParallel = TRUE
  )
  
  if (method == "repeatedcv") {
    args$repeats <- repeats
  }
  
  return(do.call(caret::trainControl, args))
}


## Internal: evaluate model on hold-out test set
#' @keywords internal
#' @noRd
.eval_test <- function(model, test_x, test_y, positive_class) {
  pred_class <- stats::predict(model, newdata = test_x)
  pred_prob <- stats::predict(model, newdata = test_x, type = "prob")
  cm <- caret::confusionMatrix(pred_class, test_y, positive = positive_class)

  auc_val <- NA_real_
  if (requireNamespace("pROC", quietly = TRUE)) {
    roc_obj <- pROC::roc(response = test_y,
                         predictor = pred_prob[[positive_class]],
                         levels = levels(test_y), quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
  }

  return(list(
    predictions = pred_class,
    probabilities = pred_prob,
    confusion = cm,
    accuracy = as.numeric(cm$overall["Accuracy"]),
    auc = auc_val,
    sensitivity = as.numeric(cm$byClass["Sensitivity"]),
    specificity = as.numeric(cm$byClass["Specificity"])
  ))
}


## Internal: S3 constructor
#' @keywords internal
#' @noRd
.new_tcm_ml <- function(method, model, importance, selected_features,
                        cv_performance, test_performance = NULL,
                        ml_data) {
  return(structure(
    list(method = method,
         model = model,
         importance = importance,
         selected_features = selected_features,
         cv_performance = cv_performance,
         test_performance = test_performance,
         ml_data = ml_data),
    class = "tcm_ml"
  ))
}


#' @export
print.tcm_ml <- function(x, ...) {
  tag <- if (!is.null(x$test_performance)) "Test" else "CV"
  perf <- if (tag == "Test") x$test_performance else x$cv_performance
  .pv <- function(v) if (is.null(v) || is.na(v)) "NA" else sprintf("%.4f", v)

  cat(sprintf("=== tcm_ml: %s ===\n", toupper(x$method)))
  cat(sprintf("  Features: %d | %s AUC: %s | Sens: %s | Spec: %s\n",
              length(x$selected_features),
              tag, .pv(perf$auc), .pv(perf$sensitivity), .pv(perf$specificity)))
  invisible(x)
}

#' @export
print.tcm_ml_list <- function(x, ...) {
  cat(sprintf("=== ML Screening (%d methods) ===\n", length(x)))
  .pv <- function(v) if (is.null(v) || is.na(v)) "NA" else sprintf("%.4f", v)
  for (m in x) {
    perf <- if (!is.null(m$test_performance)) m$test_performance else m$cv_performance
    cat(sprintf("  [%s] %d features | AUC = %s\n",
                toupper(m$method), length(m$selected_features), .pv(perf$auc)))
  }
  invisible(x)
}


#' Summarise ML screening results
#'
#' Returns a tidy data.frame with one row per method, showing the number of
#' selected features and performance metrics (CV or Test AUC, Sensitivity,
#' Specificity).
#'
#' @param object A `tcm_ml_list` produced by [run_ml_screening()].
#' @param ... Ignored.
#'
#' @return A data.frame with columns: `method`, `n_features`, `auc_type`,
#'   `auc`, `auc_sd`, `sensitivity`, `specificity`.
#' @examples
#' \dontrun{
#'   ml_data <- prepare_ml_data(expr_mat, group)
#'   res_list <- run_ml_screening(ml_data)
#'   summary(res_list)
#' }
#' @export
summary.tcm_ml_list <- function(object, ...) {
  .safe <- function(v) if (is.null(v) || length(v) == 0 || is.na(v)) NA_real_ else v

  rows <- lapply(object, function(m) {
    has_test <- !is.null(m$test_performance)
    perf <- if (has_test) m$test_performance else m$cv_performance
    cv_sd <- if (!has_test && !is.null(m$cv_performance$auc_sd))
               m$cv_performance$auc_sd else NA_real_

    data.frame(
      method = m$method,
      n_features = length(m$selected_features),
      auc_type = if (has_test) "Test" else "CV",
      auc = .safe(perf$auc),
      auc_sd = cv_sd,
      sensitivity = .safe(perf$sensitivity),
      specificity = .safe(perf$specificity),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}


#' Re-select top features from a fitted ML model
#' After running a model (especially Ridge or XGBoost where all features
#' are initially retained), inspect \code{result$importance} to decide
#' how many to keep, then call this function to trim the selection.
#'
#' @param ml_obj A \code{tcm_ml} object from any \code{ml_*} function.
#' @param top_n Integer; the number of top features to keep (ranked by importance).
#' @return A modified \code{tcm_ml} object with updated
#'   \code{$selected_features} and \code{$genes}.
#' @details The underlying model and importance table are unchanged.
#'   Only the gene selection is trimmed. Performance metrics still
#'   reflect the original run; re-run the model with \code{top_n}
#'   if you need updated metrics.
#'
#' @examples
#' \dontrun{
#'   xgb <- ml_xgboost(ml_data)        # keep all features as default
#'   head(xgb$importance, 20)           # inspect ranking
#'   xgb <- select_features(xgb, 15)   # keep top 15
#' }
#' @export
select_features <- function(ml_obj, top_n) {
  stopifnot(inherits(ml_obj, "tcm_ml"))
  stopifnot(is.numeric(top_n), length(top_n) == 1L, top_n >= 1)

  imp <- ml_obj$importance
  if (is.null(imp) || nrow(imp) == 0L)
    stop("No importance data found in this model object.")

  top_n <- min(as.integer(top_n), nrow(imp))
  new_genes <- imp$gene[seq_len(top_n)]

  ml_obj$selected_features <- new_genes
  ml_obj$genes <- new_genes

  message(sprintf("[%s] Re-selected top %d / %d features",
                  toupper(ml_obj$method), top_n, nrow(imp)))
  return(ml_obj)
}
