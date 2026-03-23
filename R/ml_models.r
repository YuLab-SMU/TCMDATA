## Elastic Net family

#' Elastic Net feature selection via cv.glmnet
#' Penalised logistic regression with configurable alpha.
#' \itemize{
#'   \item alpha = 1 : LASSO
#'   \item alpha = 0 : Ridge
#'   \item 0 < alpha < 1 : Elastic Net
#' }
#' Caution: For Ridge (\code{alpha = 0}), all coefficients are non-zero.
#' When \code{top_n = NULL} (default), all features are retained so that
#' users can inspect \code{$importance} first, then call
#' [select_features()] to trim to the desired number.
#'
#' @param ml_data A \code{tcm_ml_data} object from [prepare_ml_data()].
#' @param alpha Mixing parameter, 0 to 1. Default 0.5.
#' @param lambda_rule \code{"1se"} (parsimonious) or \code{"min"} (best AUC).
#' @param top_n For Ridge: max number of top genes to keep (by
#'   \code{abs(coefficient)}). Default \code{NULL} (keep all).
#'   Ignored for LASSO / Elastic Net.
#' @param type_measure Loss function for CV in [glmnet::cv.glmnet()].
#'   Default \code{"auc"}. Other options: \code{"deviance"},
#'   \code{"class"}, \code{"mse"}, \code{"mae"}.
#' @param cv_folds Number of CV folds. Default 10.
#' @param seed Random seed. Default 2025.
#' @param ... Passed to [glmnet::cv.glmnet()].
#'
#' @return A \code{tcm_ml} object.
#' @importFrom stats coef predict
#' @examples
#' \dontrun{
#'   ml_data <- prepare_ml_data(expr_mat, group)
#'   res <- ml_enet(ml_data)
#' }
#' @export
ml_enet <- function(ml_data,
                    alpha = 0.5,
                    lambda_rule = c("1se", "min"),
                    top_n = NULL,
                    type_measure = "auc",
                    cv_folds = 10,
                    seed = 2025,
                    ...) {

  .check_ml_deps("glmnet")
  stopifnot(inherits(ml_data, "tcm_ml_data"))
  lambda_rule <- match.arg(lambda_rule)
  stopifnot(is.numeric(alpha), alpha >= 0, alpha <= 1)

  method_name <- if (alpha == 1) "lasso" else if (alpha == 0) "ridge" else "enet"

  x <- as.matrix(ml_data$train_x)
  y <- ml_data$train_y

  ## Relevel so that glmnet treats levels[1] (our positive class) as

  ## its event class (= the second factor level).  This makes
  ## predict(type = "response") return P(levels[1]) directly and
  ## avoids scattered 1-p inversions elsewhere.
  y_glm <- stats::relevel(y, ref = ml_data$levels[2])

  set.seed(seed)
  message(sprintf("Running %s (cv.glmnet, alpha = %g) ...",
                  toupper(method_name), alpha))
  cvfit <- glmnet::cv.glmnet(x = x, y = y_glm, family = "binomial",
                              alpha = alpha, nfolds = cv_folds,
                              type.measure = type_measure,
                              keep = TRUE, ...)

  lambda_sel <- if (lambda_rule == "1se") cvfit$lambda.1se else cvfit$lambda.min

  coefs <- stats::coef(cvfit, s = lambda_sel)
  coef_df <- data.frame(gene = rownames(coefs),
                        coefficient = as.numeric(as.matrix(coefs)),
                        stringsAsFactors = FALSE)
  coef_df <- coef_df[coef_df$gene != "(Intercept)", , drop = FALSE]
  coef_df$abs_coef <- abs(coef_df$coefficient)

  ## Feature selection 
  if (alpha == 0) {
    ## Ridge: all coefs non-zero -> rank by |coef| and take top_n
    coef_df <- coef_df[order(-coef_df$abs_coef), ]
    n_keep <- if (!is.null(top_n)) min(top_n, nrow(coef_df)) else nrow(coef_df)
    selected <- coef_df$gene[seq_len(n_keep)]
  } else {
    ## LASSO / Elastic Net: non-zero coefficients
    selected <- coef_df$gene[coef_df$coefficient != 0]
  }

  imp_df <- coef_df[order(-coef_df$abs_coef), c("gene", "abs_coef")]
  names(imp_df)[2] <- "importance"
  rownames(imp_df) <- imp_df$gene

  idx_lam <- which.min(abs(cvfit$lambda - lambda_sel))

  ## CV Sens / Spec from out-of-fold predictions
  ## fit.preval is n x nlambda on the linear predictor (logit) scale;
  ## after releveling, plogis() converts to P(levels[1]) directly.
  oof_prob <- stats::plogis(cvfit$fit.preval[, idx_lam])
  oof_cls <- factor(
    ifelse(oof_prob > 0.5, ml_data$levels[1], ml_data$levels[2]),
    levels = ml_data$levels
  )

  pos <- ml_data$levels[1]; neg <- ml_data$levels[2]
  tp <- sum(oof_cls == pos & y == pos)
  fn <- sum(oof_cls == neg & y == pos)
  tn <- sum(oof_cls == neg & y == neg)
  fp <- sum(oof_cls == pos & y == neg)

  cv_perf <- list(auc = cvfit$cvm[idx_lam], auc_sd = cvfit$cvsd[idx_lam],
                  sensitivity = tp / (tp + fn), specificity = tn / (tn + fp),
                  lambda = lambda_sel, lambda_rule = lambda_rule,
                  alpha = alpha, nfolds = cv_folds)

  test_perf <- NULL
  if (!isTRUE(ml_data$full_cv) && !is.null(ml_data$test_x)) {
    .check_ml_deps("caret")
    ctrl <- .make_train_ctrl(method = "cv", number = cv_folds,
                             repeats = 1, seed = seed)
    set.seed(seed)
    fit_caret <- caret::train(
      x = ml_data$train_x, y = ml_data$train_y,
      method = "glmnet", trControl = ctrl,
      tuneGrid = data.frame(alpha = alpha, lambda = lambda_sel),
      metric = "ROC", preProcess = c("center", "scale"),
      family = "binomial"
    )
    test_perf <- .eval_test(fit_caret, ml_data$test_x,
                            ml_data$test_y, ml_data$levels[1])
  }

  obj <- .new_tcm_ml(method = method_name, model = cvfit,
                     importance = imp_df, selected_features = selected,
                     cv_performance = cv_perf, test_performance = test_perf,
                     ml_data = ml_data)
  obj$cv_fit <- cvfit
  obj$full_fit <- cvfit$glmnet.fit
  obj$lambda_used <- lambda_sel
  obj$coefficients <- coef_df
  obj$genes <- selected
  ## Store OOF probability as P(levels[1]) for ROC plotting
  ## (already P(levels[1]) after releveling — no inversion needed)
  obj$oof_prob <- oof_prob

  sel_msg <- if (alpha == 0)
    sprintf("Top %d / %d (by |coef|)", length(selected), ncol(x))
  else
    sprintf("Non-zero: %d / %d", length(selected), ncol(x))

  message(sprintf("  lambda.%s = %.5f | %s",
                  lambda_rule, lambda_sel, sel_msg))
  return(obj)
}


#' LASSO feature selection (a wrapper function from [ml_enet()])
#' Calls [ml_enet()] with \code{alpha = 1}.
#' @inheritParams ml_enet
#' @return A \code{tcm_ml} object (method = "lasso").
#' @examples
#' \dontrun{
#'   ml_data <- prepare_ml_data(expr_mat, group)
#'   res <- ml_lasso(ml_data)
#' }
#' @export
ml_lasso <- function(ml_data,
                     lambda_rule = c("1se", "min"),
                     cv_folds = 10,
                     seed = 2025,
                     ...) {
  return(ml_enet(ml_data, alpha = 1, lambda_rule = lambda_rule,
                 cv_folds = cv_folds, seed = seed, ...))
}


#' Ridge regression feature selection (a wrapper function from [ml_enet()])
#' Calls [ml_enet()] with \code{alpha = 0}.
#' Caution: Ridge does not shrink coefficients to zero; when \code{top_n = NULL}
#' (default), all features are kept. Use [select_features()] after inspecting importance to trim.
#' @inheritParams ml_enet
#' @return A \code{tcm_ml} object (method = "ridge").
#' @examples
#' \dontrun{
#'   ml_data <- prepare_ml_data(expr_mat, group)
#'   res <- ml_ridge(ml_data)
#' }
#' @export
ml_ridge <- function(ml_data,
                     lambda_rule = c("1se", "min"),
                     top_n = NULL,
                     cv_folds = 10,
                     seed = 2025,
                     ...) {
  return(ml_enet(ml_data, alpha = 0, lambda_rule = lambda_rule,
                 top_n = top_n, cv_folds = cv_folds, seed = seed, ...))
}


#' Extract ranked coefficients from LASSO / Elastic Net / Ridge
#' @param ml_obj A \code{tcm_ml} object from [ml_lasso()], [ml_enet()], or [ml_ridge()].
#' @param top_n Return only top-n genes (optional).
#' @param min_coef Minimum absolute coefficient. Default 0.
#' @param sort_by \code{"abs"} or \code{"coef"}. Default \code{"abs"}.
#' @param include_zero Logical. Whether to include zero-coefficient genes. Default \code{FALSE}.
#'
#' @return A data.frame: \code{gene}, \code{coefficient}, \code{abs_coef}.
#' @importFrom utils head
#' @examples
#' \dontrun{
#'   ml_data <- prepare_ml_data(expr_mat, group)
#'   res <- ml_lasso(ml_data)
#'   get_enet_coefs(res, top_n = 20)
#' }
#' @export
get_enet_coefs <- function(ml_obj,
                           top_n = NULL,
                           min_coef = 0,
                           sort_by = c("abs", "coef"),
                           include_zero = FALSE) {

  stopifnot(inherits(ml_obj, "tcm_ml"))
  if (!tolower(ml_obj$method) %in% c("lasso", "enet", "ridge"))
    stop("get_enet_coefs() requires a LASSO / Elastic Net / Ridge model.")
  sort_by <- match.arg(sort_by)

  df <- ml_obj$coefficients
  if (is.null(df)) stop("No coefficient data found.")

  if (!isTRUE(include_zero))
    df <- df[abs(df$coefficient) > min_coef, , drop = FALSE]

  df <- if (sort_by == "abs") df[order(-df$abs_coef), ] else df[order(-df$coefficient), ]
  if (!is.null(top_n)) df <- utils::head(df, top_n)
  rownames(df) <- NULL
  return(df)
}


## Random Forest and Boruta

#' Random Forest + Boruta feature selection for key-gene screening
#' @description
#' Supervised feature-selection pipeline designed for network-pharmacology
#' workflows: starting from PPI-derived candidate genes and RNA-seq expression
#' profiles, this function identifies a stable, interpretable set of key genes
#' rather than optimising prediction accuracy.
#'
#' Pipeline overview
#' \enumerate{
#'   \item Train a Random Forest on **all** candidate genes to obtain
#'         importance scores (MeanDecreaseGini / MeanDecreaseAccuracy).
#'   \item Run Boruta all-relevant feature selection; tentative features are
#'         resolved via \code{\link[Boruta]{TentativeRoughFix}}.
#'   \item *(Optional, default)* Re-fit a Random Forest on the **selected**
#'         genes only so that OOB and test-set metrics honestly reflect the
#'         chosen subset (\code{refit_on_selected = TRUE}).
#' }
#'
#' @section Positive class:
#' \code{ml_data$levels[1]} is treated as the positive class throughout
#' (Sensitivity, Specificity, AUC direction).
#' Set it via the \code{positive_class} argument in \code{\link{prepare_ml_data}}.
#'
#' @param ml_data A \code{tcm_ml_data} object produced by \code{\link{prepare_ml_data}}. Must contain exactly two classes.
#' @param n_trees Number of trees for both the initial RF and the refit RF. Default \code{500}.
#' @param top_n Hard cap: keep at most \code{top_n} Boruta-confirmed genes (ranked by MeanDecreaseGini). \code{NULL} = no cap (default).
#' @param boruta_fallback_n When Boruta confirms zero features, fall back to the top-\code{boruta_fallback_n} genes ranked by MeanDecreaseGini.
#'   Default \code{20}.
#' @param max_runs Maximum Boruta iterations passed to \code{\link[Boruta]{Boruta}}. Default \code{100}.
#' @param refit_on_selected Logical. Re-fit a RF on selected genes only so that OOB / test metrics reflect the actual feature subset. Default \code{TRUE} (recommended).
#' @param seed Random seed. Default \code{2025}.
#' @param ... Additional arguments forwarded to \code{\link[randomForest]{randomForest}}.
#' @return A \code{tcm_ml} object after feature selection.
#' @importFrom stats predict
#' @examples
#' \dontrun{
#'   ml_data <- prepare_ml_data(expr_mat, group)
#'   res <- ml_rf(ml_data)
#' }
#' @export
ml_rf <- function(ml_data,
                  n_trees = 500,
                  top_n = NULL,
                  boruta_fallback_n = 20,
                  max_runs = 100,
                  refit_on_selected = TRUE,
                  seed = 2025,
                  ...) {

  ## Dependencies & input validation
  .check_ml_deps(c("randomForest", "Boruta"))
  if (!inherits(ml_data, "tcm_ml_data"))
    stop("'ml_data' must be a <tcm_ml_data> object from prepare_ml_data().",
         call. = FALSE)
  if (nlevels(ml_data$train_y) != 2L)
    stop(sprintf("ml_rf() requires exactly 2 classes, got %d: %s",
                 nlevels(ml_data$train_y),
                 paste(levels(ml_data$train_y), collapse = ", ")),
         call. = FALSE)

  stopifnot(is.numeric(n_trees), length(n_trees) == 1L, n_trees >= 1L)
  stopifnot(is.numeric(max_runs), length(max_runs) == 1L, max_runs >= 1L)
  stopifnot(is.numeric(boruta_fallback_n), length(boruta_fallback_n) == 1L,
            boruta_fallback_n >= 1L)
  stopifnot(is.logical(refit_on_selected), length(refit_on_selected) == 1L)
  if (!is.null(top_n)) {
    stopifnot(is.numeric(top_n), length(top_n) == 1L, top_n >= 1L)
    top_n <- as.integer(top_n)
  }
  n_trees <- as.integer(n_trees)
  boruta_fallback_n <- as.integer(boruta_fallback_n)

  pos <- ml_data$levels[1]
  neg <- ml_data$levels[2]

  ## Full-feature Random Forest (importance source)
  set.seed(seed)
  message("Phase 1: Training Random Forest on all features ...")
  rf_full <- randomForest::randomForest(
    x = ml_data$train_x, y = ml_data$train_y,
    ntree = n_trees, importance = TRUE, ...)

  imp_raw <- randomForest::importance(rf_full)
  imp_df <- data.frame(
    gene = rownames(imp_raw),
    MeanDecreaseAccuracy = imp_raw[, "MeanDecreaseAccuracy"],
    MeanDecreaseGini = imp_raw[, "MeanDecreaseGini"],
    stringsAsFactors = FALSE)
  imp_df <- imp_df[order(-imp_df$MeanDecreaseGini), ]
  rownames(imp_df) <- imp_df$gene

  ## OOB performance of full model (always recorded)
  oob_full <- .compute_oob(rf_full, ml_data$train_y, pos, neg, n_trees = n_trees)

  ## Boruta feature selection
  set.seed(seed)
  message("Phase 2: Running Boruta feature selection ...")
  boruta_res <- Boruta::Boruta(
    x = ml_data$train_x, y = ml_data$train_y,
    maxRuns = max_runs, ntree = n_trees, doTrace = 0)
  boruta_res <- Boruta::TentativeRoughFix(boruta_res)

  confirmed <- names(
    boruta_res$finalDecision[boruta_res$finalDecision == "Confirmed"])
  ## Keep Gini order for downstream ranking
  selected <- imp_df$gene[imp_df$gene %in% confirmed]
  boruta_used_fallback <- FALSE

  if (length(selected) == 0L) {
    n_fb <- min(boruta_fallback_n, nrow(imp_df))
    warning(sprintf(
      "Boruta confirmed 0 features; falling back to top %d by MeanDecreaseGini.",
      n_fb), call. = FALSE)
    selected <- utils::head(imp_df$gene, n_fb)
    boruta_used_fallback <- TRUE
  }

  if (!is.null(top_n) && length(selected) > top_n) {
    selected <- utils::head(selected, top_n)
  }

  ## Refit on selected features
  rf_refit <- NULL
  imp_refit <- NULL
  if (refit_on_selected) {
    set.seed(seed)
    message(sprintf("Phase 3: Re-fitting RF on %d selected features ...",
                    length(selected)))
    train_sel <- ml_data$train_x[, selected, drop = FALSE]
    rf_refit <- randomForest::randomForest(
      x = train_sel, y = ml_data$train_y,
      ntree = n_trees, importance = TRUE)
    imp_raw_r <- randomForest::importance(rf_refit)
    imp_refit <- data.frame(
      gene = rownames(imp_raw_r),
      MeanDecreaseAccuracy = imp_raw_r[, "MeanDecreaseAccuracy"],
      MeanDecreaseGini = imp_raw_r[, "MeanDecreaseGini"],
      stringsAsFactors = FALSE)
    imp_refit <- imp_refit[order(-imp_refit$MeanDecreaseGini), ]
    rownames(imp_refit) <- imp_refit$gene
  }

  ## The "active" model: refit if available, else full
  rf_active <- if (!is.null(rf_refit)) rf_refit else rf_full

  ## OOB performance from active model
  oob_active <- .compute_oob(rf_active, ml_data$train_y, pos, neg,
                             n_trees = n_trees)
  ## Stored as cv_performance for S3 compatibility
  cv_perf <- oob_active

  ## Test-set evaluation (Mode B / C)
  test_perf <- NULL
  if (!isTRUE(ml_data$full_cv) && !is.null(ml_data$test_x)) {
    .check_ml_deps("caret")
    test_x_sel <- if (!is.null(rf_refit))
      ml_data$test_x[, selected, drop = FALSE]
    else
      ml_data$test_x
    pred_class <- stats::predict(rf_active, newdata = test_x_sel)
    pred_prob  <- stats::predict(rf_active, newdata = test_x_sel, type = "prob")
    cm <- caret::confusionMatrix(pred_class, ml_data$test_y, positive = pos)
    auc_val <- NA_real_
    if (requireNamespace("pROC", quietly = TRUE) && pos %in% colnames(pred_prob)) {
      roc_obj <- pROC::roc(response = ml_data$test_y,
                           predictor = pred_prob[, pos],
                           levels = ml_data$levels, quiet = TRUE)
      auc_val <- as.numeric(pROC::auc(roc_obj))
    }
    test_perf <- list(
      predictions = pred_class,
      probabilities = as.data.frame(pred_prob),
      confusion = cm,
      accuracy = as.numeric(cm$overall["Accuracy"]),
      auc = auc_val,
      sensitivity = as.numeric(cm$byClass["Sensitivity"]),
      specificity = as.numeric(cm$byClass["Specificity"]))
  }

  ## Assemble tcm_ml object
  obj <- .new_tcm_ml(
    method = "rf",
    model = rf_active,
    importance = imp_df,
    selected_features = selected,
    cv_performance = cv_perf,
    test_performance = test_perf,
    ml_data = ml_data)

  obj$boruta <- boruta_res
  obj$boruta_fallback <- boruta_used_fallback
  obj$rf_full <- rf_full
  obj$rf_fit <- rf_active
  obj$rf_refit <- rf_refit
  obj$oob_full <- oob_full
  obj$importance_refit <- imp_refit
  obj$genes <- selected

  ## Summary message
  n_conf <- sum(boruta_res$finalDecision == "Confirmed")
  refit_tag <- if (refit_on_selected) " (refit)" else " (full)"
  message(sprintf(
    "  Boruta confirmed %d | Selected %d | OOB AUC%s = %.4f",
    n_conf, length(selected), refit_tag, cv_perf$auc))
  return(obj)
}


## Internal: compute OOB performance from a randomForest object
## @keywords internal
## @noRd
.compute_oob <- function(rf_fit, train_y, pos, neg, n_trees) {
  oob_err <- rf_fit$err.rate[n_trees, "OOB"]
  oob_auc <- NA_real_
  votes <- rf_fit$votes
  if (!is.null(votes) && requireNamespace("pROC", quietly = TRUE)) {
    if (pos %in% colnames(votes)) {
      roc_obj <- pROC::roc(response = train_y,
                           predictor = votes[, pos],
                           levels = levels(train_y), quiet = TRUE)
      oob_auc <- as.numeric(pROC::auc(roc_obj))
    }
  }
  cm_oob <- rf_fit$confusion[, seq_len(nlevels(train_y))]
  oob_sens <- cm_oob[pos, pos] / sum(cm_oob[pos, ])
  oob_spec <- cm_oob[neg, neg] / sum(cm_oob[neg, ])

  list(auc = oob_auc,
       auc_sd = NA_real_,
       accuracy = 1 - oob_err,
       oob_error = oob_err,
       sensitivity = oob_sens,
       specificity = oob_spec,
       perf_type = "OOB")
}


## SVM-RFE

#' SVM-RFE feature selection
#' Recursive Feature Elimination with SVM via [caret::rfe()].
#' Designed for binary classification of gene expression data (e.g. RNA-seq
#' profiles of PPI candidate genes vs. disease/control labels).
#' @param ml_data A \code{tcm_ml_data} object. Must contain exactly two classes.
#' @param sizes Integer vector of feature-subset sizes to evaluate in RFE.
#'   Default \code{NULL}: sizes are derived automatically from the data by
#'   log-spacing 10 values from 2 to \code{min(p - 1, max(p \%/\% 2, 50))},
#'   where \code{p} is the number of input features. When provided manually,
#'   values outside \code{[2, p]} are silently dropped.
#' @param top_n Override the auto-selected optimal size: keep the top
#'   \code{top_n} genes by aggregated importance. Default \code{NULL}
#'   (use the size that maximises \code{metric} in the RFE profile).
#' @param min_size Minimum number of features to retain. The profile search
#'   is restricted to \code{Variables >= min_size}. Default \code{5}.
#'   Useful to avoid biologically trivial 2-3 gene subsets.
#' @param kernel SVM kernel to use. \code{"svmLinear"} (default, recommended
#'   for gene expression data — weights \eqn{|w_j|^2} are directly used for
#'   recursive elimination per Guyon 2002) or \code{"svmRadial"}.
#' @param metric Optimization metric for RFE. Default \code{"ROC"}.
#'   Other options include \code{"Accuracy"}, \code{"Kappa"}.
#'   When \code{"ROC"} / \code{"Sens"} / \code{"Spec"} is used,
#'   \code{twoClassSummary} is applied automatically; otherwise
#'   \code{defaultSummary} is used.
#' @param cv_folds Outer CV folds. Default 5.
#' @param cv_repeats Outer CV repeats. Default 5.
#' @param seed Random seed. Default 2025.
#' @param ... Passed to [caret::rfe()].
#'
#' @return A \code{tcm_ml} object with \code{$rfe_result}, \code{$profile}, \code{$genes}.
#' @importFrom stats aggregate
#' @examples
#' \dontrun{
#'   ml_data <- prepare_ml_data(expr_mat, group)
#'   res <- ml_svm_rfe(ml_data)
#' }
#' @export
ml_svm_rfe <- function(ml_data,
                       sizes = NULL,
                       top_n = NULL,
                       min_size = 5,
                       kernel = c("svmLinear", "svmRadial"),
                       metric = "ROC",
                       cv_folds = 5,
                       cv_repeats = 5,
                       seed = 2025,
                       ...) {

  .check_ml_deps(c("caret", "kernlab"))
  stopifnot(inherits(ml_data, "tcm_ml_data"))
  if (nlevels(ml_data$train_y) != 2L)
    stop(sprintf("ml_svm_rfe() requires exactly 2 classes, got %d.",
                 nlevels(ml_data$train_y)), call. = FALSE)
  kernel <- match.arg(kernel)

  p <- ncol(ml_data$train_x)
  n <- nrow(ml_data$train_x)
  if (p > n)
    warning(sprintf(
      "p (%d) > n (%d): SVM-RFE may be unstable; consider pre-filtering candidate genes.",
      p, n), call. = FALSE)
  if (is.null(sizes)) {
    ## Auto grid: 10 log-spaced integers from 2 to min(p-1, max(p/2, 50))
    max_s <- min(p - 1L, max(as.integer(p %/% 2L), 50L))
    if (max_s <= 2L) {
      sizes <- 2L
    } else {
      sizes <- unique(as.integer(round(exp(seq(log(2), log(max_s), length.out = 10)))))
    }
  } else {
    sizes <- sort(unique(as.integer(sizes[sizes >= 2L & sizes <= p])))
    if (length(sizes) == 0L) sizes <- min(p, 10L)
  }

  rfe_funcs <- caret::caretFuncs
  use_class_probs <- metric %in% c("ROC", "Sens", "Spec")
  summary_func <- if (use_class_probs) caret::twoClassSummary else caret::defaultSummary
  rfe_funcs$summary <- summary_func

  ## Override fit so inner train() uses the correct metric.
  ## rfe() consumes `metric` as its own arg and does NOT forward it via `...`,

  ## so the default train() falls back to "Accuracy" — wrong for twoClassSummary.
  rfe_funcs$fit <- function(x, y, first, last, ...) {
    caret::train(x, y, metric = metric, ...)
  }

  rfe_ctrl <- caret::rfeControl(
    functions = rfe_funcs,
    method = "repeatedcv", number = cv_folds, repeats = cv_repeats,
    verbose = FALSE, allowParallel = TRUE,
    saveDetails = TRUE
  )

  message(sprintf("Running SVM-RFE (%s) ...", kernel))
  set.seed(seed)
  rfe_res <- caret::rfe(
    x = ml_data$train_x, y = ml_data$train_y,
    sizes = sizes, rfeControl = rfe_ctrl,
    method = kernel, metric = metric,
    preProcess = c("center", "scale"),
    ## Inner CV: single-fold CV for SVM hyperparameter tuning within each RFE fold.
    ## Intentionally plain trainControl (not .make_train_ctrl): outer rfeControl
    ## already handles repeatedcv, so inner CV only needs one pass.
    trControl = caret::trainControl(
      method = "cv", number = cv_folds,
      classProbs = use_class_probs, summaryFunction = summary_func
    ), ...
  )

  ## aggregate importance (used for ranking & gene selection)
  agg_imp <- stats::aggregate(Overall ~ var, data = rfe_res$variables, FUN = mean)
  agg_imp <- agg_imp[order(-agg_imp$Overall), ]
  names(agg_imp) <- c("gene", "importance")
  rownames(agg_imp) <- agg_imp$gene

  ## pick optimal size
  ## caret::rfe() always adds p to the profile; if it wins, the RFE is useless.
  ## We optionally allow users to override opt_size via top_n.
  profile <- rfe_res$results
  
  if (!is.null(top_n)) {
    opt_size <- min(as.integer(top_n), p)
    selected <- utils::head(agg_imp$gene, opt_size)
  } else {
    profile_sub <- profile[profile$Variables >= min_size & profile$Variables < p, , drop = FALSE]
    if (nrow(profile_sub) > 0 && metric %in% names(profile_sub)) {
      best_idx <- which.max(profile_sub[[metric]])
      opt_size <- profile_sub$Variables[best_idx]
      selected <- utils::head(agg_imp$gene, opt_size)
    } else {
      ## fallback: use caret's own choice
      opt_size <- rfe_res$optsize
      selected <- caret::predictors(rfe_res)
    }
  }

  opt_row <- profile[profile$Variables == opt_size, , drop = FALSE]
  metric_sd <- paste0(metric, "SD")
  metric_val <- if (nrow(opt_row) > 0 && metric %in% names(opt_row))
    as.numeric(opt_row[[metric]][1]) else NA_real_
  cv_perf <- list(
    auc = if (metric == "ROC") metric_val else NA_real_,
    auc_sd = if (metric == "ROC" && nrow(opt_row) > 0 && metric_sd %in% names(opt_row))
               as.numeric(opt_row[[metric_sd]][1]) else NA_real_,
    metric_value = metric_val,
    metric_sd = if (nrow(opt_row) > 0 && metric_sd %in% names(opt_row))
                  as.numeric(opt_row[[metric_sd]][1]) else NA_real_,
    sensitivity = if (nrow(opt_row) > 0 && "Sens" %in% names(opt_row))
                    as.numeric(opt_row$Sens[1]) else NA_real_,
    specificity = if (nrow(opt_row) > 0 && "Spec" %in% names(opt_row))
                    as.numeric(opt_row$Spec[1]) else NA_real_,
    metric_name = metric,
    opt_size = opt_size
  )

  test_perf <- NULL
  if (!isTRUE(ml_data$full_cv) && !is.null(ml_data$test_x)) {
    ctrl <- .make_train_ctrl(
      method = "repeatedcv",
      number = cv_folds,
      repeats = cv_repeats,
      seed = seed,
      class_probs = use_class_probs
    )
    set.seed(seed)
    fit_final <- caret::train(
      x = ml_data$train_x[, selected, drop = FALSE],
      y = ml_data$train_y, method = kernel,
      trControl = ctrl, metric = metric,
      preProcess = c("center", "scale")
    )
    test_perf <- .eval_test(fit_final,
                            ml_data$test_x[, selected, drop = FALSE],
                            ml_data$test_y, ml_data$levels[1])
  }

  obj <- .new_tcm_ml(method = "svm_rfe", model = rfe_res,
                     importance = agg_imp, selected_features = selected,
                     cv_performance = cv_perf, test_performance = test_perf,
                     ml_data = ml_data)
  obj$rfe_result <- rfe_res
  obj$profile <- rfe_res$results
  obj$genes <- selected

  ## Extract OOF probabilities from saved CV predictions
  ## rfe_res$pred contains per-fold predictions at each subset size;
  ## we filter for the optimal size and average across repeats.
  oof_pred <- rfe_res$pred
  if (!is.null(oof_pred) && nrow(oof_pred) > 0 &&
      "Variables" %in% names(oof_pred)) {
    pos <- ml_data$levels[1]
    sub_pred <- oof_pred[oof_pred$Variables == opt_size, , drop = FALSE]
    if (nrow(sub_pred) > 0 && pos %in% names(sub_pred)) {
      ## Average probability across repeated CV folds per sample
      avg <- stats::aggregate(
        sub_pred[[pos]], by = list(rowIndex = sub_pred$rowIndex), FUN = mean
      )
      avg <- avg[order(avg$rowIndex), ]
      obj$oof_prob <- avg$x
    }
  }

  message(sprintf("  Optimal: %d features (excl. full set) | CV %s = %s",
                  opt_size, metric,
                  if (is.na(cv_perf$metric_value)) "NA" else sprintf("%.4f", cv_perf$metric_value)))
  return(obj)
}


## XGBoost

#' XGBoost feature selection
#' Trains an XGBoost classifier with \code{nrounds}-round CV, then ranks features by information gain.
#' Since XGBoost does not inherently perform variable selection, \code{top_n} controls how many top features to keep.
#' When \code{top_n = NULL} (default), all features with Gain > 0 are retained; use [select_features()] afterwards to trim.
#' @param ml_data A \code{tcm_ml_data} object.
#' @param top_n Number of top features to select by Gain. Default \code{NULL} (keep all).
#' @param eval_metric Evaluation metric for XGBoost CV. Default \code{"auc"}.
#'   Other options include \code{"error"}, \code{"logloss"}, \code{"aucpr"}.
#'   See \code{xgboost} documentation for full list.
#' @param nrounds Maximum boosting rounds. Default 200.
#' @param max_depth Maximum tree depth. Default 6.
#' @param eta Learning rate. Default 0.1.
#' @param cv_folds CV folds for early stopping. Default 5.
#' @param early_stopping_rounds Stop if no improvement for this many rounds. Default 20.
#' @param seed Random seed. Default 2025.
#' @param ... Passed to [xgboost::xgb.cv()] and [xgboost::xgb.train()].
#'
#' @return A \code{tcm_ml} object with \code{$xgb_fit}, \code{$cv_log}, \code{$genes}.
#' @importFrom stats predict
#' @examples
#' \dontrun{
#'   ml_data <- prepare_ml_data(expr_mat, group)
#'   res <- ml_xgboost(ml_data)
#' }
#' @export
ml_xgboost <- function(ml_data,
                       top_n = NULL,
                       eval_metric = "auc",
                       nrounds = 200,
                       max_depth = 6,
                       eta = 0.1,
                       cv_folds = 5,
                       early_stopping_rounds = 20,
                       seed = 2025,
                       ...) {

  .check_ml_deps("xgboost")
  stopifnot(inherits(ml_data, "tcm_ml_data"))
  if (nlevels(ml_data$train_y) != 2L) {
    stop("ml_xgboost() requires a binary outcome (exactly 2 levels).",
         call. = FALSE)
  }

  x_mat <- as.matrix(ml_data$train_x)
  y_num <- as.integer(ml_data$train_y == ml_data$levels[1])

  dtrain <- xgboost::xgb.DMatrix(data = x_mat, label = y_num)

  params <- list(
    objective = "binary:logistic",
    eval_metric = eval_metric,
    max_depth = max_depth,
    eta = eta,
    ...
  )

  ## CV for optimal nrounds
  set.seed(seed)
  message("Running XGBoost CV ...")
  cv_res <- xgboost::xgb.cv(
    params = params,
    data = dtrain,
    nrounds = nrounds,
    nfold = cv_folds,
    early_stopping_rounds = early_stopping_rounds,
    prediction = TRUE,
    verbose = 0
  )

  ## xgboost >= 2.1: best_iteration moved into $early_stop sub-list
  best_iter <- if (!is.null(cv_res$best_iteration)) {
    cv_res$best_iteration
  } else if (!is.null(cv_res$early_stop$best_iteration)) {
    cv_res$early_stop$best_iteration
  } else {
    NULL
  }
  if (is.null(best_iter) || length(best_iter) == 0L) {
    ## fallback: locate peak from evaluation_log
    cv_log_tmp <- as.data.frame(cv_res$evaluation_log)
    metric_col <- grep(paste0("test.*", eval_metric, ".*mean"),
                       names(cv_log_tmp), value = TRUE)[1]
    best_iter <- if (!is.na(metric_col)) which.max(cv_log_tmp[[metric_col]])
                 else nrounds
  }
  ## Final safety guard — which.max() may return integer(0) on all-NaN columns
  if (length(best_iter) == 0L || !is.finite(best_iter) || best_iter <= 0L)
    best_iter <- nrounds
  cv_log <- as.data.frame(cv_res$evaluation_log)

  ## Train final model
  set.seed(seed)
  message(sprintf("  Best iteration: %d | Training final model ...", best_iter))
  xgb_fit <- xgboost::xgb.train(
    params = params,
    data = dtrain,
    nrounds = best_iter,
    verbose = 0
  )

  ## Importance
  imp_mat <- xgboost::xgb.importance(model = xgb_fit)
  ## xgboost >= 2.1 renamed 'Feature' → 'Features'
  feat_col <- if ("Feature" %in% names(imp_mat)) "Feature" else "Features"
  imp_df <- data.frame(
    gene = imp_mat[[feat_col]],
    importance = imp_mat$Gain,
    stringsAsFactors = FALSE
  )
  imp_df <- imp_df[order(-imp_df$importance), ]
  rownames(imp_df) <- imp_df$gene

  if (!is.null(top_n)) {
    n_keep <- min(as.integer(top_n), nrow(imp_df))
    selected <- imp_df$gene[seq_len(n_keep)]
  } else {
    selected <- imp_df$gene
  }

  ## CV metric from cv_log
  auc_col <- grep(paste0("test.*", eval_metric, ".*mean"),
                  names(cv_log), value = TRUE)[1]
  cv_auc <- if (!is.na(auc_col)) cv_log[[auc_col]][best_iter] else NA_real_
  sd_col <- grep(paste0("test.*", eval_metric, ".*std"),
                 names(cv_log), value = TRUE)[1]
  cv_sd <- if (!is.na(sd_col)) cv_log[[sd_col]][best_iter] else NA_real_

  ## CV Sens / Spec from out-of-fold predictions
  ## xgboost >= 2.1: predictions moved into $cv_predict sub-list
  oof_prob <- if (!is.null(cv_res$pred)) {
    cv_res$pred
  } else if (!is.null(cv_res$cv_predict$pred)) {
    cv_res$cv_predict$pred
  } else {
    NULL
  }
  if (is.null(oof_prob)) {
    warning("Out-of-fold predictions unavailable; CV Sens/Spec set to NA.",
            call. = FALSE)
  }
  oof_cls <- if (!is.null(oof_prob)) {
    factor(
      ifelse(oof_prob > 0.5, ml_data$levels[1], ml_data$levels[2]),
      levels = ml_data$levels
    )
  }
  pos <- ml_data$levels[1]; neg <- ml_data$levels[2]
  if (!is.null(oof_cls)) {
    tp <- sum(oof_cls == pos & ml_data$train_y == pos)
    fn <- sum(oof_cls == neg & ml_data$train_y == pos)
    tn <- sum(oof_cls == neg & ml_data$train_y == neg)
    fp <- sum(oof_cls == pos & ml_data$train_y == neg)
    cv_sens <- tp / (tp + fn)
    cv_spec <- tn / (tn + fp)
  } else {
    cv_sens <- NA_real_
    cv_spec <- NA_real_
  }

  cv_perf <- list(auc = cv_auc, auc_sd = cv_sd,
                  sensitivity = cv_sens,
                  specificity = cv_spec,
                  best_iteration = best_iter)

  ## Test evaluation
  test_perf <- NULL
  if (!isTRUE(ml_data$full_cv) && !is.null(ml_data$test_x)) {
    .check_ml_deps("caret")
    test_mat <- as.matrix(ml_data$test_x)
    pred_prob <- stats::predict(xgb_fit, newdata = test_mat)
    pred_cls <- factor(ifelse(pred_prob > 0.5,
                               ml_data$levels[1], ml_data$levels[2]),
                        levels = ml_data$levels)
    cm <- caret::confusionMatrix(pred_cls, ml_data$test_y,
                                 positive = ml_data$levels[1])
    auc_val <- NA_real_
    if (requireNamespace("pROC", quietly = TRUE)) {
      roc_obj <- pROC::roc(response = ml_data$test_y,
                           predictor = pred_prob,
                           levels = ml_data$levels, quiet = TRUE)
      auc_val <- as.numeric(pROC::auc(roc_obj))
    }
    test_perf <- list(predictions = pred_cls,
                      probabilities = stats::setNames(
                        data.frame(pred_prob, 1 - pred_prob),
                        ml_data$levels
                      ),
                      confusion = cm,
                      accuracy = as.numeric(cm$overall["Accuracy"]),
                      auc = auc_val,
                      sensitivity = as.numeric(cm$byClass["Sensitivity"]),
                      specificity = as.numeric(cm$byClass["Specificity"]))
  }

  obj <- .new_tcm_ml(method = "xgboost", model = xgb_fit,
                     importance = imp_df, selected_features = selected,
                     cv_performance = cv_perf, test_performance = test_perf,
                     ml_data = ml_data)
  obj$xgb_fit <- xgb_fit
  obj$cv_log <- cv_log
  obj$genes <- selected
  ## Store OOF probability as P(levels[1]) for ROC plotting
  obj$oof_prob <- as.numeric(oof_prob)

  message(sprintf("  Top %d / %d features | CV %s = %s",
                  length(selected), ncol(x_mat), eval_metric,
                  if (is.na(cv_auc)) "NA" else sprintf("%.4f", cv_auc)))
  return(obj)
}


#' Run multiple ML methods for feature selection
#' Convenience wrapper that runs LASSO / Elastic Net / Ridge / RF / SVM-RFE / XGBoost.
#' By default \code{top_n = NULL}: Ridge and XGBoost retain all features, SVM-RFE auto-selects the optimal size.
#' After inspecting results, call [select_features()] on individual models to trim to the desired count.
#' @param ml_data A \code{tcm_ml_data} object.
#' @param methods Subset of \code{c("lasso", "enet", "ridge", "rf", "svm_rfe", "xgboost")}.
#' @param seed Random seed. Default 2025.
#' @param cv_folds CV folds. Default 5.
#' @param cv_repeats SVM-RFE outer CV repeats. Default 5.
#' @param alpha Elastic Net alpha. Default 0.5 (ignored for lasso/ridge).
#' @param top_n Top-n features for Ridge / XGBoost / SVM-RFE. Default \code{NULL} (keep all or auto-select).
#' @param ... Forwarded to individual model functions.
#'
#' @return A named \code{tcm_ml_list}.
#' @examples
#' \dontrun{
#'   ml_data <- prepare_ml_data(expr_mat, group)
#'   res_list <- run_ml_screening(ml_data, methods = c("lasso", "rf"))
#'   summary(res_list)
#' }
#' @export
run_ml_screening <- function(ml_data,
                             methods = c("lasso", "rf", "svm_rfe", "xgboost"),
                             seed = 2025,
                             cv_folds = 5,
                             cv_repeats = 5,
                             alpha = 0.5,
                             top_n = NULL,
                             ...) {

  stopifnot(inherits(ml_data, "tcm_ml_data"))
  methods <- match.arg(methods,
                       c("lasso", "enet", "ridge", "rf", "svm_rfe", "xgboost"),
                       several.ok = TRUE)

  results <- list()
  for (i in seq_along(methods)) {
    m <- methods[i]
    message(sprintf("\n[%d/%d] %s", i, length(methods), toupper(m)))
    results[[m]] <- switch(m,
      lasso = ml_lasso(ml_data, cv_folds = cv_folds, seed = seed, ...),
      enet = ml_enet(ml_data, alpha = alpha, cv_folds = cv_folds,
                     seed = seed, ...),
      ridge = ml_ridge(ml_data, top_n = top_n, cv_folds = cv_folds,
                       seed = seed, ...),
      rf = ml_rf(ml_data, seed = seed, ...),
      svm_rfe = ml_svm_rfe(ml_data, top_n = top_n, cv_folds = cv_folds,
                            cv_repeats = cv_repeats, seed = seed, ...),
      xgboost = ml_xgboost(ml_data, top_n = top_n, cv_folds = cv_folds,
                            seed = seed, ...)
    )
  }

  class(results) <- "tcm_ml_list"
  attr(results, "ml_data") <- ml_data
  message("\nAll methods complete.")
  return(results)
}


## Consensus
#' Consensus features across ML methods
#' Genes selected by >= \code{min_methods} models.
#'
#' @param ml_list A \code{tcm_ml_list} or plain list of \code{tcm_ml} objects.
#' @param min_methods Minimum agreement count. Default 2.
#'
#' @return Character vector of consensus gene symbols.
#' @importFrom utils combn head
#' @examples
#' \dontrun{
#'   ml_data <- prepare_ml_data(expr_mat, group)
#'   res_list <- run_ml_screening(ml_data)
#'   get_ml_consensus(res_list, min_methods = 2)
#' }
#' @export
get_ml_consensus <- function(ml_list, min_methods = 2) {

  feat_list <- lapply(ml_list, function(m) m$selected_features)
  all_genes <- table(unlist(feat_list))
  consensus <- names(sort(all_genes[all_genes >= min_methods], decreasing = TRUE))

  message(sprintf("Consensus (>= %d / %d methods): %d genes",
                  min_methods, length(ml_list), length(consensus)))

  if (length(consensus) == 0 && length(feat_list) >= 2) {
    pairs <- utils::combn(names(feat_list), 2, simplify = FALSE)
    for (pr in pairs) {
      shared <- intersect(feat_list[[pr[1]]], feat_list[[pr[2]]])
      if (length(shared) > 0)
        message(sprintf("  %s & %s share %d: %s",
                        toupper(pr[1]), toupper(pr[2]), length(shared),
                        paste(utils::head(shared, 5), collapse = ", ")))
    }
  }

  return(consensus)
}


#' Assemble individually fitted ML models into a \code{tcm_ml_list}
#' @param ... Named or unnamed \code{tcm_ml} objects. Names become the method labels used in downstream plots. At least one object must be supplied.
#' @return A \code{tcm_ml_list} object identical in structure to the output of \code{\link{run_ml_screening}}, with the \code{"ml_data"} attribute set from the first element.
#'
#' @examples
#' \dontrun{
#'   ml_data <- prepare_ml_data(expr_mat, group)
#'   res_lasso <- ml_lasso(ml_data)
#'   res_rf <- ml_rf(ml_data)
#'   res_svm <- ml_svm_rfe(ml_data)
#'
#'   ## assemble into a tcm_ml_list
#'   ml_list <- create_tcm_ml_list(
#'     lasso   = res_lasso,
#'     rf      = res_rf,
#'     svm_rfe = res_svm)
#' }
#' @export
create_tcm_ml_list <- function(...) {
  lst <- list(...)

  if (length(lst) == 0L)
    stop("Please provide at least one tcm_ml object.", call. = FALSE)

  bad <- vapply(lst, function(x) !inherits(x, "tcm_ml"), logical(1L))
  if (any(bad))
    stop(sprintf(
      "All inputs must be <tcm_ml> objects. Non-conforming element(s) at position(s): %s",
      paste(which(bad), collapse = ", ")
    ), call. = FALSE)

  nms <- names(lst)
  if (is.null(nms)) nms <- rep("", length(lst))
  unnamed_idx <- which(nms == "")
  if (length(unnamed_idx) > 0L) {
    used <- nms[nms != ""]
    for (i in unnamed_idx) {
      base_nm <- if (!is.null(lst[[i]]$method)) lst[[i]]$method else paste0("model", i)
      nm <- base_nm
      suffix <- 1L
      while (nm %in% used) {
        nm <- paste0(base_nm, "_", suffix)
        suffix <- suffix + 1L
      }
      nms[i] <- nm
      used <- c(used, nm)
    }
    names(lst) <- nms
  }

  class(lst) <- "tcm_ml_list"

  if (!is.null(lst[[1L]]$ml_data))
    attr(lst, "ml_data") <- lst[[1L]]$ml_data

  return(lst)
}