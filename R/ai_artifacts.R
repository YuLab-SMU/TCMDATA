# ai_artifacts.R
# Artifact store for TCM agent layer.
# Manages complex R objects as addressable, reusable, trackable handles.

# Internal artifact registry (session-level environment)
.tcm_artifact_env <- new.env(parent = emptyenv())
.tcm_artifact_meta_env <- new.env(parent = emptyenv())

#' Save a TCM artifact to the registry
#'
#' Stores a complex R object (igraph, enrichResult, data.frame, etc.) and
#' returns a lightweight handle that can be passed to the LLM.
#'
#' @param object The R object to store.
#' @param artifact_type Character. Type label, e.g. "enrichment_result",
#'   "herb_graph", "ppi_graph", "ranking_table".
#' @param artifact_id Character or NULL. If NULL, auto-generated.
#' @param summary Character or NULL. Human-readable summary for LLM context.
#' @param provenance List or NULL. Metadata about how this artifact was created.
#'
#' @return A list with artifact metadata (artifact_id, artifact_type, summary,
#'   preview, provenance, created_at). The actual object is stored internally.
#'
#' @examples
#' \dontrun{
#'   # Save an enrichment result
#'   enrich_res <- herb_enricher(genes = my_genes)
#'   artifact <- save_tcm_artifact(
#'     object = enrich_res,
#'     artifact_type = "enrichment_result",
#'     summary = "Herb enrichment on 126 genes; 14 terms passed cutoff."
#'   )
#'   # artifact$artifact_id can be passed to next tool call
#' }
#' @export
save_tcm_artifact <- function(object,
                              artifact_type,
                              artifact_id = NULL,
                              summary = NULL,
                              provenance = NULL) {
  # Auto-generate ID if not provided

  if (is.null(artifact_id)) {
    prefix <- switch(
      artifact_type,
      "enrichment_result" = "enrich",
      "herb_graph"        = "hgraph",
      "ppi_graph"         = "ppi",
      "ranking_table"     = "rank",
      "mcode_result"      = "mcode",
      "ml_data"           = "ml",
      "ml_result"         = "ml",
      "pubmed_evidence"   = "pubmed",
      "pubchem_result"    = "chem",
      "table_result"      = "table",
      "plot"              = "plot",
      "search_result"     = "search",
      "art"
    )
    existing_ids <- unique(c(ls(.tcm_artifact_env), ls(.tcm_artifact_meta_env)))
    pattern <- paste0("^", prefix, "_")
    existing_nums <- as.integer(
      gsub(pattern, "", grep(pattern, existing_ids, value = TRUE))
    )
    next_num <- if (length(existing_nums) == 0) 1L else max(existing_nums) + 1L
    artifact_id <- sprintf("%s_%03d", prefix, next_num)
  }

  # Generate preview based on object type
  preview <- .generate_preview(object)

  # Generate summary if not provided
  if (is.null(summary)) {
    summary <- .generate_summary(object, artifact_type)
  }

  # Build provenance if not provided
  if (is.null(provenance)) {
    provenance <- list(created_at = Sys.time())
  } else {
    provenance$created_at <- provenance$created_at %||% Sys.time()
  }
  provenance$r_class <- class(object)[1]

  metadata <- list(
    artifact_id = artifact_id,
    artifact_type = artifact_type,
    summary = summary,
    preview = preview,
    provenance = provenance,
    created_at = provenance$created_at,
    r_class = class(object)[1]
  )

  # Store the actual object
  assign(artifact_id, object, envir = .tcm_artifact_env)
  assign(artifact_id, metadata, envir = .tcm_artifact_meta_env)

  # Return handle (not the object itself)
  structure(
    metadata,
    class = c("tcm_artifact_handle", "list")
  )
}

#' Load a TCM artifact from the registry
#'
#' Retrieves the actual R object by its artifact_id.
#'
#' @param artifact_id Character. The artifact identifier.
#'
#' @return The stored R object.
#'
#' @examples
#' \dontrun{
#'   enrich_obj <- load_tcm_artifact("enrich_001")
#' }
#' @export
load_tcm_artifact <- function(artifact_id) {
  if (!exists(artifact_id, envir = .tcm_artifact_env)) {
    stop(sprintf("Artifact '%s' not found.", artifact_id), call. = FALSE)
  }
  get(artifact_id, envir = .tcm_artifact_env)
}

#' List all TCM artifacts in the registry
#'
#' Returns a data.frame summarizing all stored artifacts.
#'
#' @return A data.frame with columns: artifact_id, artifact_type, r_class,
#'   created_at.
#'
#' @examples
#' \dontrun{
#'   list_tcm_artifacts()
#' }
#' @export
list_tcm_artifacts <- function() {
  ids <- ls(.tcm_artifact_meta_env)
  if (length(ids) == 0) {
    return(data.frame(
      artifact_id   = character(0),
      artifact_type = character(0),
      r_class       = character(0),
      created_at    = character(0),
      stringsAsFactors = FALSE
    ))
  }

  # Retrieve metadata from stored objects
  info <- lapply(ids, function(id) {
    meta <- get(id, envir = .tcm_artifact_meta_env)
    list(
      artifact_id   = id,
      artifact_type = meta$artifact_type %||% "unknown",
      r_class       = meta$r_class %||% NA_character_,
      created_at    = as.character(meta$created_at %||% NA)
    )
  })

  do.call(rbind, lapply(info, as.data.frame, stringsAsFactors = FALSE))
}

#' Summarize a TCM artifact for LLM context
#'
#' Generates a concise text summary suitable for passing to the LLM.
#'
#' @param artifact_id Character. The artifact identifier.
#' @param max_lines Integer. Maximum lines for preview tables.
#'
#' @return Character string with summary.
#'
#' @examples
#' \dontrun{
#'   summarize_tcm_artifact("enrich_001")
#' }
#' @export
summarize_tcm_artifact <- function(artifact_id, max_lines = 5L) {
  obj <- load_tcm_artifact(artifact_id)
  meta <- .get_artifact_metadata(artifact_id)
  summary <- meta$summary %||% .generate_summary(obj, meta$artifact_type %||% "unknown")
  preview <- .generate_preview(obj, max_rows = max_lines)

  paste0(
    "Artifact: ", artifact_id, "\n",
    "Summary: ", summary, "\n",
    "Preview:\n", preview

)
}

#' Clear all artifacts from registry
#'
#' Removes all stored artifacts. Use with caution.
#'
#' @return Invisible NULL.
#' @export
clear_tcm_artifacts <- function() {
  rm(list = ls(.tcm_artifact_env), envir = .tcm_artifact_env)
  rm(list = ls(.tcm_artifact_meta_env), envir = .tcm_artifact_meta_env)
  invisible(NULL)
}

#' Check if artifact exists
#'
#' @param artifact_id Character. The artifact identifier.
#' @return Logical.
#' @export
artifact_exists <- function(artifact_id) {
  exists(artifact_id, envir = .tcm_artifact_env) ||
    exists(artifact_id, envir = .tcm_artifact_meta_env)
}

#' Get artifact metadata from the registry
#' @keywords internal
#' @noRd
.get_artifact_metadata <- function(artifact_id) {
  if (!exists(artifact_id, envir = .tcm_artifact_meta_env)) {
    return(NULL)
  }
  get(artifact_id, envir = .tcm_artifact_meta_env)
}

#' Build a compact preview table for ML screening results
#' @keywords internal
#' @noRd
.preview_ml_list <- function(object, max_rows = 5L) {
  methods <- names(object)
  if (length(methods) == 0) {
    return("(no machine learning methods found)")
  }

  rows <- lapply(methods, function(method_name) {
    model <- object[[method_name]]
    perf <- model$test_performance %||% model$cv_performance %||% list()
    data.frame(
      method = method_name,
      features = length(model$selected_features %||% character(0)),
      auc = perf$auc %||% NA_real_,
      sensitivity = perf$sensitivity %||% NA_real_,
      specificity = perf$specificity %||% NA_real_,
      stringsAsFactors = FALSE
    )
  })

  df <- utils::head(do.call(rbind, rows), max_rows)
  paste(utils::capture.output(print(df)), collapse = "\n")
}

#' Generate preview string from object
#' @keywords internal
#' @noRd
.generate_preview <- function(object, max_rows = 5L) {
  if (inherits(object, "tcm_pubmed")) {
    df <- object$data
    if (is.null(df) || nrow(df) == 0) {
      return("(no PubMed records)")
    }
    cols <- intersect(c("PMID", "Title", "Year", "Journal"), names(df))
    df <- utils::head(df[, cols, drop = FALSE], max_rows)
    paste(utils::capture.output(print(df)), collapse = "\n")
  } else if (inherits(object, "tcm_ml_list")) {
    .preview_ml_list(object, max_rows = max_rows)
  } else if (inherits(object, "tcm_ml")) {
    imp <- object$importance
    if (is.null(imp) || !is.data.frame(imp) || nrow(imp) == 0) {
      return("(no feature importance table)")
    }
    df <- utils::head(imp, max_rows)
    paste(utils::capture.output(print(df)), collapse = "\n")
  } else if (inherits(object, "enrichResult")) {
    df <- utils::head(as.data.frame(object), max_rows)
    if (nrow(df) == 0) return("(no significant terms)")
    cols <- intersect(c("ID", "Description", "pvalue", "p.adjust", "Count"), names(df))
    if (length(cols) > 0) df <- df[, cols, drop = FALSE]
    paste(utils::capture.output(print(df)), collapse = "\n")
  } else if (inherits(object, "ggplot")) {
    "ggplot object"
  } else if (any(class(object) %in% c("Heatmap", "HeatmapList"))) {
    paste(class(object), collapse = ", ")
  } else if (inherits(object, "igraph")) {
    sprintf("igraph: %d nodes, %d edges",
            igraph::vcount(object), igraph::ecount(object))
  } else if (is.data.frame(object)) {
    df <- utils::head(object, max_rows)
    paste(utils::capture.output(print(df)), collapse = "\n")
  } else if (is.list(object)) {
    sprintf("list with %d elements: %s",
            length(object),
            paste(utils::head(names(object), 5), collapse = ", "))
  } else {
    paste(utils::capture.output(utils::str(object, max.level = 1)), collapse = "\n")
  }
}

#' Generate summary string from object
#' @keywords internal
#' @noRd
.generate_summary <- function(object, artifact_type) {
  if (inherits(object, "tcm_pubmed")) {
    n_articles <- nrow(object$data %||% data.frame())
    years <- object$stats$Year %||% numeric(0)
    if (length(years) > 0) {
      sprintf("PubMed evidence: %d articles across %d-%d.",
              n_articles, min(years), max(years))
    } else {
      sprintf("PubMed evidence: %d article(s).", n_articles)
    }
  } else if (inherits(object, "tcm_ml_data")) {
    train_n <- nrow(object$train_x %||% data.frame())
    test_n <- nrow(object$test_x %||% data.frame())
    gene_n <- length(object$gene_names %||% character(0))
    if (isTRUE(object$full_cv) || is.null(object$test_x)) {
      sprintf("ML dataset: %d samples and %d genes for cross-validation.",
              train_n, gene_n)
    } else {
      sprintf("ML dataset: %d training samples, %d test samples, and %d genes.",
              train_n, test_n, gene_n)
    }
  } else if (inherits(object, "tcm_ml_list")) {
    best_auc <- vapply(object, function(model) {
      perf <- model$test_performance %||% model$cv_performance %||% list(auc = NA_real_)
      as.numeric(perf$auc %||% NA_real_)
    }, numeric(1))
    best_name <- names(best_auc)[which.max(best_auc)]
    if (length(best_name) == 0 || is.na(best_auc[best_name])) {
      sprintf("ML screening result: %d method(s).", length(object))
    } else {
      sprintf("ML screening result: %d method(s); best AUC %.3f from %s.",
              length(object), best_auc[best_name], toupper(best_name))
    }
  } else if (inherits(object, "tcm_ml")) {
    perf <- object$test_performance %||% object$cv_performance %||% list(auc = NA_real_)
    sprintf("ML model %s: %d selected feature(s), AUC %.3f.",
            toupper(object$method %||% "unknown"),
            length(object$selected_features %||% character(0)),
            as.numeric(perf$auc %||% NA_real_))
  } else if (inherits(object, "enrichResult")) {
    df <- as.data.frame(object)
    sig <- sum(df$p.adjust < 0.05, na.rm = TRUE)
    sprintf("Enrichment result: %d terms, %d significant (p.adjust < 0.05)",
            nrow(df), sig)
  } else if (inherits(object, "ggplot") || identical(artifact_type, "plot")) {
    "Plot object ready for rendering."
  } else if (any(class(object) %in% c("Heatmap", "HeatmapList"))) {
    "ComplexHeatmap object ready for rendering."
  } else if (inherits(object, "igraph")) {
    sprintf("Graph: %d nodes, %d edges",
            igraph::vcount(object), igraph::ecount(object))
  } else if (is.data.frame(object)) {
    sprintf("Data frame: %d rows, %d columns", nrow(object), ncol(object))
  } else {
    sprintf("Object of class %s", class(object)[1])
  }
}

#' Print method for artifact handle
#' @param x A tcm_artifact_handle object.
#' @param ... Additional arguments (ignored).
#' @export
print.tcm_artifact_handle <- function(x, ...) {
  cat("TCM Artifact Handle\n")
  cat("  ID:      ", x$artifact_id, "\n")
  cat("  Type:    ", x$artifact_type, "\n")
  cat("  Summary: ", x$summary, "\n")
  cat("  Created: ", as.character(x$created_at), "\n")
  if (!is.null(x$preview)) {
    cat("  Preview:\n")
    cat("    ", gsub("\n", "\n    ", x$preview), "\n")
  }
  invisible(x)
}
