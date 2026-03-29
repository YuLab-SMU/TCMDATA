# S3 print and as.character methods for tcm_ai_analysis and tcm_ai_draft.

#' @export
print.tcm_ai_analysis <- function(x, ...) {
  meta <- x$metadata
  cat(sprintf("=== TCM AI Analysis: %s ===\n", x$input$type))
  cat(sprintf("  Model: %s | Language: %s | Audience: %s\n",
              meta$model, meta$language, meta$audience))
  cat(sprintf("  Generated: %s\n\n", format(meta$generated_at)))

  out <- x$output
  if (is.null(out)) {
    cat("  (No output — generation may have failed)\n")
    return(invisible(x))
  }

  # Fallback: model returned plain text instead of structured JSON
  if (identical(meta$output_mode %||% "structured", "fallback_text")) {
    cat("-- Overall Interpretation --\n")
    cat(out$summary, "\n\n")
    cat("  [Fallback text output; structured fields unavailable.]\n\n")
    return(invisible(x))
  }

  cat("-- Summary --\n")
  cat(out$summary, "\n\n")

  if (length(out$key_findings) > 0) {
    cat("-- Key Findings --\n")
    for (f in out$key_findings) cat("  * ", f, "\n")
    cat("\n")
  }

  cat("-- Biological Interpretation --\n")
  cat(out$biological_interpretation, "\n\n")

  if (!is.null(out$tcm_relevance) && nzchar(out$tcm_relevance)) {
    cat("-- TCM Relevance --\n")
    cat(out$tcm_relevance, "\n\n")
  }

  if (length(out$caveats) > 0) {
    cat("-- Caveats --\n")
    for (c in out$caveats) cat("  ! ", c, "\n")
    cat("\n")
  }

  invisible(x)
}

#' @export
print.tcm_ai_draft <- function(x, ...) {
  meta <- x$metadata
  cat(sprintf("=== TCM AI Draft: %s ===\n", x$input$type))
  cat(sprintf("  Model: %s | Language: %s\n",
              meta$model, meta$language))
  cat(sprintf("  Generated: %s\n\n", format(meta$generated_at)))

  draft <- x$draft
  if (is.null(draft)) {
    cat("  (No draft — generation may have failed)\n")
    return(invisible(x))
  }

  # Fallback: model returned plain text instead of structured JSON
  if (identical(meta$output_mode %||% "structured", "fallback_text")) {
    cat("-- Result Paragraph --\n")
    cat(draft$paragraph, "\n\n")
    cat("  [Fallback text output; structured fields unavailable.]\n\n")
    return(invisible(x))
  }

  cat("-- Result Paragraph --\n")
  cat(draft$paragraph, "\n\n")

  if (!is.null(draft$figure_legend_hint) &&
      nzchar(draft$figure_legend_hint)) {
    cat("-- Figure Legend Hint --\n")
    cat(draft$figure_legend_hint, "\n\n")
  }

  invisible(x)
}

#' @export
as.character.tcm_ai_analysis <- function(x, ...) {
  x$output$summary %||% ""
}

#' @export
as.character.tcm_ai_draft <- function(x, ...) {
  x$draft$paragraph %||% ""
}

#' @export
print.tcm_ai_custom <- function(x, ...) {
  meta <- x$metadata
  cat(sprintf("=== TCM AI Custom: %s ===\n", x$input$type))
  cat(sprintf("  Model: %s | Language: %s | Audience: %s\n",
              meta$model, meta$language, meta$audience))
  cat(sprintf("  Generated: %s\n\n", format(meta$generated_at)))

  out <- x$output
  if (is.null(out) || length(out) == 0L) {
    cat("  (No output — generation may have failed)\n")
    return(invisible(x))
  }

  # Fallback: model returned plain text instead of structured JSON
  if (identical(meta$output_mode %||% "structured", "fallback_text")) {
    cat("-- Overall Interpretation --\n")
    cat(out$raw_text %||% "", "\n\n")
    cat("  [Fallback text output; structured fields unavailable.]\n\n")
    return(invisible(x))
  }

  cat(sprintf("Fields: %s\n\n", paste(names(out), collapse = ", ")))

  for (nm in names(out)) {
    val <- out[[nm]]
    cat(sprintf("-- %s --\n", nm))
    if (is.list(val)) {
      for (item in val) cat("  *", as.character(item), "\n")
    } else {
      cat(as.character(val), "\n")
    }
    cat("\n")
  }

  invisible(x)
}

#' @export
as.character.tcm_ai_custom <- function(x, ...) {
  out <- x$output
  if (is.null(out)) return("")
  # Return $summary if present, otherwise the first non-empty string field
  if (!is.null(out$summary) && nzchar(out$summary)) return(out$summary)
  for (nm in names(out)) {
    val <- out[[nm]]
    if (is.character(val) && length(val) == 1L && nzchar(val)) return(val)
  }
  ""
}
