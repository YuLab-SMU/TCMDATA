# draft_result_paragraph(): publication-ready result paragraph
# from a tcm_ai_analysis object or from raw analysis output.

#' Draft a publication-ready result paragraph
#'
#' Accepts either a \code{tcm_ai_analysis} object (from any
#' \code{interpret_*} function) or a raw analysis object. When a
#' raw object is supplied, the context is compressed on the fly.
#'
#' @param x A \code{tcm_ai_analysis} object, or a raw analysis
#'   object (enrichResult, igraph, data.frame).
#' @param type Character. Required when \code{x} is a raw object.
#'   One of \code{"enrichment"}, \code{"ppi"}, \code{"table"}.
#' @param language Character. \code{"zh"} or \code{"en"}.
#'   Default \code{"zh"}.
#' @param prompt Character or NULL. Custom drafting instruction. Replaces the
#'   default \emph{"Draft a result paragraph from this \{type\} data:"}
#'   opening; the compressed data context is always appended regardless.
#'   Use this to steer tone, focus, or word count, e.g.
#'   \code{"用不超过3句话总结以下富集分析，聚焦免疫相关通路："}.
#' @param model Model identifier or NULL.
#'
#' @return A \code{tcm_ai_draft} object.
#' @examples
#' \dontrun{
#'   ai_res <- interpret_enrichment(enrich_res)
#'   draft  <- draft_result_paragraph(ai_res)
#'   cat(draft$draft$paragraph)
#'
#'   # Custom instruction
#'   draft2 <- draft_result_paragraph(ai_res,
#'     prompt = "用不超过3句话总结，聚焦免疫相关通路，语气适合Nature子刊：")
#' }
#' @export
draft_result_paragraph <- function(
    x,
    type = c("enrichment", "ppi", "table"),
    language = c("zh", "en"),
    prompt = NULL,
    model = NULL) {
  language <- match.arg(language)
  model    <- .resolve_model(model)

  if (inherits(x, "tcm_ai_analysis")) {
    context     <- x$context
    type        <- x$input$type
    input_class <- x$input$object_class
  } else {
    type <- match.arg(type)
    compress_fn <- switch(type,
      enrichment = .compress_enrichment,
      ppi        = .compress_ppi,
      table      = .compress_table
    )
    top_n <- if (type == "table") 20 else 10
    context <- compress_fn(x, top_n)
    input_class <- paste(class(x), collapse = "/")
  }

  lang_inst <- if (language == "zh") {
    "Write entirely in Chinese (Simplified)."
  } else {
    "Write entirely in English."
  }

  sys <- paste(
    "You are a scientific writing assistant.",
    lang_inst,
    "Write a concise result paragraph (3-6 sentences).",
    "Cite specific numbers from the data.",
    "Do not invent statistics.",
    "Respond with a single valid JSON object only.",
    sep = "\n"
  )

  prompt <- paste(
    prompt %||% sprintf("Draft a result paragraph from this %s data:", type),
    context, sep = "\n\n"
  )

  result <- .call_generate_object(model, prompt, .draft_schema(), sys,
                                   temperature = 0.4)

  extracted <- .extract_object(result, "draft")
  .new_tcm_ai_draft(
    input    = list(type = type,
                    object_class = input_class),
    context  = context,
    draft    = extracted$output,
    metadata = .build_metadata(model, language,
                               audience = "paper",
                               input_class,
                               output_mode = extracted$output_mode)
  )
}
