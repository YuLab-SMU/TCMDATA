# Fixed z_object() schemas used by generate_object() calls.
# All interpretation functions share the same output contract.
#' Schema for structured analysis interpretation output
#' @keywords internal
#' @noRd
.analysis_schema <- function() {
  aisdk::z_object(
    summary = aisdk::z_string(
      "A concise 2-4 sentence overview of the analysis results"
    ),
    key_findings = aisdk::z_array(
      aisdk::z_string(),
      "List of key findings directly supported by the data"
    ),
    biological_interpretation = aisdk::z_string(
      "Biological interpretation of the results in context"
    ),
    tcm_relevance = aisdk::z_string(
      paste(
        "Relevance to TCM herb-target pharmacology or Traditional Chinese",
        "Medicine theory (e.g. which herbs or classic formulas target the",
        "identified genes/pathways). Leave empty string if not applicable."
      )
    ),
    caveats = aisdk::z_array(
      aisdk::z_string(),
      "Methodological caveats, limitations, or potential confounders"
    )
  )
}

#' Schema for result paragraph drafting
#' @keywords internal
#' @noRd
.draft_schema <- function() {
  aisdk::z_object(
    paragraph = aisdk::z_string(
      "A publication-ready result paragraph (3-6 sentences)"
    ),
    figure_legend_hint = aisdk::z_string(
      "A brief suggestion for a matching figure legend (1-2 sentences)"
    )
  )
}

# ── Public schema builder API ─────────────────────────────────────────────────
# Thin wrappers around aisdk::z_*() so callers never need to load aisdk.

#' Build a custom output schema for \code{tcm_interpret_schema()}
#'
#' A family of thin wrappers around \code{aisdk::z_object()} and the
#' associated field constructors. These functions let you define structured
#' output contracts without writing \code{aisdk::} anywhere in your code.
#'
#' \describe{
#'   \item{\code{tcm_schema(...)}}{Container: wraps \code{aisdk::z_object()}.
#'     Pass named \code{tcm_field_*()} calls as arguments.}
#'   \item{\code{tcm_field_string(description)}}{A single string field.}
#'   \item{\code{tcm_field_number(description)}}{A single numeric field.}
#'   \item{\code{tcm_field_boolean(description)}}{A single boolean field.}
#'   \item{\code{tcm_field_array(description)}}{An array of strings.}
#'   \item{\code{tcm_field_enum(choices, description)}}{One of a fixed set of
#'     string values; \code{choices} is a character vector.}
#' }
#'
#' @param ... Named field definitions created with \code{tcm_field_*()}.
#' @param description Character. Guidance for the model about this field.
#'   Default \code{""} (no description).
#' @param choices Character vector of allowed values (for \code{tcm_field_enum}
#'   only).
#'
#' @return A schema object (class \code{z_schema} from aisdk) suitable for
#'   \code{\link{tcm_interpret_schema}}.
#'
#' @examples
#' \dontrun{
#'   my_schema <- tcm_schema(
#'     summary     = tcm_field_string("2-3 sentence overview"),
#'     mechanism   = tcm_field_string("Key molecular mechanism"),
#'     key_targets = tcm_field_array("Top gene or protein targets"),
#'     score       = tcm_field_number("Confidence score 0-1"),
#'     is_tcm      = tcm_field_boolean("Whether TCM herbs are involved"),
#'     confidence  = tcm_field_enum(c("high", "medium", "low"))
#'   )
#'
#'   res <- tcm_interpret_schema(enrich_res, schema = my_schema,
#'                               language = "zh")
#'   res$output$summary
#'   res$output$key_targets   # character vector
#'   res$output$confidence    # one of "high" / "medium" / "low"
#' }
#' @name tcm_schema
#' @export
tcm_schema <- function(...) {
  .check_aisdk()
  aisdk::z_object(...)
}

#' @rdname tcm_schema
#' @export
tcm_field_string <- function(description = "") {
  .check_aisdk()
  aisdk::z_string(description)
}

#' @rdname tcm_schema
#' @export
tcm_field_number <- function(description = "") {
  .check_aisdk()
  aisdk::z_number(description)
}

#' @rdname tcm_schema
#' @export
tcm_field_boolean <- function(description = "") {
  .check_aisdk()
  aisdk::z_boolean(description)
}

#' @rdname tcm_schema
#' @export
tcm_field_array <- function(description = "") {
  .check_aisdk()
  aisdk::z_array(aisdk::z_string(), description)
}

#' @rdname tcm_schema
#' @export
tcm_field_enum <- function(choices, description = "") {
  .check_aisdk()
  aisdk::z_enum(choices)
}
