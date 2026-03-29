# ai_core.R
# Internal helpers: aisdk dependency check, model resolution, result wrappers.
# aisdk is an optional (Suggests) dependency — missing aisdk never breaks
# non-AI TCMDATA functions.

#' Check whether aisdk is available
#' @return TRUE invisibly if aisdk is installed; otherwise stops with a
#'   user-friendly message.
#' @keywords internal
#' @noRd
.check_aisdk <- function() {
  if (!requireNamespace("aisdk", quietly = TRUE)) {
    stop(
      "The 'aisdk' package is required for AI interpretation functions.\n",
      "Install it with: devtools::install_github('YuLab-SMU/aisdk')",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Resolve a model argument for AI functions
#'
#' If \code{model} is NULL, delegates to \code{aisdk::get_model()} which
#' reads the package-wide default (set via \code{aisdk::set_model()}).
#'
#' @param model NULL, a string ID like \code{"openai:gpt-4o"}, or a
#'   LanguageModelV1 object.
#' @return A resolved model suitable for \code{aisdk::generate_object()}.
#' @keywords internal
#' @noRd
.resolve_model <- function(model = NULL) {
  .check_aisdk()
  if (is.null(model)) {
    model <- aisdk::get_model()
  }
  if (is.null(model)) {
    stop("No model configured. Run tcm_setup() first.", call. = FALSE)
  }
  return(model)
}

#' Build a tcm_ai_analysis S3 object
#' @keywords internal
#' @noRd
.new_tcm_ai_analysis <- function(input, context, output, metadata) {
  structure(
    list(
      input    = input,
      context  = context,
      output   = output,
      metadata = metadata
    ),
    class = c("tcm_ai_analysis", "tcm_ai_result")
  )
}

#' Build a tcm_ai_draft S3 object
#' @keywords internal
#' @noRd
.new_tcm_ai_draft <- function(input, context, draft, metadata) {
  structure(
    list(
      input    = input,
      context  = context,
      draft    = draft,
      metadata = metadata
    ),
    class = c("tcm_ai_draft", "tcm_ai_result")
  )
}

#' Build a tcm_ai_custom S3 object
#' @keywords internal
#' @noRd
.new_tcm_ai_custom <- function(input, context, output, metadata, schema) {
  structure(
    list(
      input    = input,
      context  = context,
      output   = output,
      metadata = metadata,
      schema   = schema
    ),
    class = c("tcm_ai_custom", "tcm_ai_result")
  )
}

#' Extract the object from a generate_object() result for a custom schema,
#' with fallback parsing. Returns list(output = ..., output_mode = ...) where
#' output_mode is "structured" or "fallback_text".
#' @keywords internal
#' @noRd
.extract_custom_object <- function(result) {
  if (!is.null(result$object)) {
    return(list(output = result$object, output_mode = "structured"))
  }

  raw <- trimws(result$raw_text %||% "")
  if (!nzchar(raw)) {
    return(list(output = list(raw_text = ""), output_mode = "fallback_text"))
  }

  # Strip markdown ```json ... ``` fences if present
  clean <- gsub("^```(?:json)?[[:space:]]*|[[:space:]]*```$", "",
                raw, perl = TRUE)

  parsed <- tryCatch(
    jsonlite::fromJSON(clean, simplifyVector = FALSE),
    error = function(e) NULL
  )

  if (!is.null(parsed)) {
    return(list(output = parsed, output_mode = "structured"))
  }
  list(output = list(raw_text = raw), output_mode = "fallback_text")
}

#' Build standard metadata list
#' @keywords internal
#' @noRd
.build_metadata <- function(model, language, audience, input_class,
                            prompt_version = "1.0",
                            output_mode    = "structured") {
  model_id <- if (is.character(model)) {
    model
  } else if (inherits(model, "LanguageModelV1")) {
    prov <- if (is.null(model$provider)) "unknown" else model$provider
    mid  <- if (is.null(model$model_id)) "unknown" else model$model_id
    paste0(prov, ":", mid)
  } else {
    "unknown"
  }

  list(
    model          = model_id,
    language       = language,
    audience       = audience,
    input_class    = input_class,
    generated_at   = Sys.time(),
    prompt_version = prompt_version,
    output_mode    = output_mode
  )
}

# Hardcoded fallback list — used when aisdk is not yet installed
# (e.g. during tcm_config() which does not require aisdk).
# "custom" is intentionally excluded: aisdk provides create_custom_provider(),
# not create_custom(), so advertising it causes misleading failures in tcm_setup().
.tcm_providers_fallback <- c(
  "openai", "anthropic", "gemini", "deepseek", "deepseek_anthropic",
  "volcengine", "stepfun", "openrouter", "xai", "nvidia", "bailian",
  "aihubmix", "aihubmix_anthropic", "aihubmix_gemini"
)

# Return the validated provider whitelist.
# Dynamic scanning of aisdk create_*() functions is intentionally avoided:
# it would also match create_agent(), create_skill_registry(), etc.
.available_providers <- function() {
  .tcm_providers_fallback
}

# Providers whose create_* function does not accept base_url.
.tcm_no_base_url <- c(
  "deepseek_anthropic", "aihubmix_anthropic", "aihubmix_gemini"
)

#' Write AI provider credentials to .env
#'
#' Saves \code{TCM_PROVIDER}, \code{TCM_API_KEY}, \code{TCM_MODEL}, and
#' optionally \code{TCM_BASE_URL} to a \code{.env} file. Existing values for
#' these four keys are overwritten; all other lines are preserved.
#'
#' Supported \code{provider} values: \code{"openai"}, \code{"anthropic"},
#' \code{"gemini"}, \code{"deepseek"}, \code{"deepseek_anthropic"},
#' \code{"volcengine"}, \code{"stepfun"}, \code{"openrouter"}, \code{"xai"},
#' \code{"nvidia"}, \code{"bailian"}, \code{"aihubmix"},
#' \code{"aihubmix_anthropic"}, \code{"aihubmix_gemini"}.
#'
#' @param provider Character. Provider name (see Details).
#' @param api_key Character. Your API key.
#' @param model Character. Model name, e.g. \code{"gpt-4o-mini"},
#'   \code{"claude-3-5-haiku-20241022"}, \code{"gemini-2.0-flash"}.
#' @param base_url Character or NULL. Override the default API endpoint.
#'   Required for proxies or self-hosted endpoints.
#' @param path Character. Path to the \code{.env} file. Default \code{".env"}.
#'
#' @return The path to the \code{.env} file, invisibly.
#' @examples
#' \dontrun{
#'   tcm_config("openai",    "sk-xxx",     "gpt-4o-mini")
#'   tcm_config("anthropic", "sk-ant-xxx", "claude-3-5-haiku-20241022")
#'   tcm_config("gemini",    "AIza-xxx",   "gemini-2.0-flash")
#'   tcm_config("deepseek",  "sk-xxx",     "deepseek-chat")
#'   tcm_config("openai",    "sk-xxx",     "gpt-5-minimal",
#'              base_url = "https://www.packyapi.com/v1")
#' }
#' @export
tcm_config <- function(provider,
                       api_key,
                       model,
                       base_url = NULL,
                       path = ".env") {
  if (!provider %in% .available_providers()) {
    stop(sprintf("Unknown provider '%s'. Supported: %s",
                 provider, paste(.available_providers(), collapse = ", ")),
         call. = FALSE)
  }

  env_path <- if (dirname(path) == ".") file.path(getwd(), path) else path
  existing <- if (file.exists(env_path)) readLines(env_path) else character(0)
  managed  <- c("TCM_PROVIDER", "TCM_API_KEY", "TCM_MODEL", "TCM_BASE_URL")
  kept <- existing[
    !grepl(paste0("^(", paste(managed, collapse = "|"), ")="), existing)
  ]

  new_lines <- c(
    kept,
    sprintf("TCM_PROVIDER=%s", provider),
    sprintf("TCM_API_KEY=%s",  api_key),
    sprintf("TCM_MODEL=%s",    model)
  )
  if (!is.null(base_url) && nzchar(base_url)) {
    new_lines <- c(new_lines, sprintf("TCM_BASE_URL=%s", base_url))
  }

  writeLines(new_lines, env_path)
  message(sprintf("tcm_config: saved [provider=%s, model=%s] -> %s",
                  provider, model, env_path))
  invisible(env_path)
}

#' Initialise the AI model from .env or explicit arguments
#'
#' Loads \code{.env} (if present), resolves \code{TCM_*} variables, calls
#' the matching \code{aisdk::create_*()} function, and registers the model via
#' \code{aisdk::set_model()}. All subsequent AI functions then work without
#' further setup.
#'
#' @param provider Character or NULL. Overrides \code{TCM_PROVIDER}.
#' @param api_key Character or NULL. Overrides \code{TCM_API_KEY}.
#' @param model Character or NULL. Overrides \code{TCM_MODEL}.
#' @param base_url Character or NULL. Overrides \code{TCM_BASE_URL}.
#' @param .env Logical. Load \code{.env} before reading env vars (default TRUE).
#' @param save Logical. If TRUE, also calls \code{\link{tcm_config}()} to
#'   persist the resolved credentials to \code{.env}. Useful for first-time
#'   setup when you want a single call to both initialise and save.
#'   Default FALSE.
#' @param test Logical. If TRUE, sends a minimal test request after setup to
#'   verify the API key and endpoint are reachable. Warnings (not errors) are
#'   issued on failure so the model is still registered. Default FALSE.
#'
#' @return The model object, invisibly.
#' @examples
#' \dontrun{
#'   # Standard two-step workflow
#'   tcm_config("openai", "sk-xxx", "gpt-4o-mini")
#'   tcm_setup()
#'
#'   # One-step: configure + initialise in a single call
#'   tcm_setup("openai", "sk-xxx", "gpt-4o-mini", save = TRUE)
#'
#'   # Verify connectivity after setup
#'   tcm_setup(test = TRUE)
#'
#'   # Override at runtime without touching .env
#'   tcm_setup("deepseek", api_key = "sk-xxx", model = "deepseek-chat")
#' }
#' @export
tcm_setup <- function(provider = NULL,
                      api_key  = NULL,
                      model    = NULL,
                      base_url = NULL,
                      .env     = TRUE,
                      save     = FALSE,
                      test     = FALSE) {
  .check_aisdk()

  if (.env && requireNamespace("dotenv", quietly = TRUE)) {
    env_file <- file.path(getwd(), ".env")
    if (file.exists(env_file)) dotenv::load_dot_env(env_file)
  }

  provider <- .env_or(provider, Sys.getenv("TCM_PROVIDER", "openai"))
  api_key  <- .env_or(api_key,  Sys.getenv("TCM_API_KEY",  ""))
  model    <- .env_or(model,    Sys.getenv("TCM_MODEL",    "gpt-4o-mini"))
  base_url <- .env_or(base_url, Sys.getenv("TCM_BASE_URL", ""))

  if (!nzchar(api_key)) {
    env_path <- file.path(getwd(), ".env")
    stop(
      "No API key found.\n",
      "  Option A: pass it directly —\n",
      "    tcm_setup(provider=\"openai\", api_key=\"sk-...\",",
      " model=\"gpt-4o-mini\")\n",
      "  Option B: save to .env first —\n",
      "    tcm_config(\"openai\", \"sk-...\", \"gpt-4o-mini\")",
      " then tcm_setup()\n",
      "  .env searched at: ", env_path, "\n",
      "  (tcm_config() writes to .env in the current working directory)",
      call. = FALSE
    )
  }

  if (save) {
    tcm_config(
      provider = provider,
      api_key  = api_key,
      model    = model,
      base_url = if (nzchar(base_url)) base_url else NULL
    )
  }

  create_fn_name <- paste0("create_", provider)
  if (!exists(create_fn_name, envir = getNamespace("aisdk"),
              inherits = FALSE)) {
    stop(
      sprintf("No aisdk function '%s'. Supported providers: %s",
              create_fn_name,
              paste(.available_providers(), collapse = ", ")),
      call. = FALSE
    )
  }
  create_fn <- getFromNamespace(create_fn_name, "aisdk")

  args <- list(api_key = api_key)
  if (nzchar(base_url) && !provider %in% .tcm_no_base_url) {
    args$base_url <- base_url
  }

  model_obj <- do.call(create_fn, args)$language_model(model)
  aisdk::set_model(model_obj)

  if (test) {
    tryCatch(
      aisdk::generate_text(model = model_obj, prompt = "ping"),
      error = function(e) warning(
        "Connection test failed: ", conditionMessage(e), call. = FALSE
      )
    )
  }

  message(sprintf("tcm_setup: %s / %s ready.", provider, model))
  invisible(model_obj)
}

# null-coalescing operator (NULL only — consistent with rlang::`%||%`)
`%||%` <- function(x, y) if (is.null(x)) y else x

# env-var helper: NULL *or* empty string falls back to default.
# Used exclusively inside tcm_setup() for Sys.getenv() resolution.
.env_or <- function(x, default) if (is.null(x) || !nzchar(x)) default else x

#' Extract structured output from a generate_object result
#'
#' Tries three strategies in order:
#' 1. Native structured object from aisdk (\code{result$object}).
#' 2. JSON parsed from \code{result$raw_text} (models that return JSON text
#'    but do not support the structured-output API).
#' 3. Plain-text fallback: wraps \code{raw_text} into \code{$summary} so the
#'    result is never NULL, even for models that ignore schema instructions.
#'
#' @param result A \code{GenerateObjectResult} from
#'   \code{aisdk::generate_object()}.
#' @param type One of \code{"analysis"} or \code{"draft"}, controls the shape
#'   of the plain-text fallback object.
#' @return A list with two elements: \code{$output} (a named list matching the
#'   schema, never NULL) and \code{$output_mode} ("structured" when the model
#'   returned valid JSON, "fallback_text" otherwise).
#' @keywords internal
#' @noRd
.extract_object <- function(result, type = "analysis") {
  # Case 1: model returned a proper structured object
  if (!is.null(result$object)) {
    return(list(output = result$object, output_mode = "structured"))
  }

  raw <- trimws(result$raw_text %||% "")

  # Case 2: empty response
  if (!nzchar(raw)) {
    return(list(output = .fallback_object("", type),
                output_mode = "fallback_text"))
  }

  # Case 3: try stripping markdown code fences then parsing JSON
  clean <- gsub("^```(?:json)?[[:space:]]*|[[:space:]]*```$", "",
                raw, perl = TRUE)
  parsed <- tryCatch(
    jsonlite::fromJSON(clean, simplifyVector = FALSE),
    error = function(e) NULL
  )
  if (!is.null(parsed)) {
    return(list(output = parsed, output_mode = "structured"))
  }

  # Case 4: plain-text fallback — model ignored the schema
  list(output = .fallback_object(raw, type), output_mode = "fallback_text")
}

#' Build a minimal fallback object from plain text
#' @keywords internal
#' @noRd
.fallback_object <- function(text, type) {
  if (type == "draft") {
    list(paragraph = text, figure_legend_hint = "")
  } else {
    list(
      summary                   = text,
      key_findings              = list(),
      biological_interpretation = "",
      tcm_relevance             = "",
      caveats                   = list()
    )
  }
}
