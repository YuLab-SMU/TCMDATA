#' Create a bioinformatics interpreter agent in TCMDATA
#'
#' A thin wrapper around \code{aisdk::create_agent()}. Pass the
#' result to \code{\link{run_tcm_agent}} to call it like a regular
#' function.
#'
#' @param name Character. A short identifier for the agent.
#' @param description Character. One-line description of the agent's role.
#' @param system_prompt Character. The system prompt defining agent behaviour.
#'
#' @return An aisdk agent object.
#' @examples
#' \dontrun{
#'   bio_agent <- create_tcm_agent(
#'     name = "BioInterpreter",
#'     description = "Bioinformatics result interpreter",
#'     system_prompt = "You are a bioinformatics expert. Explain results
#'                      concisely in 3 sentences for a research report."
#'   )
#' }
#' @export
create_tcm_agent <- function(name, description, system_prompt) {
  .check_aisdk()
  aisdk::create_agent(
    name = name,
    description = description,
    system_prompt = system_prompt
  )
}

#' Create a reusable AI function from an instruction
#'
#' A function factory: supply an instruction (system prompt) once, get back a
#' plain R function. The returned function accepts a query string and behaves
#' like any ordinary R function — it works directly in \code{sapply()},
#' \code{purrr::map()}, \code{mutate()}, and other vectorised workflows.
#'
#' Compared to calling \code{\link{tcm_interpret}} with a custom \code{system}
#' parameter each time, \code{make_tcm_function} fixes the role once and
#' returns a named, reusable callable — no need to re-pass parameters on every
#' call. Compared to \code{\link{create_tcm_agent}} +
#' \code{\link{run_tcm_agent}}, the result is a plain function rather than an
#' agent object, so it composes directly with vectorised R idioms.
#'
#' @param instruction Character. The system prompt that defines the function's
#'   behaviour (its "role" and task description).
#' @param name Character. A short identifier used for the underlying agent.
#'   Useful for debugging. Default \code{"CustomFn"}.
#'
#' @return A function with signature
#'   \code{function(input, model = NULL, verbose = TRUE)} that calls the agent
#'   with \code{input} as the task. When \code{verbose = TRUE} (default) the
#'   response is printed via \code{cat()} and returned invisibly; when
#'   \code{FALSE} it is returned visibly for assignment or pipeline use.
#'
#' @examples
#' \dontrun{
#'   # One-liner function creation
#'   explain_pathway <- make_tcm_function(
#'     "你是MSigDB/KEGG通路专家，用2句中文解释通路的生物学功能，适合科研报告。"
#'   )
#'
#'   # Call like a normal R function
#'   explain_pathway("HALLMARK_INFLAMMATORY_RESPONSE")
#'
#'   # Direct drop-in for sapply() — no wrapper needed
#'   pathways <- c("HALLMARK_E2F_TARGETS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
#'   explanations <- sapply(pathways, explain_pathway)
#'   data.frame(pathway = pathways, explanation = explanations)
#'
#'   # Create multiple domain-specific functions from the same pattern
#'   draft_legend <- make_tcm_function(
#'     "Write a concise figure legend in English for a journal submission.",
#'     name = "LegendWriter"
#'   )
#'   review_method <- make_tcm_function(
#'     "作为严格的同行评审人，列出该分析段落中潜在的方法学问题。",
#'     name = "PeerReviewer"
#'   )
#' }
#' @export
make_tcm_function <- function(instruction, name = "CustomFn") {
  .check_aisdk()
  agent <- create_tcm_agent(
    name          = name,
    description   = name,
    system_prompt = instruction
  )
  function(input, model = NULL, verbose = TRUE) {
    run_tcm_agent(agent, input, model = model, verbose = verbose)
  }
}

#' Run an agent as a plain function
#'
#' Wraps \code{agent$run()} so the agent behaves like an ordinary R function.
#' Prints the response via \code{cat()} and returns it invisibly.
#'
#' @param agent An agent object from \code{\link{create_tcm_agent}}.
#' @param query Character. The query or task to send to the agent.
#' @param model A model identifier, LanguageModelV1 object, or NULL (uses the
#'   package-wide default set via \code{aisdk::set_model()}).
#'
#' @param verbose Logical. If \code{TRUE} (default), print the response to the
#'   console. Set \code{FALSE} for programmatic / batch use; in that case the
#'   response is returned visibly so it can be assigned or passed downstream.
#'
#' @return The response text (character). Returned invisibly when
#'   \code{verbose = TRUE}, visibly when \code{verbose = FALSE}.
#' @examples
#' \dontrun{
#'   bio_agent <- create_tcm_agent(
#'     name = "BioInterpreter",
#'     description = "Bioinformatics result interpreter",
#'     system_prompt = "You are a bioinformatics expert. Explain results
#'                      concisely in 3 sentences."
#'   )
#'
#'   # Interactive use — prints to console
#'   run_tcm_agent(bio_agent,
#'     "What does HALLMARK_INFLAMMATORY_RESPONSE mean in TME?")
#'
#'   # Batch with sapply() — suppress printing, collect results
#'   pathways <- c("HALLMARK_E2F_TARGETS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
#'   res <- sapply(pathways,
#'     function(p) run_tcm_agent(bio_agent, p, verbose = FALSE))
#' }
#' @export
run_tcm_agent <- function(agent, query, model = NULL, verbose = TRUE) {
  .check_aisdk()
  model <- .resolve_model(model)
  result <- agent$run(task = query, model = model)
  if (verbose) {
    cat(result$text, "\n")
    invisible(result$text)
  } else {
    result$text
  }
}
