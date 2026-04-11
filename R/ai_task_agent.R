# ai_task_agent.R
# Main entry points for TCM agent layer.
# Provides tcm_agent() and related functions for natural language task execution.

#' Create a TCM task agent
#'
#' Creates an aisdk agent configured for TCM network pharmacology analysis.
#' The agent can call TCMDATA tools natively.
#'
#' @param tools List of Tool objects. If NULL, uses default tools from
#'   \code{\link{create_tcm_tools}}.
#' @param system_prompt Character. Custom system prompt. If NULL, uses default.
#' @param skills Character vector of skill directories. If \code{NULL}, uses
#'   the active TCMDATA skill directory plus aisdk's
#'   \code{skill-creator} skill (if available). Use \code{character(0)} to
#'   disable skills entirely.
#'
#' @return An aisdk Agent object.
#'
#' @examples
#' \dontrun{
#'   agent <- create_tcm_task_agent()
#'   result <- run_tcm_task(agent, "Search targets of Astragalus")
#' }
#' @export
create_tcm_task_agent <- function(tools = NULL, system_prompt = NULL,
                                  skills = NULL) {
  .check_aisdk()

  if (is.null(tools)) {
    tools <- create_tcm_tools()
  }

  if (is.null(system_prompt)) {
    system_prompt <- .default_tcm_system_prompt()
  }

  # Resolve skills path
  skill_paths <- NULL
  if (!is.null(skills)) {
    skill_paths <- skills
  } else {
    creator_skill <- tryCatch(tcm_aisdk_skill(), error = function(e) NULL)
    skill_paths <- unique(c(tcm_skill_dir(), creator_skill))
  }

  aisdk::create_agent(
    name = "TCMTaskAgent",
    description = "Traditional Chinese Medicine network pharmacology analysis agent",
    system_prompt = system_prompt,
    tools = tools,
    skills = skill_paths
  )
}

#' Run a task with a TCM agent
#'
#' Executes a natural language task using the provided agent.
#'
#' @param agent An Agent object from \code{\link{create_tcm_task_agent}}.
#' @param task Character. The task description in natural language.
#' @param model Model object or NULL (uses package default).
#' @param verbose Logical. Print progress and results (default TRUE).
#'
#' @return A list with:
#'   \item{text}{The agent's final response text.}
#'   \item{artifacts}{List of artifact handles created during execution.}
#'   \item{tool_calls}{Log of tool calls made.}
#'
#' @examples
#' \dontrun{
#'   agent <- create_tcm_task_agent()
#'   result <- run_tcm_task(agent, "Run herb enrichment on Astragalus targets")
#' }
#' @export
run_tcm_task <- function(agent, task, model = NULL, verbose = TRUE) {
  .check_aisdk()
  model <- .resolve_model(model)

  # Snapshot artifact IDs before the run
  arts_before <- list_tcm_artifacts()$artifact_id

  # Resolve artifact references
  resolved <- resolve_artifact_references(task)
  if (length(resolved$artifact_ids) > 0 && verbose) {
    message("Resolved artifacts: ", paste(resolved$artifact_ids, collapse = ", "))
  }

  # Run the agent
  result <- agent$run(task = resolved$resolved_task, model = model)

  if (verbose) {
    cat("\n", result$text, "\n")
  }

  # Collect artifacts created during this run
  artifacts <- list_tcm_artifacts()

  # Auto-export only NEW artifacts to global env for RStudio visibility
  new_ids <- setdiff(artifacts$artifact_id, arts_before)
  for (aid in new_ids) {
    tryCatch(
      assign(aid, load_tcm_artifact(aid), envir = globalenv()),
      error = function(e) NULL
    )
  }

  invisible(list(
    text       = result$text,
    artifacts  = artifacts,
    tool_calls = result$tool_calls %||% list()
  ))
}

#' One-step TCM task execution
#'
#' The main entry point for natural language TCM analysis. Automatically
#' routes the task, creates appropriate agent, and executes the task.
#'
#' @param task Character. Natural language task description.
#' @param model Model object or NULL (uses package default).
#' @param verbose Logical. Print progress and results (default TRUE).
#' @param use_router Logical. Use task router to select tools (default TRUE).
#'
#' @return A list with text response, artifacts, and tool call log.
#'
#' @examples
#' \dontrun{
#'   # Basic usage
#'   tcm_agent("Search targets of Astragalus")
#'
#'   # Herb enrichment workflow
#'   tcm_agent("Run herb enrichment on Astragalus targets and find top hub genes")
#'
#'   # PPI analysis
#'   tcm_agent("Analyse PPI network for diabetes-related genes and find hub genes")
#'
#'   # Interpretation
#'   tcm_agent("Interpret the previous enrichment result")
#' }
#' @export
tcm_agent <- function(task,
                      model = NULL,
                      verbose = TRUE,
                      use_router = TRUE) {
  .check_aisdk()
  model <- .resolve_model(model)

  # Route task to get appropriate tools
  if (use_router) {
    routing <- route_tcm_task(task, use_llm = FALSE)
    if (verbose) {
      message(sprintf("Task type: %s (confidence: %s)",
                      routing$task_type, routing$confidence))
    }
    tools <- create_tcm_tools(tool_names = routing$tools)
  } else {
    tools <- create_tcm_tools()
  }

  # Only load skills for ambiguous / workflow-class tasks;
  # simple routed tasks (herb_lookup, enrichment, ...) don't need the
  # full network-pharmacology skill which can mislead the model.
  use_skills <- if (use_router) {
    routing$confidence == "low" || routing$task_type == "general"
  } else {
    TRUE
  }
  agent <- create_tcm_task_agent(
    tools  = tools,
    skills = if (use_skills) NULL else character(0)
  )

  # Execute task
  run_tcm_task(agent, task, model = model, verbose = verbose)
}

#' Interactive TCM analysis chat
#'
#' Starts a rich interactive chat session for multi-turn TCM analysis
#' conversations in the terminal. Each exchange is formatted with clear
#' visual structure including routing information, tool call logs, artifact
#' updates, and the agent response.
#'
#' The session automatically routes each task to the appropriate tool set,
#' so users do not need to call individual functions.
#'
#' Built-in commands (type directly at the prompt):
#' \itemize{
#'   \item \code{/help} -- show available commands
#'   \item \code{/artifacts} -- list stored analysis artifacts
#'   \item \code{/last} -- reprint the last agent response
#'   \item \code{/history} -- show conversation history
#'   \item \code{/model} -- show current model info
#'   \item \code{/stats} -- show session statistics
#'   \item \code{/stream} -- toggle streaming mode on/off
#'   \item \code{/clear} -- remove all stored artifacts
#'   \item \code{/quit} / \code{/exit} / \code{/q} -- end the session
#' }
#'
#' @param model Model object or NULL.
#' @param verbose Logical. Print agent responses (default TRUE).
#' @param stream Logical. Enable streaming output (default TRUE). Toggle
#'   at runtime with the \code{/stream} command.
#' @param skills Character vector of skill directories to load for the chat
#'   session. Default \code{NULL} uses the active TCMDATA skill directory
#'   (scanning all skills under it) plus aisdk's \code{skill-creator} skill
#'   if available. Use \code{character(0)} to disable skills, or pass an
#'   explicit vector to override the default.
#'
#' @return Invisibly returns a list with \code{history} (each turn's task and
#'   reply) and \code{artifacts} (a data.frame snapshot from
#'   \code{\link{list_tcm_artifacts}}). Assign the result to keep the full
#'   session record: \code{res <- tcm_chat()}.
#'
#' @examples
#' \dontrun{
#'   tcm_chat()
#'   # With streaming disabled
#'   tcm_chat(stream = FALSE)
#'   # Add aisdk's skill-creator alongside TCMDATA skills
#'   tcm_chat(skills = c(tcm_skill_dir(), tcm_aisdk_skill()))
#' }
#' @export
tcm_chat <- function(model = NULL, verbose = TRUE, stream = TRUE,
                     skills = NULL) {
  .check_aisdk()
  model <- .resolve_model(model)

  # -- State --
  turn <- 0L
  history <- list()
  last_result <- NULL
  total_tools_used <- 0L
  session_start <- Sys.time()
  workflow_log <- character(0)

  # -- Formatting helpers --
  ruler   <- function(ch = "\u2500", n = 62) strrep(ch, n)
  header  <- function(text) {
    cat(ruler(), "\n")
    cat(sprintf("  %s\n", text))
    cat(ruler(), "\n")
  }
  dim_text   <- function(text) paste0("\033[2m", text, "\033[0m")
  bold_text  <- function(text) paste0("\033[1m", text, "\033[0m")
  cyan_text  <- function(text) paste0("\033[36m", text, "\033[0m")
  green_text <- function(text) paste0("\033[32m", text, "\033[0m")
  yellow_text <- function(text) paste0("\033[33m", text, "\033[0m")
  red_text   <- function(text) paste0("\033[31m", text, "\033[0m")
  agent_name <- "TCM-Pharmacist"

  # -- Model info --
  model_id <- tryCatch({
    if (is.character(model)) model
    else if (!is.null(model) && !is.null(model$model_id)) model$model_id
    else if (!is.null(model) && !is.null(model$modelId)) model$modelId
    else Sys.getenv("AISDK_MODEL", unset = "(default)")
  }, error = function(e) "(unknown)")

  # -- Welcome banner --
  all_tools <- create_tcm_tools()
  n_all_tools <- length(all_tools)

  # -- Create persistent chat session --
  session_agent <- create_tcm_task_agent(tools = all_tools, skills = skills)
  chat_session <- aisdk::create_chat_session(agent = session_agent, model = model)

  cat("\n")
  cat(cyan_text(ruler("\u2550")), "\n")
  cat(cyan_text("\u2551"), bold_text(paste0("  ", agent_name, "  \u2014  TCMDATA Analysis Chat")), "\n")
  cat(cyan_text(ruler("\u2550")), "\n")
  cat(dim_text(sprintf("  Model: %s  |  Tools: %d  |  Stream: %s",
                       model_id, n_all_tools, if (stream) "on" else "off")), "\n")
  cat("\n")
  cat("  Describe your analysis task in natural language.\n")
  cat("  The agent will route, execute tools, and return results.\n")
  cat("\n")
  cat(dim_text("  Commands:  /help  /artifacts  /history  /last  /save  /model  /stats  /clear  /quit"), "\n")
  cat(cyan_text(ruler("\u2550")), "\n\n")

  repeat {
    task <- readline(prompt = green_text("> "))
    task <- trimws(task)

    if (nchar(task) == 0) next

    # -- Built-in commands (slash style) --
    cmd <- tolower(task)

    if (cmd %in% c("/quit", "/exit", "/q", "quit", "exit", "q")) {
      cat("\n")
      header("Session ended. Goodbye!")
      break
    }

    if (cmd %in% c("/help", "help")) {
      cat("\n")
      header("Available commands")
      cat("  /help        Show this help message\n")
      cat("  /artifacts   List all stored analysis artifacts\n")
      cat("  /history     Show conversation history\n")
      cat("  /last        Reprint the last agent response\n")
      cat("  /save [WxH]  Export all artifacts (e.g. /save 10x8)\n")
      cat("  /model       Show current model info\n")
      cat("  /stats       Show session statistics\n")
      cat("  /stream      Toggle streaming mode on/off\n")
      cat("  /clear       Remove all stored artifacts\n")
      cat("  /quit        End the session\n")
      cat("\n")
      cat("  Anything else is sent directly to the agent.\n")
      cat(ruler(), "\n\n")
      next
    }

    if (cmd %in% c("/artifacts", "artifacts")) {
      cat("\n")
      header("Stored Artifacts")
      art_df <- list_tcm_artifacts()
      if (nrow(art_df) == 0) {
        cat("  (none)\n")
      } else {
        print(art_df, right = FALSE)
      }
      cat(ruler(), "\n\n")
      next
    }

    if (cmd %in% c("/history", "history")) {
      cat("\n")
      header("Conversation History")
      if (length(history) == 0) {
        cat("  (no turns yet)\n")
      } else {
        for (h in history) {
          cat(sprintf("  [#%d] You:   %s\n", h$turn, h$task))
          reply_preview <- substr(h$reply, 1, 80)
          if (nchar(h$reply) > 80) reply_preview <- paste0(reply_preview, "...")
          cat(sprintf("       Agent: %s\n", reply_preview))
        }
      }
      cat(ruler(), "\n\n")
      next
    }

    if (cmd %in% c("/last", "last")) {
      if (is.null(last_result)) {
        cat(dim_text("  (no previous response)\n\n"))
      } else {
        cat("\n")
        header(sprintf("Last response (Turn #%d)", turn))
        cat(last_result$text, "\n")
        cat(ruler(), "\n\n")
      }
      next
    }

    if (grepl("^(/save|save)", cmd)) {
      # Parse optional size: "/save 10x8" or "save 10 8"
      save_args <- trimws(sub("^/?save", "", task))
      save_w <- 8; save_h <- 6  # defaults
      if (nzchar(save_args)) {
        parts <- as.numeric(strsplit(save_args, "[xX ]+")[[1]])
        if (length(parts) >= 2 && !any(is.na(parts))) {
          save_w <- parts[1]; save_h <- parts[2]
        }
      }

      art_df <- list_tcm_artifacts()
      if (nrow(art_df) == 0) {
        cat(dim_text("  No artifacts to save.\n\n"))
      } else {
        out_dir <- file.path(getwd(), "tcm_output")
        if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

        for (aid in art_df$artifact_id) {
          obj <- load_tcm_artifact(aid)
          assign(aid, obj, envir = globalenv())

          # ggplot -> PNG + PDF
          if (inherits(obj, "ggplot")) {
            png_path <- file.path(out_dir, paste0(aid, ".png"))
            pdf_path <- file.path(out_dir, paste0(aid, ".pdf"))
            tryCatch({
              ggplot2::ggsave(png_path, plot = obj, width = save_w, height = save_h, dpi = 300)
              ggplot2::ggsave(pdf_path, plot = obj, width = save_w, height = save_h)
              cat(dim_text(sprintf("    %s  (ggplot) -> %s, .pdf\n", aid, basename(png_path))))
            }, error = function(e) {
              cat(dim_text(sprintf("    %s  (ggplot) -> export to R env only (save failed: %s)\n",
                                   aid, conditionMessage(e))))
            })

          # data.frame / matrix -> CSV
          } else if (is.data.frame(obj) || is.matrix(obj)) {
            csv_path <- file.path(out_dir, paste0(aid, ".csv"))
            tryCatch({
              utils::write.csv(obj, csv_path, row.names = FALSE)
              cat(dim_text(sprintf("    %s  (%s) -> %s\n", aid, class(obj)[1],
                                   basename(csv_path))))
            }, error = function(e) {
              cat(dim_text(sprintf("    %s  (%s) -> export to R env only\n", aid, class(obj)[1])))
            })

          # igraph -> graphml
          } else if (inherits(obj, "igraph")) {
            gml_path <- file.path(out_dir, paste0(aid, ".graphml"))
            tryCatch({
              igraph::write_graph(obj, gml_path, format = "graphml")
              cat(dim_text(sprintf("    %s  (igraph) -> %s\n", aid, basename(gml_path))))
            }, error = function(e) {
              cat(dim_text(sprintf("    %s  (igraph) -> export to R env only\n", aid)))
            })

          # enrichResult -> CSV of the result table
          } else if (inherits(obj, "enrichResult") || inherits(obj, "gseaResult")) {
            csv_path <- file.path(out_dir, paste0(aid, ".csv"))
            tryCatch({
              utils::write.csv(as.data.frame(obj), csv_path, row.names = FALSE)
              cat(dim_text(sprintf("    %s  (%s) -> %s\n", aid, class(obj)[1],
                                   basename(csv_path))))
            }, error = function(e) {
              cat(dim_text(sprintf("    %s  (%s) -> export to R env only\n", aid, class(obj)[1])))
            })

          # everything else -> R env only
          } else {
            cat(dim_text(sprintf("    %s  (%s) -> export to R env only\n", aid, class(obj)[1])))
          }
        }
        cat(dim_text(sprintf("\n  %d artifact(s) exported to global env.\n", nrow(art_df))))
        cat(dim_text(sprintf("  Files saved to: %s\n\n", out_dir)))
      }
      next
    }

    if (cmd %in% c("/model", "model")) {
      cat("\n")
      header("Model Info")
      cat(sprintf("  Model:    %s\n", model_id))
      cat(sprintf("  Agent:    %s\n", agent_name))
      cat(sprintf("  Tools:    %d available\n", n_all_tools))
      cat(sprintf("  Stream:   %s\n", if (stream) "on" else "off"))
      cat(ruler(), "\n\n")
      next
    }

    if (cmd %in% c("/stats", "stats")) {
      cat("\n")
      header("Session Statistics")
      elapsed <- as.numeric(difftime(Sys.time(), session_start, units = "mins"))
      art_count <- nrow(list_tcm_artifacts())
      cat(sprintf("  Turns:          %d\n", turn))
      cat(sprintf("  Tool calls:     %d\n", total_tools_used))
      cat(sprintf("  Artifacts:      %d\n", art_count))
      cat(sprintf("  Duration:       %.1f min\n", elapsed))
      cat(sprintf("  Model:          %s\n", model_id))
      cat(ruler(), "\n\n")
      next
    }

    if (cmd %in% c("/stream", "stream")) {
      stream <- !stream
      cat(dim_text(sprintf("  Stream: %s\n\n", if (stream) "on" else "off")))
      next
    }

    if (cmd %in% c("/clear", "clear")) {
      clear_tcm_artifacts()
      workflow_log <<- character(0)
      cat(dim_text("  Artifacts and workflow log cleared.\n\n"))
      next
    }

    # -- Agent turn --
    turn <- turn + 1L
    arts_before <- nrow(list_tcm_artifacts())

    cat("\n")
    cat(cyan_text(ruler("\u2504")), "\n")
    cat(dim_text(sprintf("  Turn #%d  |  %s", turn, format(Sys.time(), "%H:%M:%S"))), "\n")

    # Route (informational only -- session has all tools)
    routing <- route_tcm_task(task, use_llm = FALSE)
    cat(dim_text(sprintf("  -> Route: %s  |  confidence: %s  |  source: %s",
                         routing$task_type, routing$confidence, routing$source_hint)), "\n")
    cat(cyan_text(ruler("\u2504")), "\n")

    # Execute via persistent ChatSession
    tryCatch({
      resolved <- resolve_artifact_references(task)
      if (length(resolved$artifact_ids) > 0) {
        cat(dim_text(sprintf("  Linked artifacts: %s",
                             paste(resolved$artifact_ids, collapse = ", "))), "\n")
      }

      prompt <- resolved$resolved_task

      # Inject workflow progress as turn context
      turn_ctx <- NULL
      if (length(workflow_log) > 0) {
        turn_ctx <- paste0(
          "<workflow_progress>\n",
          paste(workflow_log, collapse = "\n"),
          "\n</workflow_progress>"
        )
      }

      # Streaming or batch execution
      if (stream) {
        stream_buf <- character(0)
        cat(bold_text(paste0("[", agent_name, "]")), "\n")
        chat_session$send_stream(
          prompt   = prompt,
          turn_system_prompt = turn_ctx,
          callback = function(chunk, done = FALSE) {
            cat(chunk)
            stream_buf <<- c(stream_buf, chunk)
          }
        )
        cat("\n")
        # Get text from session history
        response_text <- chat_session$get_last_response()
        if (is.null(response_text) || !nzchar(response_text)) {
          response_text <- paste0(stream_buf, collapse = "")
        }
        result <- list(text = response_text, tool_calls = list())
      } else {
        result <- chat_session$send(prompt = prompt, turn_system_prompt = turn_ctx)
      }

      # Tool call log
      tool_calls <- result$tool_calls %||% list()
      total_tools_used <- total_tools_used + length(tool_calls)
      if (length(tool_calls) > 0) {
        cat(dim_text(paste0("  Tools: ",
                             paste(vapply(tool_calls, function(tc) {
                               tc$name %||% tc$tool %||% "?"
                             }, character(1)), collapse = " -> "))), "\n")
      }

      # Artifact diff
      arts_after_df <- list_tcm_artifacts()
      arts_after <- nrow(arts_after_df)
      new_count <- arts_after - arts_before
      if (new_count > 0) {
        new_ids <- utils::tail(arts_after_df$artifact_id, new_count)
        for (nid in new_ids) {
          tryCatch(
            assign(nid, load_tcm_artifact(nid), envir = globalenv()),
            error = function(e) NULL
          )
        }
        cat(yellow_text(paste0("  + ", new_count, " artifact(s)  (total: ",
                             arts_after, ")  [exported: ",
                             paste(new_ids, collapse = ", "), "]")), "\n")
      }

      cat(cyan_text(ruler("\u2504")), "\n")

      # Non-streaming: print response now
      if (!stream) {
        cat("\n")
        cat(bold_text(paste0("[", agent_name, "]")), "\n")
        cat(result$text, "\n")
      }
      cat("\n")

      # Update workflow progress log
      tool_names_used <- vapply(tool_calls, function(tc) {
        tc$name %||% tc$tool %||% "?"
      }, character(1))
      log_parts <- sprintf("[Turn %d] %s", turn, routing$task_type)
      if (length(tool_names_used) > 0) {
        log_parts <- paste0(log_parts, " | tools: ",
                            paste(tool_names_used, collapse = " -> "))
      }
      if (new_count > 0) {
        log_parts <- paste0(log_parts, " | artifacts: ",
                            paste(new_ids, collapse = ", "))
      }
      workflow_log <<- c(workflow_log, log_parts)

      # -- Verification turn for complex workflows --
      # If multiple tools were called and artifacts were created,
      # run a self-review pass to check result quality
      if (length(tool_calls) >= 3 && new_count >= 2) {
        tryCatch({
          # Collect quality warnings from newly created artifacts
          qw <- character(0)
          for (nid in new_ids) {
            art_meta <- tryCatch(
              attr(load_tcm_artifact(nid), "quality_warnings"),
              error = function(e) NULL
            )
            if (!is.null(art_meta)) qw <- c(qw, art_meta)
          }

          verify_prompt <- paste0(
            "<verification_request>\n",
            "You just completed a multi-step analysis (turn #", turn, "). ",
            "Briefly self-review the results:\n",
            "1. Were all requested steps completed?\n",
            "2. Are there any quality warnings to highlight?\n",
            if (length(qw) > 0) paste0("Known warnings: ", paste(qw, collapse = "; "), "\n") else "",
            "3. What are the logical next steps the user might want?\n",
            "Keep the review to 2-3 sentences appended to your response. ",
            "Do NOT repeat the full analysis. Start with '---' on a new line.\n",
            "</verification_request>"
          )

          if (stream) {
            cat(dim_text("  [verification]"), "\n")
            chat_session$send_stream(
              prompt = verify_prompt,
              turn_system_prompt = turn_ctx,
              callback = function(chunk, done = FALSE) {
                cat(chunk)
                stream_buf <<- c(stream_buf, chunk)
              }
            )
            cat("\n")
          } else {
            verify_result <- chat_session$send(
              prompt = verify_prompt,
              turn_system_prompt = turn_ctx
            )
            cat("\n", verify_result$text, "\n")
          }
          workflow_log <<- c(workflow_log,
                             sprintf("[Turn %d] verification pass completed", turn))
        }, error = function(e) {
          # Silent failure for verification -- non-critical
          NULL
        })
      }

      last_result <- result
      history[[length(history) + 1L]] <- list(
        turn  = turn,
        task  = task,
        reply = result$text
      )

    }, error = function(e) {
      cat(cyan_text(ruler("\u2504")), "\n")
      cat(red_text(paste0("  Error: ", conditionMessage(e))), "\n\n")
      last_result <<- NULL
      history[[length(history) + 1L]] <<- list(
        turn  = turn,
        task  = task,
        reply = paste("Error:", conditionMessage(e))
      )
    })
  }

  invisible(list(
    history   = history,
    artifacts = list_tcm_artifacts()
  ))
}

#' Default system prompt for TCM task agent
#' @keywords internal
#' @noRd
.default_tcm_system_prompt <- function() {
  paste(
    "<identity>",
    "You are TCM-Pharmacist, a network pharmacology and bioinformatics analysis assistant",
    "built on the TCMDATA R package. You combine traditional Chinese medicine (TCM)",
    "domain knowledge with modern computational pharmacology methods to help researchers",
    "conduct rigorous, reproducible analyses.",
    "</identity>",
    "",
    "<capabilities>",
    "You have access to the following tool categories. Use them directly to perform analyses:",
    "- Herb & target search: query the TCMDATA herb-compound-target database",
    "- Disease-target search: DisGeNET-based disease-gene associations and reverse lookup",
    "- Target intersection: compute overlapping targets between herb and disease gene sets",
    "- Functional enrichment: GO, KEGG, and herb enrichment via clusterProfiler",
    "- PPI network analysis: STRING-based PPI retrieval, topological metrics, hub gene ranking, and MCODE module detection",
    "- Machine learning screening: multi-algorithm feature selection (LASSO, RF, SVM, XGBoost, etc.) and consensus gene identification",
    "- Literature mining: PubMed evidence retrieval and publication trend analysis",
    "- Compound annotation: PubChem CID resolution, molecular properties, and structural similarity",
    "- Visualization: publication-ready plots (lollipop, heatmap, sankey, ROC, Venn, UpSet, radar, docking heatmap)",
    "- Result interpretation: AI-assisted biological interpretation of stored artifacts",
    "</capabilities>",
    "",
    "<execution_policy>",
    "Execute tasks immediately using tools. Use sensible defaults for unspecified parameters.",
    "Do NOT ask for confirmation -- call the tool(s) and return results right away.",
    "",
    "SCOPE RULE: Complete EVERYTHING the user asked for -- no more, no less.",
    "- Single-step: 'search Astragalus targets' -> ONLY search, do NOT enrich or plot.",
    "- Multi-step: 'find top 50 targets of Astragalus, then do GO enrichment' -> do BOTH in one turn.",
    "- Full pipeline: 'do network pharmacology analysis' -> chain all relevant steps.",
    "Key signals for multi-step: '\u7136\u540e', '\u63a5\u7740', '\u5e76\u4e14', 'and then', 'then do', comma-separated actions.",
    "When in doubt whether user wants step B, do it if it was explicitly mentioned.",
    "",
    "Only ask a question when truly essential information is missing AND has no reasonable default.",
    "Example of NO need to ask: user says 'search targets of huangqi' -- just call the tool.",
    "Example of MUST ask: user says 'run enrichment' but gives no gene list and no previous artifact.",
    "</execution_policy>",
    "",
    "<tool_use_guidelines>",
    "- Call tools directly instead of suggesting R code for the user to run.",
    "- When a tool returns an artifact_id, reuse it in subsequent tool calls to build on previous results.",
    "- Use the provided tools and ground your responses in their outputs. Cite specific numbers from tool results.",
    "- Prefer incremental steps over skipping ahead. Each tool call should produce a verifiable intermediate result.",
    "</tool_use_guidelines>",
    "",
    "<response_style>",
    "- Provide concise result summaries highlighting key quantitative findings.",
    "- After completing the requested task, briefly mention 1-2 possible next steps but do NOT execute them.",
    "- When interpreting results, focus on biological relevance, pathway crosstalk, and analytical limitations.",
    "- Detect the language of the user's input and respond in the same language.",
    "  For example, if the user writes in Chinese, respond in Chinese; if in English, respond in English.",
    "- Keep tool-generated identifiers (gene symbols, pathway IDs, compound names) in their original form.",
    "</response_style>",
    sep = "\n"
  )
}

#' Create a TCM analysis workflow
#'
#' Defines a multi-step analysis workflow that can be executed automatically.
#'
#' @param name Character. Workflow name.
#' @param steps List of step definitions, each with 'tool' and 'params'.
#'
#' @return A workflow object.
#'
#' @examples
#' \dontrun{
#'   workflow <- create_tcm_workflow(
#'     name = "herb_to_hub",
#'     steps = list(
#'       list(tool = "search_herb_records",
#'            params = list(herb = "{{herb}}", type = "Herb_pinyin_name")),
#'       list(tool = "run_herb_enrichment",
#'            params = list(genes = "{{from_previous}}")),
#'       list(tool = "interpret_artifact",
#'            params = list(artifact_id = "{{from_previous}}"))
#'     )
#'   )
#'   run_tcm_workflow(workflow, herb = "Astragalus")
#' }
#' @export
create_tcm_workflow <- function(name, steps) {
  structure(
    list(
      name  = name,
      steps = steps
    ),
    class = "tcm_workflow"
  )
}

#' Run a TCM analysis workflow
#'
#' Executes a predefined workflow step by step.
#'
#' @param workflow A workflow object from \code{\link{create_tcm_workflow}}.
#' @param ... Named parameters to substitute into step definitions.
#' @param verbose Logical. Print progress (default TRUE).
#'
#' @return A list of step results.
#' @export
run_tcm_workflow <- function(workflow, ..., verbose = TRUE) {
  params <- list(...)
  results <- list()
  previous_result <- NULL

  for (i in seq_along(workflow$steps)) {
    step <- workflow$steps[[i]]
    if (verbose) {
      message(sprintf("Step %d: %s", i, step$tool))
    }

    # Resolve parameters
    step_params <- step$params
    for (pname in names(step_params)) {
      pval <- step_params[[pname]]
      if (is.character(pval)) {
        # Substitute {{param_name}}
        for (key in names(params)) {
          pval <- gsub(sprintf("\\{\\{%s\\}\\}", key), params[[key]], pval)
        }
        # Substitute {{from_previous}}
        if (grepl("\\{\\{from_previous\\}\\}", pval) && !is.null(previous_result)) {
          # Try to get artifact_id or genes from previous result
          if (!is.null(previous_result$artifact_id)) {
            pval <- gsub("\\{\\{from_previous\\}\\}", previous_result$artifact_id, pval)
          }
        }
        step_params[[pname]] <- pval
      }
    }

    # Execute tool
    tools <- create_tcm_tools()
    tool <- Filter(function(t) t$name == step$tool, tools)
    if (length(tool) == 0) {
      stop(sprintf("Tool '%s' not found.", step$tool))
    }
    tool <- tool[[1]]

    result <- tryCatch({
      do.call(tool$execute, step_params)
    }, error = function(e) {
      list(ok = FALSE, error = conditionMessage(e))
    })

    results[[i]] <- result
    previous_result <- result

    if (!isTRUE(result$ok)) {
      warning(sprintf("Step %d failed: %s", i, result$error %||% "unknown error"))
      break
    }
  }

  results
}

#' Print method for TCM workflow
#' @param x A tcm_workflow object.
#' @param ... Additional arguments (ignored).
#' @export
print.tcm_workflow <- function(x, ...) {
  cat("TCM Workflow:", x$name, "\n")
  cat("Steps:\n")
  for (i in seq_along(x$steps)) {
    cat(sprintf("  %d. %s\n", i, x$steps[[i]]$tool))
  }
  invisible(x)
}

#' Create a generic aisdk agent
#'
#' A thin wrapper around \code{aisdk::create_agent()}.
#' Unlike \code{\link{create_tcm_task_agent}}, this creates a generic agent
#' without TCM-specific tools or system prompt.
#' Pass the result to \code{\link{run_tcm_agent}} to call it.
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
#' like any ordinary R function -- it works directly in \code{sapply()},
#' \code{purrr::map()}, \code{mutate()}, and other vectorised workflows.
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
#'   explain_pathway <- make_tcm_function(
#'     "You are an MSigDB/KEGG pathway expert. Explain the biological
#'      function of the pathway in 2 sentences for a research report."
#'   )
#'   explain_pathway("HALLMARK_INFLAMMATORY_RESPONSE")
#'
#'   pathways <- c("HALLMARK_E2F_TARGETS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION")
#'   explanations <- sapply(pathways, explain_pathway)
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
#' @param verbose Logical. If \code{TRUE} (default), print the response to the
#'   console. Set \code{FALSE} for programmatic / batch use.
#'
#' @return The response text (character). Returned invisibly when
#'   \code{verbose = TRUE}, visibly when \code{verbose = FALSE}.
#' @examples
#' \dontrun{
#'   bio_agent <- create_tcm_agent(
#'     name = "BioInterpreter",
#'     description = "Bioinformatics result interpreter",
#'     system_prompt = "You are a bioinformatics expert."
#'   )
#'   run_tcm_agent(bio_agent, "Explain HALLMARK_E2F_TARGETS")
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
