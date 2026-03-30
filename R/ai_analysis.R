# Public interpret_*() functions: enrichment, PPI, table.
# Each function: validate input → compress context → generate_object() → wrap.

# Context compression helpers 
#' Compress enrichment result to a prompt-ready string
#' @param x An enrichResult object or a data.frame with enrichment results.
#' @param top_n Integer. Number of top terms to include.
#' @param max_genes Integer. Maximum number of genes to include per term.
#' @keywords internal
#' @noRd
.compress_enrichment <- function(x, top_n = 10L, max_genes = 5L) {
  if (inherits(x, "enrichResult")) {
    df <- as.data.frame(x)
  } else if (is.data.frame(x)) {
    df <- x
  } else {
    stop("interpret_enrichment() expects an enrichResult or data.frame.", call. = FALSE)
  }

  # Standardise column names
  id_col   <- intersect(c("ID", "Description", "Term"), names(df))[1]
  padj_col <- intersect(c("p.adjust", "qvalue", "padj", "FDR"), names(df))[1]
  ratio_col <- intersect(c("GeneRatio", "gene_ratio", "Count"), names(df))[1]
  gene_col  <- intersect(c("geneID", "core_enrichment", "Genes"), names(df))[1]

  if (is.na(id_col)) stop("Cannot find term/ID column in enrichment result.", call. = FALSE)

  # Sort by significance ascending before selecting top_n rows, so that
  # "top_n" reliably means the most significant terms regardless of input order.
  if (!is.na(padj_col)) {
    df <- df[order(df[[padj_col]]), ]
  }

  n <- min(top_n, nrow(df))
  df <- utils::head(df, n)

  lines <- vapply(seq_len(n), function(i) {
    term <- df[[id_col]][i]
    padj <- if (!is.na(padj_col)) sprintf("%.2e", df[[padj_col]][i]) else "NA"
    ratio <- if (!is.na(ratio_col)) as.character(df[[ratio_col]][i]) else "NA"
    genes <- if (!is.na(gene_col)) {
      g <- strsplit(as.character(df[[gene_col]][i]), "/|,|;")[[1]]
      paste(utils::head(g, max_genes), collapse = ", ")
    } else {
      "NA"
    }
    sprintf("%d. %s | p.adjust=%s | ratio=%s | genes: %s", i, term, padj, ratio, genes)
  }, character(1))

  paste(lines, collapse = "\n")
}

#' Compress PPI node metrics to a prompt-ready string
#' @param x An igraph with vertex attributes, a node-metric data.frame, or the list output of rank_ppi_nodes().
#' @param top_n Integer. Number of top nodes to include.
#' @keywords internal
#' @noRd
.compress_ppi <- function(x, top_n = 10L) {
  if (inherits(x, "igraph")) {
    vattrs <- igraph::as_data_frame(x, what = "vertices")
  } else if (is.data.frame(x)) {
    vattrs <- x
  } else if (is.list(x) && length(x) >= 2 && is.data.frame(x[[2]])) {
    # rank_ppi_nodes() returns list(igraph, rank_df)
    vattrs <- x[[2]]
  } else {
    stop("interpret_ppi() expects an igraph, data.frame, or rank_ppi_nodes() output.",
         call. = FALSE)
  }

  metric_cols <- intersect(
    c(# Basic topology
      "degree", "strength",
      # Betweenness / closeness (unweighted + weighted)
      "betweenness", "betweenness_w",
      "closeness",   "closeness_w",
      # Centrality variants
      "eigen_centrality", "pagerank",
      "coreness", "clustering_coef",
      "eccentricity", "radiality",
      # Topology-based hub metrics (cytoHubba-style)
      "MCC", "MNC", "DMNC", "BN", "Stress", "EPC",
      # rank_ppi_nodes() composite score
      "Score_network", "Rank_network"),
    names(vattrs)
  )

  # Sort priority: Score_network (composite, higher = better) > degree.
  # rank_score kept for backward compat with any externally-built data.frames.
  sort_col <- intersect(c("Score_network", "rank_score", "degree"),
                        names(vattrs))[1]
  if (!is.na(sort_col)) {
    vattrs <- vattrs[order(vattrs[[sort_col]], decreasing = TRUE), ]
  }

  n <- min(top_n, nrow(vattrs))
  vattrs <- utils::head(vattrs, n)

  name_col <- intersect(c("name", "node", "gene", "protein"), names(vattrs))[1]
  if (is.na(name_col)) name_col <- 1

  lines <- vapply(seq_len(n), function(i) {
    node <- vattrs[[name_col]][i]
    metrics <- vapply(metric_cols, function(m) {
      sprintf("%s=%.3g", m, as.numeric(vattrs[[m]][i]))
    }, character(1))
    sprintf("%d. %s | %s", i, node, paste(metrics, collapse = " | "))
  }, character(1))

  header <- sprintf("PPI network: %d nodes shown (top %d by %s)",
                     nrow(vattrs), n, if (is.na(sort_col)) "order" else sort_col)
  paste(c(header, lines), collapse = "\n")
}

#' Compress a data.frame to a prompt-ready string
#'
#' Serialises the first \code{top_n} rows into a numbered key=value format
#' that any LLM can parse regardless of domain. Numeric values are formatted
#' with 4 significant figures; all other values are coerced to character.
#' @param x A data.frame.
#' @param top_n Integer. Number of top rows to include.
#' @keywords internal
#' @noRd
.compress_table <- function(x, top_n) {
  if (!is.data.frame(x)) {
    stop("interpret_table() expects a data.frame.", call. = FALSE)
  }

  n    <- min(top_n, nrow(x))
  x    <- utils::head(x, n)
  cols <- names(x)

  lines <- vapply(seq_len(n), function(i) {
    vals <- vapply(cols, function(col) {
      v <- x[[col]][i]
      if (is.numeric(v)) sprintf("%s=%.4g", col, v)
      else                sprintf("%s=%s",   col, as.character(v))
    }, character(1))
    sprintf("%d. %s", i, paste(vals, collapse = " | "))
  }, character(1))

  header <- sprintf("Table: %d rows x %d cols, showing top %d",
                    nrow(x), length(cols), n)
  paste(c(header, lines), collapse = "\n")
}

# System prompt builder

#' Build the system prompt for interpretation
#' @keywords internal
#' @noRd
.build_interpret_system <- function(type, audience, language, role = NULL) {
  # Role: who the AI is. Custom string replaces the default identity line.
  role_inst <- if (!is.null(role) && nzchar(role)) {
    role
  } else {
    sprintf("You are a bioinformatics expert interpreting %s results.", type)
  }

  lang_inst <- if (language == "zh") {
    "Respond entirely in Chinese (Simplified)."
  } else {
    "Respond entirely in English."
  }

  presets <- c(
    researcher = "Target audience: bioinformatics researchers. Use precise terminology.",
    wetlab     = "Target audience: wet-lab biologists. Explain terms clearly, avoid jargon.",
    paper      = "Target audience: journal manuscript. Use formal scientific writing style."
  )
  audience_inst <- if (audience %in% names(presets)) {
    presets[[audience]]
  } else {
    sprintf("Target audience: %s.", audience)
  }

  paste(
    role_inst,
    audience_inst,
    lang_inst,
    "Separate data-supported findings from hypotheses.",
    "Never invent statistics not present in the input.",
    "Be concise and precise.",
    "Respond with a single valid JSON object only.",
    sep = "\n"
  )
}


# Public API funciton
#' Interpret enrichment analysis results with AI
#' Uses a large language model to generate a structured interpretation of enrichment results (GO, KEGG, MSigDB, herb enrichment, etc.).
#'
#' @param x An \code{enrichResult} object or a data.frame with enrichment
#'   results.
#' @param top_n Integer. Number of top terms to include. Default 10.
#' @param max_genes Integer. Max genes shown per term in the compressed
#'   context. Default 5. Increase (e.g. \code{20}) when you want the model
#'   to reason about the full gene list of specific terms.
#' @param audience Character. Preset: \code{"researcher"},
#'   \code{"wetlab"}, \code{"paper"}. Or any free-text description
#'   like \code{"clinical doctor unfamiliar with bioinformatics"}.
#' @param language Character. \code{"zh"} or \code{"en"}.
#'   Default \code{"zh"}.
#' @param role Character or NULL. The AI's identity / expertise description.
#'   Replaces the default \emph{"You are a bioinformatics expert..."} opening
#'   line while keeping all other system-prompt constraints intact. Ignored
#'   when \code{system} is also provided (full override takes precedence).
#'   Example: \code{"You are a TCM pharmacologist specialising in
#'   herb-target interactions."}.
#' @param system Character or NULL. Full system prompt override. Replaces
#'   the entire generated system prompt when provided (including \code{role},
#'   \code{audience}, and \code{language} instructions).
#' @param prompt Character or NULL. Custom user prompt instruction.
#'   Replaces the default instruction; the compressed data context
#'   is always appended.
#' @param model A model identifier, LanguageModelV1 object, or NULL.
#'
#' @return A \code{tcm_ai_analysis} object.
#' @examples
#' \dontrun{
#'   enrich_res <- herb_enricher(genes = my_genes)
#'   interpret_enrichment(enrich_res)
#'
#'   # Custom role — keeps audience/language/constraints
#'   interpret_enrichment(enrich_res,
#'     role = "You are a TCM pharmacologist specialising in herb-target networks.")
#'
#'   # Show more genes per term
#'   interpret_enrichment(enrich_res, max_genes = 20)
#' }
#' @export
interpret_enrichment <- function(x,
                                 top_n = 10,
                                 max_genes = 5L,
                                 audience = "researcher",
                                 language = c("zh", "en"),
                                 role = NULL,
                                 system = NULL,
                                 prompt = NULL,
                                 model = NULL) {
  language <- match.arg(language)
  model <- .resolve_model(model)

  context <- .compress_enrichment(x, top_n, max_genes)
  system  <- system %||%
    .build_interpret_system("enrichment", audience, language, role)
  prompt  <- paste(
    prompt %||% "Interpret the following enrichment analysis results:",
    context, sep = "\n\n"
  )

  result <- .call_generate_object(model, prompt, .analysis_schema(), system,
                                   temperature = 0.3)

  input_class <- paste(class(x), collapse = "/")
  extracted   <- .extract_object(result, "analysis")
  .new_tcm_ai_analysis(
    input    = list(type = "enrichment", top_n = top_n,
                    max_genes = max_genes, object_class = input_class),
    context  = context,
    output   = extracted$output,
    metadata = .build_metadata(model, language, audience, input_class,
                               output_mode = extracted$output_mode)
  )
}

#' Interpret PPI network analysis results with AI
#'
#' @param x An igraph object with vertex attributes from
#'   \code{compute_nodeinfo()}, a node-metric data.frame, or the list
#'   output of \code{rank_ppi_nodes()}.
#' @param top_n Integer. Number of top nodes. Default 10.
#' @param audience Character. Default \code{"researcher"}.
#'   Also accepts free-text descriptions.
#' @param language Character. Default \code{"zh"}.
#' @param role Character or NULL. The AI's identity description. Replaces the
#'   default opening line; other constraints are preserved. Ignored when
#'   \code{system} is provided.
#' @param system Character or NULL. Full system prompt override.
#' @param prompt Character or NULL. Custom user prompt instruction.
#' @param model Model identifier or NULL.
#'
#' @return A \code{tcm_ai_analysis} object.
#' @examples
#' \dontrun{
#'   ppi <- compute_nodeinfo(demo_ppi)
#'   interpret_ppi(ppi)
#'   interpret_ppi(ppi, role = "You are a systems biology expert focusing on hub gene identification.")
#' }
#' @export
interpret_ppi <- function(x,
                          top_n = 10,
                          audience = "researcher",
                          language = c("zh", "en"),
                          role = NULL,
                          system = NULL,
                          prompt = NULL,
                          model = NULL) {
  language <- match.arg(language)
  model <- .resolve_model(model)

  context <- .compress_ppi(x, top_n)
  system <- system %||%
    .build_interpret_system("PPI network", audience, language, role)
  prompt <- paste(
    prompt %||% "Interpret the following PPI network centrality analysis:",
    context, sep = "\n\n"
  )

  result <- .call_generate_object(model, prompt, .analysis_schema(), system,
                                   temperature = 0.3)

  input_class <- paste(class(x), collapse = "/")
  extracted   <- .extract_object(result, "analysis")
  .new_tcm_ai_analysis(
    input    = list(type = "ppi", top_n = top_n, object_class = input_class),
    context  = context,
    output   = extracted$output,
    metadata = .build_metadata(model, language, audience, input_class,
                               output_mode = extracted$output_mode)
  )
}

#' Interpret a data.frame (table) with AI
#'
#' Accepts any data.frame and serialises it row-by-row into a compact
#' key=value format before sending to the model.
#'
#' @param x A data.frame.
#' @param top_n Integer. Number of top rows to include. Default 20.
#' @param audience Character. Default \code{"researcher"}.
#'   Also accepts free-text descriptions.
#' @param language Character. Default \code{"zh"}.
#' @param role Character or NULL. The AI's identity description. Replaces the
#'   default opening line; other constraints are preserved. Ignored when
#'   \code{system} is provided.
#' @param system Character or NULL. Full system prompt override.
#' @param prompt Character or NULL. Custom user prompt instruction.
#' @param model Model identifier or NULL.
#'
#' @return A \code{tcm_ai_analysis} object.
#' @examples
#' \dontrun{
#'   interpret_table(my_df)
#'   interpret_table(my_df, role = "You are a clinical data analyst.")
#' }
#' @export
interpret_table <- function(x,
                            top_n = 20,
                            audience = "researcher",
                            language = c("zh", "en"),
                            role = NULL,
                            system = NULL,
                            prompt = NULL,
                            model = NULL) {
  language <- match.arg(language)
  model <- .resolve_model(model)

  context <- .compress_table(x, top_n)
  system <- system %||%
    .build_interpret_system("tabular data", audience, language, role)
  prompt <- paste(
    prompt %||% "Interpret the following table:",
    context, sep = "\n\n"
  )

  result <- .call_generate_object(model, prompt, .analysis_schema(), system,
                                   temperature = 0.3)

  input_class <- paste(class(x), collapse = "/")
  extracted   <- .extract_object(result, "analysis")
  .new_tcm_ai_analysis(
    input    = list(type = "table", top_n = top_n,
                    object_class = input_class),
    context  = context,
    output   = extracted$output,
    metadata = .build_metadata(model, language, audience, input_class,
                               output_mode = extracted$output_mode)
  )
}

#' Interpret analysis results or free text with AI
#'
#' A unified entry point that dispatches to the appropriate
#' \code{interpret_*()} function based on the class of \code{x}.
#' When \code{x} is a character string (or vector), it routes to
#' an agent-based text interpreter instead.
#'
#' Supported input types:
#' \itemize{
#'   \item \code{character} — free-text interpretation
#'   \item \code{enrichResult} or enrichment data.frame
#'   \item \code{igraph} — PPI network
#'   \item \code{data.frame} — any table
#' }
#'
#' @param x Analysis result object, or a character string / vector
#'   for free-text interpretation.
#' @param type Character or NULL. Overrides auto-detection.
#'   One of \code{"enrichment"}, \code{"ppi"}, \code{"table"}.
#'   Takes priority over class-based dispatch, so a \code{data.frame}
#'   enrichment table can be forced with \code{type = "enrichment"}.
#' @param top_n Integer. Passed to the underlying interpret
#'   function (ignored for text).
#' @param max_genes Integer. Max genes shown per term when
#'   \code{type = "enrichment"}. Ignored for text, PPI, and generic tables.
#' @param audience Character. Preset: \code{"researcher"},
#'   \code{"wetlab"}, \code{"paper"}. Or any free-text description
#'   like \code{"clinical doctor, no bioinformatics background"}.
#' @param language Character. \code{"zh"} or \code{"en"}.
#'   Default \code{"zh"}.
#' @param role Character or NULL. The AI's identity description. Passed to
#'   the underlying \code{interpret_*()} function. Ignored for free-text
#'   input and when \code{system} is provided.
#' @param system Character or NULL. Full system prompt override.
#' @param prompt Character or NULL. Custom user prompt instruction.
#' @param model Model identifier or NULL.
#' @param verbose Logical. If \code{TRUE} (default), print the response to
#'   the console. Set \code{FALSE} for programmatic / batch use to suppress
#'   \code{cat()} output. Only affects the free-text path; structured results
#'   are always returned silently and printed via the S3 \code{print} method.
#'
#' @return For structured types: a \code{tcm_ai_analysis} object.
#'   For text: a character vector of responses (visibly when
#'   \code{verbose = FALSE}, invisibly otherwise).
#' @examples
#' \dontrun{
#'   tcm_interpret(enrich_res)
#'   tcm_interpret(ppi_graph,
#'     role = "You are a systems biologist focusing on drug target identification.")
#'   # Force enrichment path for a plain data.frame enrichment table
#'   tcm_interpret(enrich_df, type = "enrichment")
#'   tcm_interpret(my_df, type = "table")
#'   # Suppress printing for batch use
#'   results <- tcm_interpret("IL6 logFC=3.2, TNF logFC=2.8", verbose = FALSE)
#' }
#' @export
tcm_interpret <- function(x,
                          type = NULL,
                          top_n = NULL,
                          max_genes = 5L,
                          audience = "researcher",
                          language = c("zh", "en"),
                          role = NULL,
                          system = NULL,
                          prompt = NULL,
                          model = NULL,
                          verbose = TRUE) {
  language <- match.arg(language)

  if (is.character(x)) {
    return(
      .interpret_text(x, language, system, model, verbose = verbose,
                      audience = audience, role = role, prompt = prompt)
    )
  }

  detected <- if (!is.null(type)) {
    type
  } else if (inherits(x, "enrichResult")) {
    "enrichment"
  } else if (inherits(x, "igraph")) {
    "ppi"
  } else if (is.data.frame(x)) {
    "table"
  } else {
    stop(
      "Cannot auto-detect type for object of class '",
      paste(class(x), collapse = "/"), "'.",
      " Pass type = \"enrichment\"|\"ppi\"|\"table\".",
      call. = FALSE
    )
  }


  fn <- switch(detected,
    enrichment = interpret_enrichment,
    ppi        = interpret_ppi,
    table      = interpret_table,
    stop("Unknown type: ", detected, call. = FALSE)
  )

  default_top_n <- if (detected == "table") 20L else 10L
  top_n <- if (is.null(top_n)) default_top_n else as.integer(top_n)

  args <- list(
    x        = x,
    top_n    = top_n,
    audience = audience,
    language = language,
    role     = role,
    system   = system,
    prompt   = prompt,
    model    = model
  )
  if (detected == "enrichment") args$max_genes <- max_genes
  do.call(fn, args)
}

#' Interpret analysis results with a user-defined output schema
#'
#' A low-level entry point that keeps all existing context-compression and
#' prompt logic intact, but replaces the fixed output contract with a schema
#' supplied by the caller. Useful when you need a different set of fields than
#' the default \code{\link{tcm_interpret}} produces.
#'
#' The \code{schema} argument accepts any \code{aisdk::z_object()} constructed
#' from \code{aisdk::z_string()}, \code{aisdk::z_array()},
#' \code{aisdk::z_enum()}, etc. The model's JSON response is automatically
#' parsed into a named R list accessible via \code{result$output}.
#'
#' When the model does not honour structured output (e.g. older endpoints),
#' the function falls back to best-effort JSON extraction from the raw text;
#' if that also fails, \code{result$output} will be \code{list(raw_text = ...)}
#' so you always get a named list.
#'
#' @param x Analysis result object. Supported: \code{enrichResult},
#'   enrichment \code{data.frame}, \code{igraph}, any \code{data.frame}.
#' @param schema An \code{aisdk::z_object()} defining the output structure.
#' @param type Character or NULL. Override auto-detection.
#'   One of \code{"enrichment"}, \code{"ppi"}, \code{"table"}.
#' @param top_n Integer. Rows / nodes to include in the compressed context.
#'   Default 10 for enrichment/ppi, 20 for table.
#' @param max_genes Integer. Max genes per term (enrichment only). Default 5.
#' @param audience Character. Passed to the default system prompt builder.
#'   Ignored when \code{system} is provided.
#' @param language Character. \code{"zh"} or \code{"en"}. Default \code{"zh"}.
#' @param role Character or NULL. Replaces the AI identity line in the default
#'   system prompt. Ignored when \code{system} is provided.
#' @param system Character or NULL. Full system prompt override.
#' @param prompt Character or NULL. Instruction prepended to the context block.
#'   Defaults to \code{"Interpret the following <type> data:"}.
#' @param model Model identifier or NULL.
#'
#' @return A \code{tcm_ai_custom} object with fields:
#'   \code{$input}, \code{$context}, \code{$output} (the parsed schema),
#'   \code{$metadata}, \code{$schema}.
#'
#' @examples
#' \dontrun{
#'   # Define a custom schema
#'   my_schema <- aisdk::z_object(
#'     summary     = aisdk::z_string("Brief overview"),
#'     mechanism   = aisdk::z_string("Molecular mechanism"),
#'     key_targets = aisdk::z_array(aisdk::z_string(), "Top targets"),
#'     confidence  = aisdk::z_enum(c("high", "medium", "low"))
#'   )
#'
#'   res <- tcm_interpret_schema(
#'     enrich_res,
#'     schema   = my_schema,
#'     language = "zh",
#'     prompt   = "请重点关注抗炎通路："
#'   )
#'
#'   # Access fields by name
#'   res$output$summary
#'   res$output$key_targets   # character vector
#'   res$output$confidence
#'
#'   print(res)               # shows all fields generically
#' }
#' @export
tcm_interpret_schema <- function(x,
                                 schema,
                                 type      = NULL,
                                 top_n     = NULL,
                                 max_genes = 5L,
                                 audience  = "researcher",
                                 language  = c("zh", "en"),
                                 role      = NULL,
                                 system    = NULL,
                                 prompt    = NULL,
                                 model     = NULL) {
  .check_aisdk()
  language <- match.arg(language)
  model    <- .resolve_model(model)

  # Type detection — identical logic to tcm_interpret()
  detected <- if (!is.null(type)) {
    type
  } else if (inherits(x, "enrichResult")) {
    "enrichment"
  } else if (inherits(x, "igraph")) {
    "ppi"
  } else if (is.data.frame(x)) {
    "table"
  } else {
    stop(
      "Cannot auto-detect type for object of class '",
      paste(class(x), collapse = "/"), "'.",
      " Pass type = \"enrichment\"|\"ppi\"|\"table\".",
      call. = FALSE
    )
  }

  # Compress context via existing helpers
  default_top_n <- if (detected == "table") 20L else 10L
  top_n <- as.integer(top_n %||% default_top_n)

  context <- switch(detected,
    enrichment = .compress_enrichment(x, top_n, max_genes),
    ppi        = .compress_ppi(x, top_n),
    table      = .compress_table(x, top_n)
  )

  # Build system prompt via existing helper (role + audience + language)
  system <- system %||% .build_interpret_system(detected, audience, language, role)

  user_prompt <- paste(
    prompt %||% sprintf("Interpret the following %s data:", detected),
    context,
    sep = "\n\n"
  )

  result <- .call_generate_object(model, user_prompt, schema, system)

  extracted <- .extract_custom_object(result)
  .new_tcm_ai_custom(
    input    = list(type         = detected,
                    object_class = paste(class(x), collapse = "/")),
    context  = context,
    output   = extracted$output,
    metadata = .build_metadata(model, language, audience,
                               paste(class(x), collapse = "/"),
                               output_mode = extracted$output_mode),
    schema   = schema
  )
}

# Agent-based free-text interpretation (internal)
# @keywords internal
# @noRd
.interpret_text <- function(x, language, system, model, verbose = TRUE,
                            audience = "researcher", role = NULL,
                            prompt = NULL) {
  .check_aisdk()
  model <- .resolve_model(model)

  role_inst <- if (!is.null(role) && nzchar(role)) {
    role
  } else {
    "You are a bioinformatics interpretation assistant."
  }

  audience_presets <- c(
    researcher = paste("Target audience: bioinformatics researchers.",
                       "Use precise terminology."),
    wetlab     = paste("Target audience: wet-lab biologists.",
                       "Explain terms clearly; avoid jargon."),
    paper      = paste("Target audience: journal manuscript.",
                       "Use formal scientific writing style.")
  )
  audience_inst <- if (audience %in% names(audience_presets)) {
    audience_presets[[audience]]
  } else {
    sprintf("Target audience: %s.", audience)
  }

  lang_inst <- if (language == "zh") {
    "Respond entirely in Chinese (Simplified)."
  } else {
    "Respond entirely in English."
  }

  system <- system %||% paste(
    role_inst,
    audience_inst,
    lang_inst,
    "Explain results clearly and concisely.",
    "Do not invent statistics not in the input.",
    sep = "\n"
  )

  agent <- create_tcm_agent(
    name          = "TextExplainer",
    description   = "Bioinformatics text interpreter",
    system_prompt = system
  )

  # Prepend custom prompt instruction to each query if provided
  make_task <- if (!is.null(prompt) && nzchar(prompt)) {
    function(q) paste(prompt, q, sep = "\n\n")
  } else {
    function(q) q
  }

  out <- vapply(x, function(query) {
    res <- agent$run(task = make_task(query), model = model)
    txt <- res$text
    if (verbose) {
      if (length(x) > 1L) {
        cat(sprintf("--- [%s] ---\n", substr(query, 1, 50)))
      }
      cat(txt, "\n\n")
    }
    txt
  }, character(1), USE.NAMES = FALSE)

  if (verbose) invisible(out) else out
}
