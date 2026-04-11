# ai_tools.R
# Wrapper tools for the TCM agent layer.
# Wraps deterministic TCMDATA and clusterProfiler functions as aisdk tools.

#' Create all TCM analysis tools
#'
#' Returns a list of aisdk Tool objects that wrap TCMDATA analysis functions.
#' These tools can be passed to an aisdk agent for native tool calling.
#'
#' @param task_type Character or NULL. If provided, only returns tools relevant
#'   to this task type (for example, "enrichment", "ppi_analysis",
#'   or "ml_screening"). If NULL, returns all tools.
#' @param tool_names Character vector or NULL. If provided, overrides task_type
#'   and returns only tools whose names are in this vector. Useful when the
#'   router merges tools from multiple matched task types.
#'
#' @return A list of aisdk Tool objects.
#'
#' @examples
#' \dontrun{
#'   tools <- create_tcm_tools()
#'   tools <- create_tcm_tools(task_type = "enrichment")
#' }
#' @export
create_tcm_tools <- function(task_type = NULL, tool_names = NULL) {
  .check_aisdk()

  all_tools <- list(
    # --- Search Tools ---
    tool_search_herb(),
    tool_search_target(),
    tool_search_disease(),
    tool_search_gene_disease(),

    # --- Intersection Tool ---
    tool_compute_target_intersection(),

    # --- Enrichment Tools ---
    tool_run_herb_enrichment(),
    tool_run_go_enrichment(),
    tool_run_kegg_enrichment(),

    # --- PPI Tools ---
    tool_get_ppi_network(),
    tool_subset_ppi_network(),
    tool_compute_ppi_metrics(),
    tool_rank_ppi_nodes(),
    tool_run_mcode_clustering(),

    # --- Machine Learning Tools ---
    tool_prepare_ml_dataset(),
    tool_run_ml_screening(),
    tool_get_ml_consensus(),

    # --- Literature Tools ---
    tool_get_pubmed_evidence(),
    tool_extract_pubmed_table(),

    # --- GEO Search Tools ---
    tool_search_geo_datasets(),

    # --- Compound Tools ---
    tool_resolve_compound_cid(),
    tool_get_compound_properties(),
    tool_run_compound_similarity(),

    # --- Visualization Tools ---
    tool_plot_enrichment_result(),
    tool_plot_ppi_result(),
    tool_plot_ml_result(),
    tool_plot_pubmed_result(),
    tool_plot_herb_sankey(),
    tool_plot_docking_heatmap(),

    # --- Artifact Management Tools ---
    tool_list_artifacts(),
    tool_load_artifact_summary(),

    # --- Environment Access Tools ---
    tool_read_variable(),
    tool_eval_r_code(),

    # --- Cross-Validation Tools ---
    tool_generate_verification_urls(),

    # --- Interpretation Tool ---
    tool_interpret_artifact()
  )

  if (!is.null(task_type)) {
    tool_map <- list(
      herb_lookup = c(
        "search_herb_records",
        "search_targets",
        "plot_herb_sankey",
        "generate_verification_urls"
      ),
      disease_lookup = c(
        "search_disease_targets",
        "search_gene_diseases",
        "search_geo_datasets",
        "compute_target_intersection",
        "run_go_enrichment",
        "run_kegg_enrichment",
        "get_ppi_network",
        "interpret_artifact"
      ),
      enrichment = c(
        "run_herb_enrichment",
        "run_go_enrichment",
        "run_kegg_enrichment",
        "plot_enrichment_result",
        "interpret_artifact"
      ),
      ppi_analysis = c(
        "get_ppi_network",
        "subset_ppi_network",
        "compute_ppi_metrics",
        "rank_ppi_nodes",
        "run_mcode_clustering",
        "plot_ppi_result",
        "interpret_artifact"
      ),
      ml_screening = c(
        "prepare_ml_dataset",
        "run_ml_screening",
        "get_ml_consensus",
        "plot_ml_result",
        "interpret_artifact"
      ),
      pubmed = c(
        "get_pubmed_evidence",
        "extract_pubmed_table",
        "plot_pubmed_result",
        "interpret_artifact"
      ),
      geo_search = c(
        "search_geo_datasets",
        "interpret_artifact"
      ),
      compound = c(
        "resolve_compound_cid",
        "get_compound_properties",
        "run_compound_similarity",
        "plot_docking_heatmap",
        "interpret_artifact"
      ),
      visualization = c(
        "plot_enrichment_result",
        "plot_ppi_result",
        "plot_ml_result",
        "plot_pubmed_result",
        "plot_herb_sankey",
        "plot_docking_heatmap"
      ),
      interpretation = c(
        "list_artifacts",
        "load_artifact_summary",
        "interpret_artifact",
        "generate_verification_urls"
      )
    )

    allowed <- tool_map[[task_type]]
    if (!is.null(allowed)) {
      all_tools <- Filter(function(tool) tool$name %in% allowed, all_tools)
    }
  }

  # Direct tool name filtering (overrides task_type when provided)
  if (!is.null(tool_names)) {
    all_tools <- Filter(function(tool) tool$name %in% tool_names, all_tools)
  }

  all_tools
}

#' @keywords internal
#' @noRd
.require_namespace_for_tool <- function(package_name, feature_name = package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop(sprintf(
      "Package '%s' is required for %s(). Please install it first.",
      package_name,
      feature_name
    ), call. = FALSE)
  }
}

#' @keywords internal
#' @noRd
.clean_character_vector <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x[nzchar(x)]
}

#' @keywords internal
#' @noRd
.empty_string_to_null <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NULL)
  }
  x <- as.character(x)
  x <- trimws(x)
  if (all(!nzchar(x))) {
    return(NULL)
  }
  x[nzchar(x)]
}

#' @keywords internal
#' @noRd
.normalize_optional_integer <- function(x, default = NULL) {
  if (is.null(x) || length(x) == 0 || is.na(x) || x <= 0) {
    return(default)
  }
  as.integer(x)
}

#' @keywords internal
#' @noRd
.load_orgdb_object <- function(orgdb_package) {
  orgdb_package <- .empty_string_to_null(orgdb_package)
  if (is.null(orgdb_package)) {
    stop("An OrgDb package name is required, for example 'org.Hs.eg.db'.",
         call. = FALSE)
  }
  .require_namespace_for_tool(orgdb_package, "OrgDb lookup")
  get(orgdb_package, envir = asNamespace(orgdb_package))
}

#' @keywords internal
#' @noRd
.save_tool_artifact <- function(object,
                                artifact_type,
                                function_name,
                                params = list(),
                                summary = NULL,
                                next_actions = character(),
                                extra = list()) {
  # Quality gate: auto-diagnose common issues and append warnings
  warnings <- .quality_gate(object, artifact_type)

  if (length(warnings) > 0 && !is.null(summary)) {
    summary <- paste0(summary, " WARNING: ", paste(warnings, collapse = "; "))
  }

  artifact <- save_tcm_artifact(
    object = object,
    artifact_type = artifact_type,
    summary = summary,
    provenance = list(
      function_name = function_name,
      params = params
    )
  )

  base <- list(
    ok = TRUE,
    artifact_id = artifact$artifact_id,
    artifact_type = artifact$artifact_type,
    summary = artifact$summary,
    preview = artifact$preview
  )

  if (length(warnings) > 0) {
    base$quality_warnings <- warnings
  }

  if (length(next_actions) > 0) {
    base$next_actions <- next_actions
  }

  if (length(extra) > 0) {
    base <- utils::modifyList(base, extra)
  }

  base
}

#' @keywords internal
#' @noRd
.quality_gate <- function(object, artifact_type) {
  warnings <- character(0)

  # -- Search results --
  if (artifact_type == "search_result") {
    if (is.data.frame(object) && nrow(object) == 0) {
      warnings <- c(warnings, "Search returned 0 records. Try broader search terms or check name spelling.")
    }
  }

  # -- Intersection results --
  if (artifact_type == "intersection_result") {
    n_genes <- 0L
    if (is.list(object) && !is.null(object$intersection)) {
      n_genes <- length(object$intersection)
    }
    if (n_genes == 0) {
      warnings <- c(warnings, "Intersection is empty. The gene sets have no overlap. Try different search terms or data sources.")
    } else if (n_genes < 5) {
      warnings <- c(warnings, sprintf("Only %d intersection genes. Results may be underpowered. Consider broadening search criteria.", n_genes))
    } else if (n_genes > 500) {
      warnings <- c(warnings, sprintf("%d intersection genes is very large. Consider narrowing search criteria for more specific results.", n_genes))
    }
  }

  # -- Enrichment results --
  if (artifact_type == "enrichment_result") {
    df <- tryCatch(as.data.frame(object), error = function(e) data.frame())
    if (nrow(df) == 0) {
      warnings <- c(warnings, "No enrichment terms found. Try relaxing pvalueCutoff/qvalueCutoff or using a different gene set.")
    } else {
      sig_count <- sum(df$p.adjust < 0.05, na.rm = TRUE)
      if (sig_count == 0) {
        warnings <- c(warnings, "No terms reached significance (p.adjust < 0.05). Consider relaxing cutoffs or checking gene list quality.")
      }
    }
  }

  # -- PPI graph --
  if (artifact_type == "ppi_graph" && inherits(object, "igraph")) {
    n_nodes <- igraph::vcount(object)
    n_edges <- igraph::ecount(object)
    if (n_edges == 0) {
      warnings <- c(warnings, "PPI network has 0 edges. The genes may not have known interactions at this confidence level. Try lowering score_threshold.")
    } else {
      n_components <- igraph::count_components(object)
      isolated <- sum(igraph::degree(object) == 0)
      if (isolated > n_nodes * 0.5) {
        warnings <- c(warnings, sprintf("%d of %d nodes are isolated (no interactions). Consider lowering score_threshold for better connectivity.", isolated, n_nodes))
      }
      if (n_components > 1 && n_components > n_nodes * 0.3) {
        warnings <- c(warnings, sprintf("Network is fragmented into %d components. Hub gene analysis may be less meaningful.", n_components))
      }
    }
  }

  # -- ML results --
  if (artifact_type == "ml_result" && inherits(object, "tcm_ml_list")) {
    methods_run <- length(object)
    if (methods_run < 2) {
      warnings <- c(warnings, "Only 1 ML method was run. Consensus analysis requires at least 2 methods for reliable feature selection.")
    }
  }

  warnings
}

#' @keywords internal
#' @noRd
.build_expression_matrix <- function(expr_matrix, gene_names, sample_names = character(0)) {
  if (length(expr_matrix) == 0) {
    stop("expr_matrix must contain at least one gene row.", call. = FALSE)
  }

  gene_names <- .clean_character_vector(gene_names)
  if (length(gene_names) == 0) {
    stop("gene_names must contain at least one gene symbol.", call. = FALSE)
  }

  row_lengths <- lengths(expr_matrix)
  if (length(unique(row_lengths)) != 1L) {
    stop("Each row in expr_matrix must have the same number of sample values.",
         call. = FALSE)
  }

  mat <- do.call(rbind, lapply(expr_matrix, as.numeric))
  if (nrow(mat) != length(gene_names)) {
    stop("The number of gene_names must match the number of rows in expr_matrix.",
         call. = FALSE)
  }

  sample_names <- .empty_string_to_null(sample_names)
  if (is.null(sample_names)) {
    sample_names <- sprintf("Sample_%03d", seq_len(ncol(mat)))
  }
  if (length(sample_names) != ncol(mat)) {
    stop("sample_names must match the number of columns in expr_matrix.",
         call. = FALSE)
  }

  rownames(mat) <- gene_names
  colnames(mat) <- sample_names
  storage.mode(mat) <- "numeric"
  mat
}

#' @keywords internal
#' @noRd
.standardize_ppi_graph <- function(graph) {
  if (!inherits(graph, "igraph")) {
    stop("The returned PPI object is not an igraph object.", call. = FALSE)
  }

  edge_attrs <- igraph::edge_attr_names(graph)
  candidate_attrs <- c("score", "combined_score", "combinedScore", "weight")
  available <- intersect(candidate_attrs, edge_attrs)

  if (length(available) == 0) {
    return(graph)
  }

  score_values <- igraph::edge_attr(graph, available[1])
  if (is.numeric(score_values) && any(score_values > 1, na.rm = TRUE) &&
      max(score_values, na.rm = TRUE) <= 1000) {
    score_values <- score_values / 1000
  }

  graph <- igraph::set_edge_attr(graph, "score", value = score_values)
  graph
}

#' @keywords internal
#' @noRd
.summarize_enrichment <- function(result, label) {
  df <- as.data.frame(result)
  sig_count <- sum(df$p.adjust < 0.05, na.rm = TRUE)
  sprintf(
    "%s: %d total terms, %d significant (p.adjust < 0.05).",
    label,
    nrow(df),
    sig_count
  )
}

#' @keywords internal
#' @noRd
.extract_cids_from_artifact <- function(artifact_id) {
  artifact_id <- .empty_string_to_null(artifact_id)
  if (is.null(artifact_id)) {
    return(character(0))
  }

  object <- load_tcm_artifact(artifact_id)
  if (!is.data.frame(object) || !"cid" %in% names(object)) {
    stop("The artifact does not contain a 'cid' column.", call. = FALSE)
  }

  .clean_character_vector(object$cid)
}

#' @keywords internal
#' @noRd
.ml_result_summary_table <- function(ml_list) {
  if (!inherits(ml_list, "tcm_ml_list")) {
    stop("Expected a tcm_ml_list object.", call. = FALSE)
  }

  rows <- lapply(names(ml_list), function(method_name) {
    model <- ml_list[[method_name]]
    perf <- model$test_performance %||% model$cv_performance %||% list()
    data.frame(
      method = method_name,
      features = length(model$selected_features %||% character(0)),
      auc = as.numeric(perf$auc %||% NA_real_),
      sensitivity = as.numeric(perf$sensitivity %||% NA_real_),
      specificity = as.numeric(perf$specificity %||% NA_real_),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

#' @keywords internal
#' @noRd
.as_rank_table <- function(object) {
  if (is.data.frame(object)) {
    return(object)
  }
  if (is.list(object) && !is.null(object$table) && is.data.frame(object$table)) {
    return(object$table)
  }
  if (inherits(object, "igraph")) {
    return(rank_ppi_nodes(object)$table)
  }
  stop("Expected a ranked PPI table, a rank_ppi_nodes() result, or an igraph object.",
       call. = FALSE)
}

#' @keywords internal
#' @noRd
.coerce_docking_matrix <- function(object) {
  if (is.matrix(object)) {
    return(object)
  }
  if (!is.data.frame(object)) {
    stop("The docking artifact must be a matrix or data.frame.", call. = FALSE)
  }

  if ("target" %in% names(object)) {
    rownames(object) <- object$target
    object$target <- NULL
  }

  object
}

#' Tool: Search herb records
#' @keywords internal
#' @noRd
tool_search_herb <- function() {
  aisdk::tool(
    name = "search_herb_records",
    description = paste(
      "Search herb records, including linked molecules and targets.",
      "Use this when the user asks about herb composition, compounds, or targets."
    ),
    parameters = aisdk::z_object(
      herb = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Vector of herb names to search"
      ),
      type = aisdk::z_enum(
        values = c("Herb_cn_name", "Herb_pinyin_name", "Herb_en_name"),
        description = "Herb name type"
      )
    ),
    execute = function(herb, type) {
      tryCatch({
        herb <- .clean_character_vector(herb)
        result <- search_herb(herb = herb, type = type)

        .save_tool_artifact(
          object = result,
          artifact_type = "search_result",
          function_name = "search_herb",
          params = list(herb = herb, type = type),
          summary = sprintf(
            "Found %d herb-molecule-target records for %d herb(s).",
            nrow(result),
            length(herb)
          ),
          next_actions = list(
            list(tool = "run_herb_enrichment",
                 reason = "Identify which herbs are significantly associated with these targets"),
            list(tool = "run_go_enrichment",
                 reason = "Perform GO functional annotation on the retrieved targets"),
            list(tool = "plot_herb_sankey",
                 reason = "Visualize herb-compound-target relationships")
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Search target records
#' @keywords internal
#' @noRd
tool_search_target <- function() {
  aisdk::tool(
    name = "search_targets",
    description = paste(
      "Search target-centric herb records.",
      "Use this when the user asks which herbs or compounds are associated with specific genes or proteins."
    ),
    parameters = aisdk::z_object(
      target = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Vector of target gene symbols"
      )
    ),
    execute = function(target) {
      tryCatch({
        target <- .clean_character_vector(target)
        result <- search_target(gene_list = target)

        .save_tool_artifact(
          object = result,
          artifact_type = "search_result",
          function_name = "search_target",
          params = list(target = target),
          summary = sprintf(
            "Found %d herb-molecule-target records for %d target(s).",
            nrow(result),
            length(target)
          ),
          next_actions = c(
            "run_go_enrichment",
            "run_kegg_enrichment",
            "get_ppi_network"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Search disease targets (disease -> genes)
#' @keywords internal
#' @noRd
tool_search_disease <- function() {
  aisdk::tool(
    name = "search_disease_targets",
    description = paste(
      "Search for genes associated with a disease using DisGeNET data.",
      "Accepts disease names (e.g. 'sepsis', 'diabetes') or UMLS CUI IDs (e.g. 'C0243026').",
      "Returns disease_id, disease_name, gene_id, and gene symbol.",
      "Data source: DisGeNET via DOSE package (30170 diseases, 21671 genes)."
    ),
    parameters = aisdk::z_object(
      disease = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Disease name(s) or UMLS CUI ID(s)"
      )
    ),
    execute = function(disease) {
      tryCatch({
        disease <- .clean_character_vector(disease)
        result <- search_disease(disease = disease, readable = TRUE)
        if (is.null(result) || nrow(result) == 0) {
          return(list(ok = FALSE, error = paste("No disease found for:", paste(disease, collapse = ", "))))
        }

        .save_tool_artifact(
          object = result,
          artifact_type = "search_result",
          function_name = "search_disease",
          params = list(disease = disease),
          summary = sprintf(
            "Found %d gene(s) for %d disease(s): %s.",
            length(unique(result$gene_id)),
            length(unique(result$disease_id)),
            paste(unique(result$disease_name), collapse = ", ")
          ),
          next_actions = list(
            list(tool = "run_go_enrichment",
                 reason = "Annotate disease gene functions via GO enrichment"),
            list(tool = "run_kegg_enrichment",
                 reason = "Identify disease-related KEGG pathways"),
            list(tool = "get_ppi_network",
                 reason = "Build PPI network for disease gene interaction analysis"),
            list(tool = "search_gene_disease",
                 reason = "Reverse lookup to find additional disease associations")
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Search gene-associated diseases (gene -> diseases)
#' @keywords internal
#' @noRd
tool_search_gene_disease <- function() {
  aisdk::tool(
    name = "search_gene_diseases",
    description = paste(
      "Reverse lookup: find diseases associated with given genes.",
      "Accepts gene symbols (e.g. 'TNF', 'IL6') or Entrez IDs.",
      "Returns disease_id, disease_name, gene_id, and gene symbol.",
      "Data source: DisGeNET via DOSE package."
    ),
    parameters = aisdk::z_object(
      gene = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Gene symbol(s) or Entrez ID(s)"
      )
    ),
    execute = function(gene) {
      tryCatch({
        gene <- .clean_character_vector(gene)
        result <- search_gene_disease(gene = gene, readable = TRUE)
        if (is.null(result) || nrow(result) == 0) {
          return(list(ok = FALSE, error = paste("No diseases found for:", paste(gene, collapse = ", "))))
        }

        .save_tool_artifact(
          object = result,
          artifact_type = "search_result",
          function_name = "search_gene_disease",
          params = list(gene = gene),
          summary = sprintf(
            "Found %d disease(s) for %d gene(s): %s.",
            length(unique(result$disease_id)),
            length(unique(result$gene_id)),
            paste(unique(result$symbol[!is.na(result$symbol)]), collapse = ", ")
          ),
          next_actions = c(
            "search_disease_targets", "run_go_enrichment",
            "get_ppi_network"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Compute target intersection (herb targets \u2229 disease targets)
#' @keywords internal
#' @noRd
tool_compute_target_intersection <- function() {
  aisdk::tool(
    name = "compute_target_intersection",
    description = paste(
      "Compute the intersection of two or more gene/target lists.",
      "Typically used to find overlapping targets between herb targets and disease targets.",
      "Pass artifact_ids from previous search results, or raw gene vectors.",
      "Returns the intersection genes and a Venn diagram data frame."
    ),
    parameters = aisdk::z_object(
      gene_lists = aisdk::z_array(
        items = aisdk::z_array(items = aisdk::z_string()),
        description = "Two or more gene symbol vectors to intersect"
      ),
      set_names = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Names for each gene set (e.g. 'Herb Targets', 'Disease Targets')"
      )
    ),
    execute = function(gene_lists, set_names = NULL) {
      tryCatch({
        gene_lists <- lapply(gene_lists, .clean_character_vector)
        if (length(gene_lists) < 2 || length(gene_lists) > 4) {
          return(list(ok = FALSE, error = "Provide 2-4 gene lists."))
        }
        if (is.null(set_names)) {
          set_names <- paste0("Set", seq_along(gene_lists))
        } else {
          set_names <- .clean_character_vector(set_names)
        }

        # Build Venn data
        venn_df <- do.call(getvenndata, c(gene_lists, list(set_names = set_names)))

        # Compute intersection (all TRUE)
        logical_cols <- venn_df[, -1, drop = FALSE]
        all_true <- apply(logical_cols, 1, all)
        intersection_genes <- venn_df$Element[all_true]

        # Per-set counts
        set_counts <- vapply(set_names, function(nm) sum(venn_df[[nm]]), integer(1))
        count_str <- paste(sprintf("%s=%d", names(set_counts), set_counts), collapse = ", ")

        result <- list(
          intersection = intersection_genes,
          venn_df = venn_df,
          set_sizes = set_counts
        )

        .save_tool_artifact(
          object = result,
          artifact_type = "intersection_result",
          function_name = "compute_target_intersection",
          params = list(set_names = set_names),
          summary = sprintf(
            "Intersection of %d sets (%s): %d common targets. Sets: %s.",
            length(gene_lists), paste(set_names, collapse = " & "),
            length(intersection_genes), count_str
          ),
          next_actions = list(
            list(tool = "get_ppi_network",
                 reason = "Build PPI network from intersection genes to find hub genes"),
            list(tool = "run_go_enrichment",
                 reason = "Annotate intersection gene functions"),
            list(tool = "run_kegg_enrichment",
                 reason = "Identify enriched pathways in intersection genes"),
            list(tool = "interpret_artifact",
                 reason = "Generate biological interpretation of the intersection")
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Run herb enrichment analysis
#' @keywords internal
#' @noRd
tool_run_herb_enrichment <- function() {
  aisdk::tool(
    name = "run_herb_enrichment",
    description = paste(
      "Perform herb-target enrichment analysis with clusterProfiler.",
      "Use this to identify herbs significantly associated with a gene list."
    ),
    parameters = aisdk::z_object(
      genes = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Vector of gene symbols"
      ),
      type = aisdk::z_enum(
        values = c("Herb_pinyin_name", "Herb_cn_name", "Herb_en_name"),
        description = "Herb label column used in enrichment terms"
      ),
      pvalueCutoff = aisdk::z_number(
        description = "P-value cutoff"
      ),
      qvalueCutoff = aisdk::z_number(
        description = "Q-value cutoff"
      )
    ),
    execute = function(genes,
                       type = "Herb_pinyin_name",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2) {
      tryCatch({
        genes <- unique(.clean_character_vector(genes))
        result <- herb_enricher(
          genes = genes,
          type = type,
          pvalueCutoff = pvalueCutoff,
          qvalueCutoff = qvalueCutoff
        )

        .save_tool_artifact(
          object = result,
          artifact_type = "enrichment_result",
          function_name = "herb_enricher",
          params = list(
            genes = genes,
            type = type,
            pvalueCutoff = pvalueCutoff,
            qvalueCutoff = qvalueCutoff
          ),
          summary = .summarize_enrichment(
            result,
            sprintf("Herb enrichment on %d gene(s)", length(genes))
          ),
          next_actions = list(
            list(tool = "plot_enrichment_result",
                 reason = "Visualize herb enrichment results as lollipop plot"),
            list(tool = "interpret_artifact",
                 reason = "AI interpretation of enriched herbs and their significance"),
            list(tool = "get_ppi_network",
                 reason = "Build PPI network for the input genes")
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Run GO enrichment with clusterProfiler
#' @keywords internal
#' @noRd
tool_run_go_enrichment <- function() {
  aisdk::tool(
    name = "run_go_enrichment",
    description = paste(
      "Run GO enrichment with clusterProfiler::enrichGO().",
      "Use this for Biological Process, Cellular Component, Molecular Function, or ALL GO categories."
    ),
    parameters = aisdk::z_object(
      genes = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Vector of gene identifiers"
      ),
      orgdb_package = aisdk::z_string(
        description = "OrgDb package name, for example org.Hs.eg.db"
      ),
      key_type = aisdk::z_string(
        description = "Gene identifier type, for example SYMBOL or ENTREZID"
      ),
      ont = aisdk::z_enum(
        values = c("BP", "CC", "MF", "ALL"),
        description = "GO ontology branch"
      ),
      readable = aisdk::z_boolean(
        description = "Whether to convert gene IDs to readable symbols"
      ),
      pvalueCutoff = aisdk::z_number(
        description = "P-value cutoff"
      ),
      qvalueCutoff = aisdk::z_number(
        description = "Q-value cutoff"
      )
    ),
    execute = function(genes,
                       orgdb_package = "org.Hs.eg.db",
                       key_type = "SYMBOL",
                       ont = "BP",
                       readable = TRUE,
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2) {
      tryCatch({
        .require_namespace_for_tool("clusterProfiler", "run_go_enrichment")
        orgdb <- .load_orgdb_object(orgdb_package)
        genes <- unique(.clean_character_vector(genes))

        result <- clusterProfiler::enrichGO(
          gene = genes,
          OrgDb = orgdb,
          ont = ont,
          keyType = key_type,
          readable = isTRUE(readable),
          pvalueCutoff = pvalueCutoff,
          qvalueCutoff = qvalueCutoff
        )

        .save_tool_artifact(
          object = result,
          artifact_type = "enrichment_result",
          function_name = "clusterProfiler::enrichGO",
          params = list(
            genes = genes,
            orgdb_package = orgdb_package,
            key_type = key_type,
            ont = ont,
            readable = readable,
            pvalueCutoff = pvalueCutoff,
            qvalueCutoff = qvalueCutoff
          ),
          summary = .summarize_enrichment(
            result,
            sprintf("GO enrichment (%s) on %d gene(s)", ont, length(genes))
          ),
          next_actions = list(
            list(tool = "plot_enrichment_result",
                 reason = "Visualize GO enrichment results"),
            list(tool = "interpret_artifact",
                 reason = "AI interpretation of enriched GO terms and biological processes")
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Run KEGG enrichment with clusterProfiler
#' @keywords internal
#' @noRd
tool_run_kegg_enrichment <- function() {
  aisdk::tool(
    name = "run_kegg_enrichment",
    description = paste(
      "Run KEGG enrichment with clusterProfiler::enrichKEGG().",
      "Use this for pathway-level over-representation analysis."
    ),
    parameters = aisdk::z_object(
      genes = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Vector of gene identifiers"
      ),
      organism = aisdk::z_string(
        description = "KEGG organism code, for example hsa or mmu"
      ),
      key_type = aisdk::z_enum(
        values = c("kegg", "ncbi-geneid", "ncbi-proteinid", "uniprot"),
        description = "Identifier type used by enrichKEGG()"
      ),
      readable = aisdk::z_boolean(
        description = "Whether to attempt readable gene symbols when possible"
      ),
      orgdb_package = aisdk::z_string(
        description = "Optional OrgDb package name used by setReadable(); pass an empty string to skip"
      ),
      pvalueCutoff = aisdk::z_number(
        description = "P-value cutoff"
      ),
      qvalueCutoff = aisdk::z_number(
        description = "Q-value cutoff"
      )
    ),
    execute = function(genes,
                       organism = "hsa",
                       key_type = "kegg",
                       readable = TRUE,
                       orgdb_package = "",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.2) {
      tryCatch({
        .require_namespace_for_tool("clusterProfiler", "run_kegg_enrichment")
        genes <- unique(.clean_character_vector(genes))

        result <- clusterProfiler::enrichKEGG(
          gene = genes,
          organism = organism,
          keyType = key_type,
          pvalueCutoff = pvalueCutoff,
          qvalueCutoff = qvalueCutoff
        )

        if (isTRUE(readable) && key_type %in% c("kegg", "ncbi-geneid")) {
          orgdb_package <- .empty_string_to_null(orgdb_package)
          if (!is.null(orgdb_package)) {
            orgdb <- .load_orgdb_object(orgdb_package)
            result <- tryCatch(
              clusterProfiler::setReadable(
                result,
                OrgDb = orgdb,
                keyType = "ENTREZID"
              ),
              error = function(e) result
            )
          }
        }

        .save_tool_artifact(
          object = result,
          artifact_type = "enrichment_result",
          function_name = "clusterProfiler::enrichKEGG",
          params = list(
            genes = genes,
            organism = organism,
            key_type = key_type,
            readable = readable,
            orgdb_package = orgdb_package,
            pvalueCutoff = pvalueCutoff,
            qvalueCutoff = qvalueCutoff
          ),
          summary = .summarize_enrichment(
            result,
            sprintf("KEGG enrichment on %d gene(s)", length(genes))
          ),
          next_actions = list(
            list(tool = "plot_enrichment_result",
                 reason = "Visualize KEGG pathway enrichment results"),
            list(tool = "interpret_artifact",
                 reason = "AI interpretation of enriched KEGG pathways")
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Retrieve a PPI network with clusterProfiler
#' @keywords internal
#' @noRd
tool_get_ppi_network <- function() {
  aisdk::tool(
    name = "get_ppi_network",
    description = paste(
      "Retrieve a STRING-based PPI network with get_ppi().",
      "Use this when the user provides a gene list and wants a network for downstream topology analysis."
    ),
    parameters = aisdk::z_object(
      genes = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Vector of gene symbols"
      ),
      tax_id = aisdk::z_integer(
        description = "NCBI taxonomy ID, for example 9606 for human"
      )
    ),
    execute = function(genes, tax_id = 9606L) {
      tryCatch({
        genes <- unique(.clean_character_vector(genes))
        result <- get_ppi(genes, taxID = as.integer(tax_id))
        result <- .standardize_ppi_graph(result)

        .save_tool_artifact(
          object = result,
          artifact_type = "ppi_graph",
          function_name = "get_ppi",
          params = list(genes = genes, tax_id = tax_id),
          summary = sprintf(
            "Retrieved a PPI network with %d node(s) and %d edge(s).",
            igraph::vcount(result),
            igraph::ecount(result)
          ),
          next_actions = list(
            list(tool = "compute_ppi_metrics",
                 reason = "Calculate topological metrics (degree, betweenness, closeness) for hub gene identification"),
            list(tool = "run_mcode_clustering",
                 reason = "Detect densely connected modules in the PPI network"),
            list(tool = "subset_ppi_network",
                 reason = "Extract a subnetwork for focused analysis")
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Subset a PPI network
#' @keywords internal
#' @noRd
tool_subset_ppi_network <- function() {
  aisdk::tool(
    name = "subset_ppi_network",
    description = paste(
      "Filter a stored PPI network by edge score and optionally keep the top-degree nodes.",
      "Use this after retrieving or creating a PPI graph artifact."
    ),
    parameters = aisdk::z_object(
      artifact_id = aisdk::z_string(
        description = "Artifact ID of a PPI graph"
      ),
      score_cutoff = aisdk::z_number(
        description = "Minimum edge score to keep"
      ),
      n = aisdk::z_integer(
        description = "Optional number of top-degree nodes to keep; pass 0 to keep all nodes"
      )
    ),
    execute = function(artifact_id, score_cutoff = 0.7, n = 0L) {
      tryCatch({
        ppi_obj <- load_tcm_artifact(artifact_id)
        if (!inherits(ppi_obj, "igraph")) {
          stop("Artifact is not an igraph object.", call. = FALSE)
        }

        n_arg <- .normalize_optional_integer(n, default = NULL)
        result <- ppi_subset(
          ppi_obj = ppi_obj,
          score_cutoff = score_cutoff,
          n = n_arg
        )

        .save_tool_artifact(
          object = result,
          artifact_type = "ppi_graph",
          function_name = "ppi_subset",
          params = list(
            input_artifact = artifact_id,
            score_cutoff = score_cutoff,
            n = n_arg
          ),
          summary = sprintf(
            "PPI subnetwork: %d node(s), %d edge(s), score cutoff %.2f.",
            igraph::vcount(result),
            igraph::ecount(result),
            score_cutoff
          ),
          next_actions = c(
            "compute_ppi_metrics",
            "run_mcode_clustering",
            "plot_ppi_result"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Compute PPI network metrics
#' @keywords internal
#' @noRd
tool_compute_ppi_metrics <- function() {
  aisdk::tool(
    name = "compute_ppi_metrics",
    description = paste(
      "Compute topological metrics for a stored PPI graph.",
      "Use this to annotate the network before ranking hub genes."
    ),
    parameters = aisdk::z_object(
      artifact_id = aisdk::z_string(
        description = "Artifact ID of a PPI graph"
      )
    ),
    execute = function(artifact_id) {
      tryCatch({
        ppi_graph <- load_tcm_artifact(artifact_id)
        if (!inherits(ppi_graph, "igraph")) {
          stop("Artifact is not an igraph object.", call. = FALSE)
        }

        result <- compute_nodeinfo(ppi_graph)

        .save_tool_artifact(
          object = result,
          artifact_type = "ppi_graph",
          function_name = "compute_nodeinfo",
          params = list(input_artifact = artifact_id),
          summary = sprintf(
            "Computed node metrics for %d network node(s).",
            igraph::vcount(result)
          ),
          next_actions = c(
            "rank_ppi_nodes",
            "plot_ppi_result",
            "interpret_artifact"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Rank PPI nodes
#' @keywords internal
#' @noRd
tool_rank_ppi_nodes <- function() {
  aisdk::tool(
    name = "rank_ppi_nodes",
    description = paste(
      "Rank nodes in a PPI network by integrated hub score.",
      "Use this after computing network metrics."
    ),
    parameters = aisdk::z_object(
      artifact_id = aisdk::z_string(
        description = "Artifact ID of a scored PPI graph"
      )
    ),
    execute = function(artifact_id) {
      tryCatch({
        ppi_graph <- load_tcm_artifact(artifact_id)
        if (!inherits(ppi_graph, "igraph")) {
          stop("Artifact is not an igraph object.", call. = FALSE)
        }

        result <- rank_ppi_nodes(ppi_graph)
        rank_df <- .as_rank_table(result)
        top_genes <- paste(utils::head(rank_df$name %||% rank_df[[1]], 5), collapse = ", ")

        .save_tool_artifact(
          object = rank_df,
          artifact_type = "ranking_table",
          function_name = "rank_ppi_nodes",
          params = list(input_artifact = artifact_id),
          summary = sprintf(
            "Ranked %d node(s). Top hits: %s.",
            nrow(rank_df),
            top_genes
          ),
          next_actions = c(
            "plot_ppi_result",
            "interpret_artifact"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Run MCODE clustering
#' @keywords internal
#' @noRd
tool_run_mcode_clustering <- function() {
  aisdk::tool(
    name = "run_mcode_clustering",
    description = paste(
      "Run the MCODE clustering algorithm on a stored PPI graph.",
      "Use this to identify dense functional modules."
    ),
    parameters = aisdk::z_object(
      artifact_id = aisdk::z_string(
        description = "Artifact ID of a PPI graph"
      )
    ),
    execute = function(artifact_id) {
      tryCatch({
        ppi_graph <- load_tcm_artifact(artifact_id)
        if (!inherits(ppi_graph, "igraph")) {
          stop("Artifact is not an igraph object.", call. = FALSE)
        }

        result <- runMCODE(ppi_graph)
        n_clusters <- if (is.list(result)) length(result) else 1L

        .save_tool_artifact(
          object = result,
          artifact_type = "mcode_result",
          function_name = "runMCODE",
          params = list(input_artifact = artifact_id),
          summary = sprintf("MCODE identified %d cluster(s).", n_clusters),
          next_actions = c(
            "interpret_artifact"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Prepare an ML dataset
#' @keywords internal
#' @noRd
tool_prepare_ml_dataset <- function() {
  aisdk::tool(
    name = "prepare_ml_dataset",
    description = paste(
      "Prepare an expression matrix and class labels for machine learning screening.",
      "Use this before running LASSO, Random Forest, SVM-RFE, or XGBoost."
    ),
    parameters = aisdk::z_object(
      expr_matrix = aisdk::z_array(
        items = aisdk::z_array(
          items = aisdk::z_number(),
          description = "A numeric row for one gene across all samples"
        ),
        description = "Nested numeric array with genes in rows and samples in columns"
      ),
      gene_names = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Vector of gene names matching expr_matrix rows"
      ),
      sample_names = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Optional vector of sample names; pass an empty array to auto-generate names"
      ),
      group = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Sample class labels in column order"
      ),
      positive_class = aisdk::z_string(
        description = "Positive class label; pass an empty string to use the default ordering"
      ),
      genes = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Optional candidate gene subset; pass an empty array to keep all genes"
      ),
      split = aisdk::z_boolean(
        description = "Whether to create an internal train/test split"
      ),
      train_ratio = aisdk::z_number(
        description = "Training set fraction when split is TRUE"
      ),
      seed = aisdk::z_integer(
        description = "Random seed"
      )
    ),
    execute = function(expr_matrix,
                       gene_names,
                       sample_names = character(0),
                       group,
                       positive_class = "",
                       genes = character(0),
                       split = FALSE,
                       train_ratio = 0.7,
                       seed = 2025L) {
      tryCatch({
        expr_mat <- .build_expression_matrix(expr_matrix, gene_names, sample_names)
        genes_arg <- .empty_string_to_null(genes)
        positive_class_arg <- .empty_string_to_null(positive_class)

        result <- prepare_ml_data(
          expr_mat = expr_mat,
          group = .clean_character_vector(group),
          positive_class = positive_class_arg,
          genes = genes_arg,
          split = isTRUE(split),
          train_ratio = train_ratio,
          seed = as.integer(seed)
        )

        .save_tool_artifact(
          object = result,
          artifact_type = "ml_data",
          function_name = "prepare_ml_data",
          params = list(
            n_genes = nrow(expr_mat),
            n_samples = ncol(expr_mat),
            positive_class = positive_class_arg,
            genes = genes_arg,
            split = split,
            train_ratio = train_ratio,
            seed = seed
          ),
          summary = .generate_summary(result, "ml_data"),
          next_actions = c(
            "run_ml_screening",
            "interpret_artifact"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Run ML screening
#' @keywords internal
#' @noRd
tool_run_ml_screening <- function() {
  aisdk::tool(
    name = "run_ml_screening",
    description = paste(
      "Run multiple machine learning feature-selection methods on a prepared ML dataset.",
      "Supported methods include lasso, elastic net, ridge, random forest, SVM-RFE, and XGBoost."
    ),
    parameters = aisdk::z_object(
      artifact_id = aisdk::z_string(
        description = "Artifact ID of a prepared ML dataset"
      ),
      methods = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Vector of methods to run; pass an empty array to use the default method set"
      ),
      seed = aisdk::z_integer(
        description = "Random seed"
      ),
      cv_folds = aisdk::z_integer(
        description = "Number of cross-validation folds"
      ),
      cv_repeats = aisdk::z_integer(
        description = "Number of SVM-RFE outer repeats"
      ),
      alpha = aisdk::z_number(
        description = "Elastic Net alpha"
      ),
      top_n = aisdk::z_integer(
        description = "Optional top-n feature cap for methods that support it; pass 0 to keep defaults"
      )
    ),
    execute = function(artifact_id,
                       methods = character(0),
                       seed = 2025L,
                       cv_folds = 5L,
                       cv_repeats = 5L,
                       alpha = 0.5,
                       top_n = 0L) {
      tryCatch({
        ml_data <- load_tcm_artifact(artifact_id)
        if (!inherits(ml_data, "tcm_ml_data")) {
          stop("Artifact is not a tcm_ml_data object.", call. = FALSE)
        }

        methods_arg <- .empty_string_to_null(methods)
        top_n_arg <- .normalize_optional_integer(top_n, default = NULL)
        result <- run_ml_screening(
          ml_data = ml_data,
          methods = methods_arg %||% c("lasso", "rf", "svm_rfe", "xgboost"),
          seed = as.integer(seed),
          cv_folds = as.integer(cv_folds),
          cv_repeats = as.integer(cv_repeats),
          alpha = alpha,
          top_n = top_n_arg
        )

        summary_table <- .ml_result_summary_table(result)
        best_row <- summary_table[which.max(summary_table$auc), , drop = FALSE]
        summary_text <- if (nrow(best_row) == 1 && !is.na(best_row$auc)) {
          sprintf(
            "ML screening completed with %d method(s); best AUC %.3f from %s.",
            nrow(summary_table),
            best_row$auc,
            toupper(best_row$method)
          )
        } else {
          sprintf("ML screening completed with %d method(s).", nrow(summary_table))
        }

        .save_tool_artifact(
          object = result,
          artifact_type = "ml_result",
          function_name = "run_ml_screening",
          params = list(
            input_artifact = artifact_id,
            methods = methods_arg,
            seed = seed,
            cv_folds = cv_folds,
            cv_repeats = cv_repeats,
            alpha = alpha,
            top_n = top_n_arg
          ),
          summary = summary_text,
          next_actions = c(
            "get_ml_consensus",
            "plot_ml_result",
            "interpret_artifact"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Get ML consensus genes
#' @keywords internal
#' @noRd
tool_get_ml_consensus <- function() {
  aisdk::tool(
    name = "get_ml_consensus",
    description = paste(
      "Find consensus genes selected by multiple ML methods.",
      "Use this after ML screening to identify robust candidate markers."
    ),
    parameters = aisdk::z_object(
      artifact_id = aisdk::z_string(
        description = "Artifact ID of an ML screening result"
      ),
      min_methods = aisdk::z_integer(
        description = "Minimum number of methods that must agree on a gene"
      )
    ),
    execute = function(artifact_id, min_methods = 2L) {
      tryCatch({
        ml_list <- load_tcm_artifact(artifact_id)
        if (!inherits(ml_list, "tcm_ml_list")) {
          stop("Artifact is not a tcm_ml_list object.", call. = FALSE)
        }

        consensus <- get_ml_consensus(ml_list, min_methods = as.integer(min_methods))
        feature_counts <- sort(table(unlist(lapply(ml_list, function(model) {
          model$selected_features %||% character(0)
        }))), decreasing = TRUE)
        result <- data.frame(
          gene = names(feature_counts),
          agreement = as.integer(feature_counts),
          stringsAsFactors = FALSE
        )
        result <- result[result$agreement >= as.integer(min_methods), , drop = FALSE]

        .save_tool_artifact(
          object = result,
          artifact_type = "table_result",
          function_name = "get_ml_consensus",
          params = list(input_artifact = artifact_id, min_methods = min_methods),
          summary = sprintf(
            "Found %d consensus gene(s) shared by at least %d method(s).",
            length(consensus),
            as.integer(min_methods)
          ),
          next_actions = c(
            "interpret_artifact"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Retrieve PubMed evidence
#' @keywords internal
#' @noRd
tool_get_pubmed_evidence <- function() {
  aisdk::tool(
    name = "get_pubmed_evidence",
    description = paste(
      "Retrieve PubMed literature for a TCM and disease keyword pair.",
      "Use this to collect supporting evidence and publication trends."
    ),
    parameters = aisdk::z_object(
      tcm_name = aisdk::z_string(
        description = "TCM keyword used in the PubMed query"
      ),
      disease_name = aisdk::z_string(
        description = "Disease keyword used in the PubMed query"
      ),
      start_year = aisdk::z_integer(
        description = "Start year; pass 0 to use the default rolling window"
      ),
      end_year = aisdk::z_integer(
        description = "End year; pass 0 to use the default rolling window"
      ),
      email = aisdk::z_string(
        description = "Email address required by NCBI E-utilities"
      ),
      retmax = aisdk::z_integer(
        description = "Maximum number of records to retrieve"
      )
    ),
    execute = function(tcm_name,
                       disease_name,
                       start_year = 0L,
                       end_year = 0L,
                       email,
                       retmax = 100L) {
      tryCatch({
        year_range <- NULL
        if (start_year > 0 && end_year > 0) {
          year_range <- c(as.integer(start_year), as.integer(end_year))
        }

        result <- get_pubmed_data(
          tcm_name = tcm_name,
          disease_name = disease_name,
          year_range = year_range,
          email = email,
          retmax = as.integer(retmax)
        )

        if (is.null(result)) {
          return(list(
            ok = TRUE,
            summary = "No PubMed records were found for the requested query.",
            articles = 0
          ))
        }

        .save_tool_artifact(
          object = result,
          artifact_type = "pubmed_evidence",
          function_name = "get_pubmed_data",
          params = list(
            tcm_name = tcm_name,
            disease_name = disease_name,
            year_range = year_range,
            retmax = retmax
          ),
          summary = .generate_summary(result, "pubmed_evidence"),
          next_actions = c(
            "extract_pubmed_table",
            "plot_pubmed_result",
            "interpret_artifact"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Extract a PubMed table
#' @keywords internal
#' @noRd
tool_extract_pubmed_table <- function() {
  aisdk::tool(
    name = "extract_pubmed_table",
    description = paste(
      "Extract a sorted article table from a stored PubMed evidence artifact.",
      "Use this when the user wants the top articles as a table."
    ),
    parameters = aisdk::z_object(
      artifact_id = aisdk::z_string(
        description = "Artifact ID of a PubMed evidence object"
      ),
      n = aisdk::z_integer(
        description = "Number of rows to return; pass 0 to return all rows"
      )
    ),
    execute = function(artifact_id, n = 20L) {
      tryCatch({
        pubmed_obj <- load_tcm_artifact(artifact_id)
        if (!inherits(pubmed_obj, "tcm_pubmed")) {
          stop("Artifact is not a tcm_pubmed object.", call. = FALSE)
        }

        n_arg <- .normalize_optional_integer(n, default = NULL)
        result <- get_pubmed_table(pubmed_obj, n = n_arg)

        .save_tool_artifact(
          object = result,
          artifact_type = "table_result",
          function_name = "get_pubmed_table",
          params = list(input_artifact = artifact_id, n = n_arg),
          summary = sprintf("Extracted %d PubMed row(s).", nrow(result)),
          next_actions = c(
            "interpret_artifact"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Search GEO for expression / single-cell datasets
#' @keywords internal
#' @noRd
tool_search_geo_datasets <- function() {
  aisdk::tool(
    name = "search_geo_datasets",
    description = paste(
      "Search NCBI GEO for relevant expression datasets (RNA-seq, microarray, scRNA-seq).",
      "Use this at the START of analysis to help users discover public datasets",
      "for their disease of interest. Returns GSE accessions with titles, summaries,",
      "platforms, and sample counts. The user can then choose which dataset to use",
      "for DEG analysis, WGCNA, or single-cell validation."
    ),
    parameters = aisdk::z_object(
      disease = aisdk::z_string(
        description = "Disease or condition to search for (e.g. 'sepsis', 'diabetic nephropathy')"
      ),
      organism = aisdk::z_string(
        description = "Organism filter; default 'Homo sapiens'"
      ),
      dataset_type = aisdk::z_string(
        description = paste(
          "Dataset type filter: 'rnaseq', 'microarray', 'single-cell', or empty string for all.",
          "Use 'single-cell' to find scRNA-seq datasets for single-cell validation."
        )
      ),
      email = aisdk::z_string(
        description = "Email address required by NCBI E-utilities"
      ),
      max_results = aisdk::z_integer(
        description = "Maximum number of datasets to return; default 20"
      )
    ),
    execute = function(disease, organism = "Homo sapiens", dataset_type = "",
                       email, max_results = 20L) {
      tryCatch({
        dt <- if (nzchar(dataset_type)) dataset_type else NULL

        result <- search_geo_datasets(
          query        = disease,
          organism     = organism,
          dataset_type = dt,
          email        = email,
          retmax       = as.integer(max_results)
        )

        if (is.null(result) || nrow(result) == 0L) {
          return(list(
            ok = TRUE,
            summary = sprintf(
              "No GEO datasets found for '%s' (%s). Try broader search terms.",
              disease, organism
            ),
            datasets = 0L
          ))
        }

        .save_tool_artifact(
          object = result,
          artifact_type = "geo_search_result",
          function_name = "search_geo_datasets",
          params = list(
            disease = disease, organism = organism,
            dataset_type = dataset_type, max_results = max_results
          ),
          summary = sprintf(
            "Found %d GEO dataset(s) for '%s'. Top entries:\n%s",
            nrow(result),
            disease,
            paste(
              utils::head(sprintf("  - %s: %s (n=%s)", result$gse, result$title, result$n_samples), 5),
              collapse = "\n"
            )
          ),
          next_actions = list(
            list(tool = "compute_target_intersection",
                 reason = "After obtaining DEGs from a chosen dataset, compute multi-way intersection"),
            list(tool = "prepare_ml_dataset",
                 reason = "Use expression data for ML-based biomarker screening"),
            list(tool = "interpret_artifact",
                 reason = "Summarize the GEO search results for the user")
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Resolve compound identifiers to PubChem CIDs
#' @keywords internal
#' @noRd
tool_resolve_compound_cid <- function() {
  aisdk::tool(
    name = "resolve_compound_cid",
    description = paste(
      "Resolve compound identifiers to PubChem CIDs.",
      "Use this before requesting compound properties or similarity results."
    ),
    parameters = aisdk::z_object(
      compounds = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Vector of compound identifiers"
      ),
      from = aisdk::z_enum(
        values = c("cid", "smiles", "inchi", "inchikey", "name"),
        description = "Identifier type"
      )
    ),
    execute = function(compounds, from = "name") {
      tryCatch({
        compounds <- .clean_character_vector(compounds)
        cids <- resolve_cid(compounds, from = from)
        result <- data.frame(
          compound = compounds,
          cid = as.character(cids),
          stringsAsFactors = FALSE
        )

        .save_tool_artifact(
          object = result,
          artifact_type = "pubchem_result",
          function_name = "resolve_cid",
          params = list(compounds = compounds, from = from),
          summary = sprintf(
            "Resolved %d CID value(s) from %d compound identifier(s).",
            sum(!is.na(result$cid) & nzchar(result$cid)),
            nrow(result)
          ),
          next_actions = c(
            "get_compound_properties",
            "run_compound_similarity"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Get PubChem compound properties
#' @keywords internal
#' @noRd
tool_get_compound_properties <- function() {
  aisdk::tool(
    name = "get_compound_properties",
    description = paste(
      "Retrieve PubChem compound properties for one or more CIDs.",
      "Use this after resolving CIDs or when the user already has CID values."
    ),
    parameters = aisdk::z_object(
      artifact_id = aisdk::z_string(
        description = "Optional artifact ID from resolve_compound_cid(); pass an empty string to use cids directly"
      ),
      cids = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Vector of PubChem CIDs; pass an empty array to read CIDs from artifact_id"
      ),
      properties = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Vector of PubChem property names; pass an empty array to use the package defaults"
      )
    ),
    execute = function(artifact_id = "", cids = character(0), properties = character(0)) {
      tryCatch({
        cid_values <- unique(c(.clean_character_vector(cids), .extract_cids_from_artifact(artifact_id)))
        if (length(cid_values) == 0) {
          stop("No PubChem CIDs were provided.", call. = FALSE)
        }

        properties_arg <- .empty_string_to_null(properties)
        result <- getprops(cid = cid_values, properties = properties_arg %||% c(
          "MolecularFormula",
          "MolecularWeight",
          "IUPACName",
          "CanonicalSMILES",
          "InChIKey",
          "XLogP"
        ))

        .save_tool_artifact(
          object = result,
          artifact_type = "pubchem_result",
          function_name = "getprops",
          params = list(
            artifact_id = .empty_string_to_null(artifact_id),
            cids = cid_values,
            properties = properties_arg
          ),
          summary = sprintf("Retrieved properties for %d compound(s).", nrow(result)),
          next_actions = c(
            "interpret_artifact"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Run compound similarity search
#' @keywords internal
#' @noRd
tool_run_compound_similarity <- function() {
  aisdk::tool(
    name = "run_compound_similarity",
    description = paste(
      "Run a PubChem 2D similarity search for a query compound.",
      "Use this to find structurally related compounds."
    ),
    parameters = aisdk::z_object(
      query = aisdk::z_string(
        description = "Query identifier"
      ),
      from = aisdk::z_enum(
        values = c("smiles", "cid", "inchikey", "name"),
        description = "Query identifier type"
      ),
      threshold = aisdk::z_integer(
        description = "Similarity threshold from 0 to 100"
      ),
      topn = aisdk::z_integer(
        description = "Maximum number of results"
      ),
      compute_score = aisdk::z_boolean(
        description = "Whether to compute local Tanimoto scores with rcdk"
      )
    ),
    execute = function(query,
                       from = "name",
                       threshold = 90L,
                       topn = 10L,
                       compute_score = TRUE) {
      tryCatch({
        result <- compound_similarity(
          query = query,
          from = from,
          threshold = as.integer(threshold),
          topn = as.integer(topn),
          compute_score = isTRUE(compute_score)
        )

        .save_tool_artifact(
          object = result,
          artifact_type = "pubchem_result",
          function_name = "compound_similarity",
          params = list(
            query = query,
            from = from,
            threshold = threshold,
            topn = topn,
            compute_score = compute_score
          ),
          summary = sprintf("Retrieved %d compound similarity hit(s).", nrow(result)),
          next_actions = c(
            "get_compound_properties",
            "interpret_artifact"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Plot enrichment results
#' @keywords internal
#' @noRd
tool_plot_enrichment_result <- function() {
  aisdk::tool(
    name = "plot_enrichment_result",
    description = paste(
      "Create a plot from a stored enrichment result.",
      "Supported plot types are lollipop and GO bar plot."
    ),
    parameters = aisdk::z_object(
      artifact_id = aisdk::z_string(
        description = "Artifact ID of an enrichment result"
      ),
      plot_type = aisdk::z_enum(
        values = c("lollipop", "go_bar"),
        description = "Plot type"
      ),
      top_n = aisdk::z_integer(
        description = "Number of top terms to plot"
      ),
      plot_title = aisdk::z_string(
        description = "Optional plot title; pass an empty string for no title"
      )
    ),
    execute = function(artifact_id,
                       plot_type = "lollipop",
                       top_n = 10L,
                       plot_title = "") {
      tryCatch({
        enrich_obj <- load_tcm_artifact(artifact_id)
        if (!inherits(enrich_obj, "enrichResult")) {
          stop("Artifact is not an enrichResult object.", call. = FALSE)
        }

        title_arg <- .empty_string_to_null(plot_title)
        plot_obj <- switch(
          plot_type,
          lollipop = gglollipop(
            enrich_obj,
            top_n = as.integer(top_n),
            plot_title = title_arg
          ),
          go_bar = go_barplot(
            enrich_obj,
            top_n = as.integer(top_n),
            plot_title = title_arg
          )
        )

        .save_tool_artifact(
          object = plot_obj,
          artifact_type = "plot",
          function_name = paste0("plot_enrichment_result:", plot_type),
          params = list(
            input_artifact = artifact_id,
            plot_type = plot_type,
            top_n = top_n,
            plot_title = title_arg
          ),
          summary = sprintf("Created an enrichment %s plot.", plot_type),
          next_actions = c(
            "interpret_artifact"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Plot PPI results
#' @keywords internal
#' @noRd
tool_plot_ppi_result <- function() {
  aisdk::tool(
    name = "plot_ppi_result",
    description = paste(
      "Create a plot from ranked PPI results.",
      "Supported plot types are heatmap and radar profile."
    ),
    parameters = aisdk::z_object(
      artifact_id = aisdk::z_string(
        description = "Artifact ID of a PPI graph or ranked PPI table"
      ),
      plot_type = aisdk::z_enum(
        values = c("heatmap", "radar"),
        description = "Plot type"
      ),
      node_name = aisdk::z_string(
        description = "Node name for radar plots; pass an empty string to use the top-ranked node"
      ),
      select_metrics = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Optional metric subset; pass an empty array to use function defaults"
      ),
      plot_title = aisdk::z_string(
        description = "Optional plot title; pass an empty string for the default"
      )
    ),
    execute = function(artifact_id,
                       plot_type = "heatmap",
                       node_name = "",
                       select_metrics = character(0),
                       plot_title = "") {
      tryCatch({
        source_obj <- load_tcm_artifact(artifact_id)
        rank_df <- .as_rank_table(source_obj)
        metrics_arg <- .empty_string_to_null(select_metrics)
        title_arg <- .empty_string_to_null(plot_title)

        plot_obj <- if (plot_type == "heatmap") {
          plot_node_heatmap(
            rank_df,
            select_cols = metrics_arg
          )
        } else {
          if (!"name" %in% names(rank_df)) {
            stop("The ranked PPI table does not contain a 'name' column.", call. = FALSE)
          }
          node_arg <- .empty_string_to_null(node_name) %||% rank_df$name[1]
          profile_df <- get_node_profile(
            rank_df,
            node_name = node_arg,
            metrics = metrics_arg %||% c(
              "degree",
              "betweenness",
              "closeness",
              "MCC",
              "MNC",
              "DMNC",
              "coreness",
              "EPC"
            )
          )
          radar_plot(profile_df, title = title_arg %||% node_arg)
        }

        .save_tool_artifact(
          object = plot_obj,
          artifact_type = "plot",
          function_name = paste0("plot_ppi_result:", plot_type),
          params = list(
            input_artifact = artifact_id,
            plot_type = plot_type,
            node_name = .empty_string_to_null(node_name),
            select_metrics = metrics_arg,
            plot_title = title_arg
          ),
          summary = sprintf("Created a PPI %s plot.", plot_type),
          next_actions = c(
            "interpret_artifact"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Plot ML results
#' @keywords internal
#' @noRd
tool_plot_ml_result <- function() {
  aisdk::tool(
    name = "plot_ml_result",
    description = paste(
      "Create a plot from a stored ML screening result.",
      "Supported plot types are ROC curves, Venn diagrams, and UpSet plots of selected genes."
    ),
    parameters = aisdk::z_object(
      artifact_id = aisdk::z_string(
        description = "Artifact ID of an ML screening result"
      ),
      plot_type = aisdk::z_enum(
        values = c("roc", "venn", "upset"),
        description = "Plot type"
      ),
      set_names = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Optional custom set names for the Venn plot; pass an empty array to infer names automatically"
      )
    ),
    execute = function(artifact_id,
                       plot_type = "roc",
                       set_names = character(0)) {
      tryCatch({
        ml_obj <- load_tcm_artifact(artifact_id)
        if (!inherits(ml_obj, "tcm_ml_list")) {
          stop("Artifact is not a tcm_ml_list object.", call. = FALSE)
        }

        set_names_arg <- .empty_string_to_null(set_names)
        plot_obj <- switch(
          plot_type,
          roc = plot_ml_roc(ml_obj),
          venn = plot_ml_venn(ml_obj, set_names = set_names_arg),
          upset = upsetplot(get_ml_gene_sets(ml_obj, set_names = set_names_arg))
        )

        if (is.null(plot_obj)) {
          stop("No plot data could be generated for the requested ML plot.",
               call. = FALSE)
        }

        .save_tool_artifact(
          object = plot_obj,
          artifact_type = "plot",
          function_name = paste0("plot_ml_result:", plot_type),
          params = list(
            input_artifact = artifact_id,
            plot_type = plot_type,
            set_names = set_names_arg
          ),
          summary = sprintf("Created an ML %s plot.", plot_type),
          next_actions = c(
            "interpret_artifact"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Plot PubMed results
#' @keywords internal
#' @noRd
tool_plot_pubmed_result <- function() {
  aisdk::tool(
    name = "plot_pubmed_result",
    description = paste(
      "Create a publication-trend or journal-frequency plot from a stored PubMed result.",
      "Use this to visualize the PubMed evidence artifact."
    ),
    parameters = aisdk::z_object(
      artifact_id = aisdk::z_string(
        description = "Artifact ID of a PubMed evidence object"
      ),
      plot_type = aisdk::z_enum(
        values = c("trend", "journal"),
        description = "Plot type"
      ),
      top_n = aisdk::z_integer(
        description = "Number of journals to show for the journal plot"
      )
    ),
    execute = function(artifact_id,
                       plot_type = "trend",
                       top_n = 10L) {
      tryCatch({
        pubmed_obj <- load_tcm_artifact(artifact_id)
        if (!inherits(pubmed_obj, "tcm_pubmed")) {
          stop("Artifact is not a tcm_pubmed object.", call. = FALSE)
        }

        plot_obj <- plot(pubmed_obj, type = plot_type, N = as.integer(top_n))
        if (is.null(plot_obj)) {
          stop("No plot could be generated from the PubMed artifact.", call. = FALSE)
        }

        .save_tool_artifact(
          object = plot_obj,
          artifact_type = "plot",
          function_name = paste0("plot_pubmed_result:", plot_type),
          params = list(input_artifact = artifact_id, plot_type = plot_type, top_n = top_n),
          summary = sprintf("Created a PubMed %s plot.", plot_type),
          next_actions = c(
            "interpret_artifact"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Plot a herb Sankey diagram
#' @keywords internal
#' @noRd
tool_plot_herb_sankey <- function() {
  aisdk::tool(
    name = "plot_herb_sankey",
    description = paste(
      "Create a herb-molecule-target Sankey diagram from a search result table.",
      "Use this after search_herb_records() or search_targets()."
    ),
    parameters = aisdk::z_object(
      artifact_id = aisdk::z_string(
        description = "Artifact ID of a search result table"
      ),
      axis_order = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Axis order for the Sankey diagram; pass an empty array to use herb-molecule-target"
      )
    ),
    execute = function(artifact_id, axis_order = character(0)) {
      tryCatch({
        search_df <- load_tcm_artifact(artifact_id)
        if (!is.data.frame(search_df)) {
          stop("Artifact is not a data.frame.", call. = FALSE)
        }

        required_cols <- c("herb", "molecule", "target")
        missing_cols <- setdiff(required_cols, names(search_df))
        if (length(missing_cols) > 0) {
          stop(sprintf(
            "The artifact is missing required columns: %s.",
            paste(missing_cols, collapse = ", ")
          ), call. = FALSE)
        }

        axis_order_arg <- .empty_string_to_null(axis_order) %||% required_cols
        plot_obj <- tcm_sankey(search_df, axis_order = axis_order_arg)

        .save_tool_artifact(
          object = plot_obj,
          artifact_type = "plot",
          function_name = "tcm_sankey",
          params = list(input_artifact = artifact_id, axis_order = axis_order_arg),
          summary = "Created a herb-molecule-target Sankey diagram.",
          next_actions = c(
            "interpret_artifact"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Plot a docking heatmap or dot heatmap
#' @keywords internal
#' @noRd
tool_plot_docking_heatmap <- function() {
  aisdk::tool(
    name = "plot_docking_heatmap",
    description = paste(
      "Create a docking affinity plot from a numeric matrix or data frame artifact.",
      "Supported plot types are dot and tile."
    ),
    parameters = aisdk::z_object(
      artifact_id = aisdk::z_string(
        description = "Artifact ID of a docking matrix or numeric data frame"
      ),
      plot_type = aisdk::z_enum(
        values = c("dot", "tile"),
        description = "Plot type"
      ),
      order = aisdk::z_enum(
        values = c("none", "median", "mean", "max", "min"),
        description = "Reordering method"
      ),
      label = aisdk::z_boolean(
        description = "Whether to display numeric affinity labels"
      )
    ),
    execute = function(artifact_id,
                       plot_type = "dot",
                       order = "none",
                       label = FALSE) {
      tryCatch({
        docking_obj <- .coerce_docking_matrix(load_tcm_artifact(artifact_id))
        plot_obj <- ggdock(
          dock_data = docking_obj,
          type = plot_type,
          order = if (identical(order, "none")) NULL else order,
          label = isTRUE(label)
        )

        .save_tool_artifact(
          object = plot_obj,
          artifact_type = "plot",
          function_name = paste0("ggdock:", plot_type),
          params = list(
            input_artifact = artifact_id,
            plot_type = plot_type,
            order = order,
            label = label
          ),
          summary = sprintf("Created a docking %s plot.", plot_type),
          next_actions = c(
            "interpret_artifact"
          )
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: List artifacts
#' @keywords internal
#' @noRd
tool_list_artifacts <- function() {
  aisdk::tool(
    name = "list_artifacts",
    description = "List all stored analysis artifacts created in the current R session.",
    parameters = aisdk::z_object(
      dummy = aisdk::z_string(
        description = "Unused placeholder; pass an empty string"
      )
    ),
    execute = function(dummy = "") {
      tryCatch({
        df <- list_tcm_artifacts()
        list(
          ok = TRUE,
          summary = sprintf("%d artifact(s) currently stored.", nrow(df)),
          artifacts = df
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Load artifact summary
#' @keywords internal
#' @noRd
tool_load_artifact_summary <- function() {
  aisdk::tool(
    name = "load_artifact_summary",
    description = "Return a text summary and preview for a stored artifact.",
    parameters = aisdk::z_object(
      artifact_id = aisdk::z_string(
        description = "Artifact ID to summarize"
      )
    ),
    execute = function(artifact_id) {
      tryCatch({
        summary_text <- summarize_tcm_artifact(artifact_id)
        list(
          ok = TRUE,
          artifact_id = artifact_id,
          summary = summary_text
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Interpret an artifact with AI
#' @keywords internal
#' @noRd
tool_interpret_artifact <- function() {
  aisdk::tool(
    name = "interpret_artifact",
    description = paste(
      "Interpret a stored artifact with the AI layer.",
      "The tool automatically dispatches to enrichment, PPI, or table interpretation when possible."
    ),
    parameters = aisdk::z_object(
      artifact_id = aisdk::z_string(
        description = "Artifact ID to interpret"
      ),
      language = aisdk::z_enum(
        values = c("en", "zh"),
        description = "Output language"
      )
    ),
    execute = function(artifact_id, language = "en") {
      tryCatch({
        obj <- load_tcm_artifact(artifact_id)

        result <- if (inherits(obj, "enrichResult")) {
          interpret_enrichment(obj, language = language)
        } else if (inherits(obj, "igraph") || (is.list(obj) && !is.null(obj$table))) {
          interpret_ppi(obj, language = language)
        } else if (inherits(obj, "tcm_pubmed")) {
          interpret_table(get_pubmed_table(obj, n = 10), language = language)
        } else if (inherits(obj, "tcm_ml_list")) {
          interpret_table(.ml_result_summary_table(obj), language = language)
        } else if (inherits(obj, "tcm_ml")) {
          interpret_table(utils::head(obj$importance %||% data.frame(), 10), language = language)
        } else if (is.data.frame(obj)) {
          interpret_table(obj, language = language)
        } else {
          tcm_interpret(summarize_tcm_artifact(artifact_id), language = language)
        }

        interpretation <- if (is.character(result)) {
          paste(result, collapse = "\n")
        } else {
          result$output$summary %||%
            result$output$paragraph %||%
            as.character(result)
        }

        list(
          ok = TRUE,
          artifact_id = artifact_id,
          interpretation = interpretation
        )
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Read an R variable from the user's global environment
#' @keywords internal
#' @noRd
tool_read_variable <- function() {
  aisdk::tool(
    name = "read_r_variable",
    description = paste(
      "Read and summarize an R variable from the user's current R session.",
      "Returns type, dimensions, and a preview of the variable's contents.",
      "Use this when the user refers to a variable name in their request",
      "(e.g. 'use my_genes for enrichment', 'analyze the deg_list').",
      "For data.frames: shows columns, types, and first 5 rows.",
      "For vectors: shows length, type, and first 10 elements.",
      "For lists: shows element names and types."
    ),
    parameters = aisdk::z_object(
      name = aisdk::z_string(
        description = "The variable name to read from the R environment"
      )
    ),
    execute = function(name) {
      tryCatch({
        name <- trimws(as.character(name))
        if (!nzchar(name)) {
          return(list(ok = FALSE, error = "Variable name cannot be empty."))
        }

        envir <- globalenv()
        if (!exists(name, envir = envir, inherits = FALSE)) {
          # Also check package namespace
          ns <- asNamespace("TCMDATA")
          if (exists(name, envir = ns, inherits = FALSE)) {
            envir <- ns
          } else {
            available <- ls(globalenv())
            hint <- if (length(available) > 0) {
              paste("Available variables:", paste(head(available, 20), collapse = ", "))
            } else {
              "No variables found in global environment."
            }
            return(list(ok = FALSE, error = paste0(
              "Variable '", name, "' not found. ", hint
            )))
          }
        }

        summary_text <- aisdk::get_r_context(name, envir = envir)
        list(ok = TRUE, variable = name, summary = summary_text)
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}

#' Tool: Evaluate R code in a sandboxed environment
#' @keywords internal
#' @noRd
tool_eval_r_code <- function() {
  # Create a sandbox once per tool instance
  sandbox <- NULL

  aisdk::tool(
    name = "eval_r_code",
    description = paste(
      "Execute a short R code snippet and return the output.",
      "The code runs in a sandboxed environment with access to base R, dplyr, and purrr.",
      "Variables from the user's globalenv() are accessible read-only.",
      "Use this for data transformations the user requests that no existing tool handles,",
      "such as subsetting, filtering, basic statistics, or preparing data for other tools.",
      "Keep code short (< 20 lines). Do NOT use this to install packages or write files."
    ),
    parameters = aisdk::z_object(
      code = aisdk::z_string(
        description = "R code to evaluate (single string, can contain newlines)"
      )
    ),
    execute = function(code) {
      tryCatch({
        code <- as.character(code)
        if (!nzchar(trimws(code))) {
          return(list(ok = FALSE, error = "Code cannot be empty."))
        }

        # Lazy-init sandbox with globalenv as parent for read access
        if (is.null(sandbox)) {
          sandbox <<- aisdk::SandboxManager$new(
            parent_env = globalenv(),
            max_output_chars = 8000
          )
        }

        output <- sandbox$execute(code)
        list(ok = TRUE, output = output)
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}


#' Tool: Generate cross-validation URLs for external TCM databases
#' @keywords internal
#' @noRd
tool_generate_verification_urls <- function() {
  aisdk::tool(
    name = "generate_verification_urls",
    description = paste(
      "Generate verification URLs for HERB database and ShenNong Alpha to cross-check",
      "herb-target or herb-disease findings. Use this after key analysis steps",
      "(target retrieval, PPI hub genes, enrichment) to provide the user with",
      "external cross-validation links. This helps reduce hallucination by grounding",
      "findings in independent databases."
    ),
    parameters = aisdk::z_object(
      herbs = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Herb names (Chinese or English) to generate verification links for"
      ),
      targets = aisdk::z_array(
        items = aisdk::z_string(),
        description = "Target gene symbols to include in verification (optional, can be empty)"
      ),
      disease = aisdk::z_string(
        description = "Disease name for context (optional, can be empty string)"
      )
    ),
    execute = function(herbs = character(0), targets = character(0), disease = "") {
      tryCatch({
        herbs <- .clean_character_vector(herbs)
        targets <- .clean_character_vector(targets)
        disease <- as.character(disease)

        if (length(herbs) == 0L) {
          return(list(ok = FALSE, error = "At least one herb name is required."))
        }

        urls <- lapply(herbs, function(herb) {
          herb_encoded <- utils::URLencode(herb, reserved = TRUE)
          list(
            herb = herb,
            herb_url = sprintf("http://47.92.70.12/#/Browse?keyword=%s", herb_encoded),
            herb_2_url = sprintf("http://47.92.70.12/#/Browse?keyword=%s", herb_encoded),
            shennong_url = sprintf(
              "https://shennongalpha.westlake.edu.cn/#/herb?keyword=%s",
              herb_encoded
            )
          )
        })

        target_urls <- if (length(targets) > 0L) {
          lapply(targets[seq_len(min(10L, length(targets)))], function(gene) {
            gene_encoded <- utils::URLencode(gene, reserved = TRUE)
            list(
              gene = gene,
              herb_target_url = sprintf(
                "http://47.92.70.12/#/Browse?keyword=%s",
                gene_encoded
              )
            )
          })
        }

        evidence_note <- paste(
          "Cross-validation protocol:",
          "1) Check HERB database for experiment-based (gene expression) evidence",
          "2) Check HERB 2.0 for clinical trial / meta-analysis support",
          "3) Check ShenNong for classical TCM text attributions",
          "Evidence hierarchy: Clinical trial > Experimental > Computational > Text-mining"
        )

        result <- list(
          ok = TRUE,
          herb_links = urls,
          target_links = target_urls,
          disease = if (nzchar(disease)) disease else NULL,
          evidence_note = evidence_note,
          summary = sprintf(
            "Generated verification URLs for %d herb(s) and %d target(s) across HERB and ShenNong databases.",
            length(herbs),
            length(targets)
          )
        )

        result
      }, error = function(e) {
        list(ok = FALSE, error = conditionMessage(e))
      })
    }
  )
}