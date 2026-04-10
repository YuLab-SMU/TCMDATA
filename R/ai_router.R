# ai_router.R
# Task router for the TCM agent layer.
# Routes user tasks to appropriate tool sets using rules with an optional LLM fallback.

.tcm_task_types <- c(
  "herb_lookup",
  "target_lookup",
  "disease_lookup",
  "enrichment",
  "ppi_analysis",
  "ml_screening",
  "pubmed",
  "compound",
  "visualization",
  "interpretation",
  "general"
)

#' Route a task to the appropriate tool set
#'
#' Analyzes the user's task description and determines which tools should be
#' made available to the agent. Uses rule-based routing first, then optionally
#' falls back to LLM routing for ambiguous cases.
#'
#' @param task Character. The user's task description.
#' @param use_llm Logical. Whether to use LLM for ambiguous routing.
#' @param model Model object for LLM routing when \code{use_llm = TRUE}.
#'
#' @return A list with components:
#' \describe{
#'   \item{task_type}{Character. The identified task type.}
#'   \item{tools}{Character vector of tool names to enable.}
#'   \item{confidence}{Character. "high" for rule-based or "medium" for LLM-assisted.}
#'   \item{source_hint}{Character. Preferred source: "local", "api", or "artifact".}
#' }
#'
#' @examples
#' \dontrun{
#'   route_tcm_task("Run GO enrichment on the target genes")
#'   route_tcm_task("Retrieve a PPI network and rank hub genes")
#' }
#' @export
route_tcm_task <- function(task, use_llm = FALSE, model = NULL) {
  task_lower <- tolower(task)

  # --- Analysis task types (checked first) ---
  rules <- list(
    pubmed = list(
      patterns = c(
        "pubmed", "literature", "pmid", "journal", "paper",
        "evidence retrieval", "publication trend",
        "\u6587\u732e", "\u8bba\u6587", "\u53d1\u8868"
      ),
      tools = c(
        "get_pubmed_evidence", "extract_pubmed_table",
        "plot_pubmed_result", "interpret_artifact"
      ),
      source_hint = "api"
    ),
    compound = list(
      patterns = c(
        "pubchem", "cid", "compound", "smiles", "inchi",
        "inchikey", "similarity", "molecular weight", "iupac"
      ),
      tools = c(
        "resolve_compound_cid", "get_compound_properties",
        "run_compound_similarity", "plot_docking_heatmap",
        "interpret_artifact"
      ),
      source_hint = "api"
    ),
    ml_screening = list(
      patterns = c(
        "machine.*learning", "feature.*selection", "biomarker",
        "lasso", "elastic.*net", "ridge", "random.*forest",
        "svm", "xgboost", "ml screening"
      ),
      tools = c(
        "prepare_ml_dataset", "run_ml_screening", "get_ml_consensus",
        "plot_ml_result", "interpret_artifact"
      ),
      source_hint = "local"
    ),
    ppi_analysis = list(
      patterns = c(
        "getppi", "get ppi", "ppi", "protein.*interaction",
        "string network", "hub gene", "centrality", "mcode",
        "module", "subnetwork", "cluster",
        "\u86cb\u767d.*\u4e92\u4f5c", "\u6838\u5fc3\u57fa\u56e0"
      ),
      tools = c(
        "get_ppi_network", "subset_ppi_network", "compute_ppi_metrics",
        "rank_ppi_nodes", "run_mcode_clustering", "plot_ppi_result",
        "interpret_artifact"
      ),
      source_hint = "api"
    ),
    enrichment = list(
      patterns = c(
        "enrichment", "enrich", "clusterprofiler", "go\\b",
        "kegg", "pathway", "ora",
        "\u5bcc\u96c6", "\u901a\u8def", "\u4fe1\u53f7\u901a\u8def"
      ),
      tools = c(
        "run_herb_enrichment", "run_go_enrichment", "run_kegg_enrichment",
        "plot_enrichment_result", "interpret_artifact"
      ),
      source_hint = "local"
    ),
    target_lookup = list(
      patterns = c(
        "search.*target", "target.*search", "gene.*target",
        "protein.*target", "target.*herb", "gene.*info",
        "\u57fa\u56e0.*\u67e5", "\u67e5.*\u57fa\u56e0"
      ),
      tools = c("search_targets", "list_artifacts", "load_artifact_summary",
                "interpret_artifact"),
      source_hint = "local"
    ),
    herb_lookup = list(
      patterns = c(
        "search.*herb", "herb.*search", "herb.*component",
        "herb.*target", "herb record", "compound profile",
        "\u4e2d\u836f", "\u8349\u836f", "\u7ec4\u65b9", "\u65b9\u5242",
        "\u67e5.*\u836f", "\u836f.*\u67e5",
        "\u9ec4\u82aa", "\u9ec4\u82a9", "\u5f53\u5f52", "\u4eba\u53c2",
        "\u7518\u8349", "\u767d\u672f", "\u832f\u82d3", "\u67f4\u80e1",
        "\u5ddd\u828e", "\u9ec4\u8fde", "\u91d1\u94f6\u82b1",
        "\u8fde\u7fd8", "\u9752\u84bf", "\u5730\u9ec4", "\u7ea2\u82b1",
        "\u6842\u679d", "\u9644\u5b50", "\u5927\u9ec4", "\u9ebb\u9ec4",
        "\u7ec4\u5206|\u6210\u5206", "\u9776\u70b9|\u9776\u6807"
      ),
      tools = c("search_herb_records", "plot_herb_sankey"),
      source_hint = "local"
    ),
    disease_lookup = list(
      patterns = c(
        "disease", "disea", "disgenet", "gene.*disease",
        "disease.*gene", "disease.*target", "target.*disease",
        "\\bsepsis\\b", "\\bdiabetes\\b", "\\basthma\\b", "\\bcancer\\b",
        "\\btumor\\b", "\\btumour\\b",
        "umls", "cui", "\\bC\\d{7}\\b",
        "\u75be\u75c5", "\u75c5\u75c7", "\u7cd6\u5c3f\u75c5",
        "\u764c", "\u80bf\u7624", "\u51a0\u5fc3\u75c5",
        "\u809d\u708e", "\u54ee\u5598"
      ),
      tools = c(
        "search_disease_targets", "search_gene_diseases",
        "run_go_enrichment", "run_kegg_enrichment",
        "get_ppi_network", "interpret_artifact"
      ),
      source_hint = "local"
    ),
    # visualization and interpretation checked LAST
    geo_search = list(
      patterns = c(
        "geo", "gse\\d", "dataset", "\\brna.?seq\\b",
        "transcriptom", "microarray", "expression.*profil",
        "single.?cell", "scrna", "wgcna", "co.?expression"
      ),
      tools = c(
        "search_geo_datasets", "interpret_artifact"
      ),
      source_hint = "api"
    ),
    visualization = list(
      patterns = c(
        "visuali[sz]e", "plot", "figure", "chart", "heatmap",
        "sankey", "lollipop", "roc", "venn", "upset", "radar",
        "\u53ef\u89c6\u5316", "\u753b.*\u56fe", "\u7ed8\u56fe", "\u56fe\u8868",
        "\u70ed\u56fe", "\u6851\u57fa\u56fe", "\u97e6\u6069\u56fe"
      ),
      tools = c(
        "plot_enrichment_result", "plot_ppi_result", "plot_ml_result",
        "plot_pubmed_result", "plot_herb_sankey", "plot_docking_heatmap"
      ),
      source_hint = "artifact"
    ),
    interpretation = list(
      patterns = c(
        "interpret", "explain", "summari[sz]e", "write.*result",
        "biological.*meaning", "draft",
        "\u89e3\u8bfb", "\u89e3\u91ca", "\u603b\u7ed3", "\u5206\u6790.*\u7ed3\u679c"
      ),
      tools = c("list_artifacts", "load_artifact_summary", "interpret_artifact"),
      source_hint = "artifact"
    )
  )

  # Collect ALL matching task types, not just the first
  matched <- list()
  for (task_type in names(rules)) {
    rule <- rules[[task_type]]
    if (any(vapply(rule$patterns, function(pattern) {
      grepl(pattern, task_lower, ignore.case = TRUE)
    }, logical(1)))) {
      matched[[task_type]] <- rule
    }
  }

  if (length(matched) > 0) {
    # Primary = first matched (by priority order above)
    primary_type <- names(matched)[1]
    # Merge tools from all matched types
    all_tools <- unique(unlist(lapply(matched, function(r) r$tools)))
    return(list(
      task_type   = primary_type,
      tools       = all_tools,
      confidence  = "high",
      source_hint = matched[[primary_type]]$source_hint
    ))
  }

  if (use_llm && !is.null(model)) {
    return(.route_with_llm(task, model))
  }

  list(
    task_type = "general",
    tools = c(
      "search_herb_records", "search_targets", "run_go_enrichment",
      "get_ppi_network", "prepare_ml_dataset", "get_pubmed_evidence",
      "resolve_compound_cid", "plot_enrichment_result", "interpret_artifact"
    ),
    confidence = "low",
    source_hint = "unknown"
  )
}

#' LLM-based routing for ambiguous tasks
#' @keywords internal
#' @noRd
.route_with_llm <- function(task, model) {
  .check_aisdk()

  schema <- aisdk::z_object(
    task_type = aisdk::z_enum(
      values = .tcm_task_types,
      description = "The detected task type"
    ),
    source_hint = aisdk::z_enum(
      values = c("local", "api", "artifact"),
      description = "Preferred source type"
    )
  )

  system_prompt <- paste(
    "<role>You are a task classifier for a TCM network pharmacology analysis system.</role>",
    "",
    "<task>Classify the user request into exactly one task_type and one source_hint.</task>",
    "",
    "<task_types>",
    "herb_lookup     - herb, molecule, or target lookup from the local database",
    "target_lookup   - target- or gene-centric search",
    "disease_lookup  - disease-gene association queries using DisGeNET",
    "enrichment      - GO, KEGG, or herb functional enrichment analysis",
    "ppi_analysis    - PPI retrieval, topology metrics, hub gene ranking, or MCODE clustering",
    "ml_screening    - machine learning feature selection or consensus biomarkers",
    "pubmed          - PubMed evidence retrieval or publication trend analysis",
    "compound        - PubChem CID, molecular properties, or structural similarity",
    "visualization   - plot, chart, or figure generation from stored artifacts",
    "interpretation  - explanation, summary, or result paragraph writing",
    "general         - anything that does not fit the above categories",
    "</task_types>",
    "",
    "<source_hints>",
    "local    - data available in the local database or package",
    "api      - requires an external API call (STRING, PubMed, PubChem)",
    "artifact - operates on previously stored analysis artifacts",
    "</source_hints>"
  )

  result <- tryCatch({
    aisdk::generate_object(
      model = model,
      prompt = task,
      schema = schema,
      system = system_prompt
    )
  }, error = function(e) {
    list(object = list(task_type = "general", source_hint = "unknown"))
  })

  obj <- result$object %||% list(task_type = "general", source_hint = "unknown")

  tool_map <- list(
    herb_lookup = c("search_herb_records", "plot_herb_sankey"),
    target_lookup = c("search_targets", "list_artifacts", "load_artifact_summary",
                       "interpret_artifact"),
    enrichment = c(
      "run_herb_enrichment", "run_go_enrichment", "run_kegg_enrichment",
      "plot_enrichment_result", "interpret_artifact"
    ),
    ppi_analysis = c(
      "get_ppi_network", "subset_ppi_network", "compute_ppi_metrics",
      "rank_ppi_nodes", "run_mcode_clustering", "plot_ppi_result",
      "interpret_artifact"
    ),
    ml_screening = c(
      "prepare_ml_dataset", "run_ml_screening", "get_ml_consensus",
      "plot_ml_result", "interpret_artifact"
    ),
    pubmed = c(
      "get_pubmed_evidence", "extract_pubmed_table",
      "plot_pubmed_result", "interpret_artifact"
    ),
    compound = c(
      "resolve_compound_cid", "get_compound_properties",
      "run_compound_similarity", "plot_docking_heatmap",
      "interpret_artifact"
    ),
    visualization = c(
      "plot_enrichment_result", "plot_ppi_result", "plot_ml_result",
      "plot_pubmed_result", "plot_herb_sankey", "plot_docking_heatmap"
    ),
    interpretation = c("list_artifacts", "load_artifact_summary", "interpret_artifact"),
    general = c(
      "search_herb_records", "search_targets", "run_go_enrichment",
      "get_ppi_network", "prepare_ml_dataset", "get_pubmed_evidence",
      "resolve_compound_cid", "plot_enrichment_result", "interpret_artifact"
    )
  )

  list(
    task_type = obj$task_type,
    tools = tool_map[[obj$task_type]] %||% tool_map$general,
    confidence = "medium",
    source_hint = obj$source_hint %||% "unknown"
  )
}

#' Resolve artifact references in task text
#'
#' Parses task text for explicit or natural-language references to previous
#' artifacts and resolves them to artifact IDs when possible.
#'
#' @param task Character. The task text.
#'
#' @return A list with:
#'   \\item{resolved_task}{Task text with references resolved.}
#'   \\item{artifact_ids}{Character vector of resolved artifact IDs.}
#'
#' @keywords internal
#' @noRd
resolve_artifact_references <- function(task) {
  artifact_ids <- character(0)

  explicit_pattern <- "(enrich|hgraph|ppi|rank|mcode|ml|pubmed|chem|table|plot|search|art)_\\d{3}"
  explicit_matches <- regmatches(task, gregexpr(explicit_pattern, task))[[1]]
  if (length(explicit_matches) > 0) {
    for (artifact_id in explicit_matches) {
      if (artifact_exists(artifact_id)) {
        artifact_ids <- c(artifact_ids, artifact_id)
      }
    }
  }

  nl_patterns <- c(
    "the (previous|last) (enrichment )?(result|analysis)" = "enrichment_result",
    "the (previous|last) ppi (network|graph|result)" = "ppi_graph",
    "the (previous|last) ml (result|screening)" = "ml_result",
    "the (previous|last) pubmed (result|search)" = "pubmed_evidence",
    "the (previous|last) plot" = "plot",
    "the (previous|last) table" = "table_result"
  )

  for (pattern in names(nl_patterns)) {
    if (grepl(pattern, task, ignore.case = TRUE)) {
      target_type <- nl_patterns[[pattern]]
      artifacts <- list_tcm_artifacts()
      if (nrow(artifacts) > 0) {
        matches <- artifacts[artifacts$artifact_type == target_type, , drop = FALSE]
        if (nrow(matches) > 0) {
          artifact_ids <- c(artifact_ids, matches$artifact_id[nrow(matches)])
        }
      }
    }
  }

  list(
    resolved_task = task,
    artifact_ids = unique(artifact_ids)
  )
}

