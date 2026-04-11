# tests/testthat/test-ai.R
# Unit tests for the aisdk integration layer.
# All tests mock aisdk::generate_object() — no live API calls.

# ── Helper: mock generate_object ─────────────────────────────────────────────

mock_generate_object <- function(...) {
  list(
    object = list(
      summary = "Mock summary of analysis results.",
      key_findings = list("Finding 1", "Finding 2"),
      biological_interpretation = "Mock biological interpretation.",
      tcm_relevance = "Mock TCM relevance note.",
      caveats = list("Caveat 1")
    )
  )
}

mock_generate_object_draft <- function(...) {
  list(
    object = list(
      paragraph = "Mock result paragraph for publication.",
      figure_legend_hint = "Mock figure legend."
    )
  )
}

# ── Tests: dependency check ──────────────────────────────────────────────────

test_that(".check_aisdk errors when aisdk is missing", {
  # Only run if aisdk is NOT installed
  skip_if(requireNamespace("aisdk", quietly = TRUE),
          "aisdk is installed — skip missing-package test")
  expect_error(.check_aisdk(), "aisdk")
})

# ── Tests: context compression ───────────────────────────────────────────────

test_that(".compress_enrichment handles data.frame", {
  df <- data.frame(
    ID = c("GO:0001", "GO:0002"),
    p.adjust = c(0.001, 0.05),
    GeneRatio = c("3/100", "5/100"),
    geneID = c("TP53/BRCA1/EGFR", "AKT1/MTOR"),
    stringsAsFactors = FALSE
  )
  ctx <- .compress_enrichment(df, top_n = 2)
  expect_type(ctx, "character")
  expect_true(grepl("GO:0001", ctx))
  expect_true(grepl("GO:0002", ctx))
})

test_that(".compress_ppi handles data.frame", {
  df <- data.frame(
    name = c("TP53", "AKT1", "EGFR"),
    degree = c(15, 12, 10),
    betweenness = c(0.5, 0.3, 0.2),
    stringsAsFactors = FALSE
  )
  ctx <- .compress_ppi(df, top_n = 2)
  expect_type(ctx, "character")
  expect_true(grepl("TP53", ctx))
})

test_that(".compress_table handles generic data.frame", {
  df <- data.frame(
    gene = c("TP53", "BRCA1", "EGFR"),
    logFC = c(2.5, -1.8, 3.1),
    p.adjust = c(0.001, 0.01, 0.005),
    stringsAsFactors = FALSE
  )
  ctx <- .compress_table(df, top_n = 3)
  expect_type(ctx, "character")
  expect_true(grepl("TP53", ctx))
  expect_true(grepl("logFC", ctx))
})

# ── Tests: interpret functions (mocked) ──────────────────────────────────────

test_that("interpret_enrichment returns tcm_ai_analysis", {
  skip_if_not_installed("aisdk")
  local_mocked_bindings(
    generate_object = mock_generate_object,
    .package = "aisdk"
  )

  df <- data.frame(
    ID = c("GO:0001", "GO:0002"),
    p.adjust = c(0.001, 0.05),
    GeneRatio = c("3/100", "5/100"),
    geneID = c("TP53/BRCA1", "AKT1/MTOR"),
    stringsAsFactors = FALSE
  )

  res <- interpret_enrichment(df, top_n = 2)
  expect_s3_class(res, "tcm_ai_analysis")
  expect_equal(res$input$type, "enrichment")
  expect_true(!is.null(res$output$summary))
  expect_true(!is.null(res$metadata$model))
})

test_that("interpret_table returns tcm_ai_analysis", {
  skip_if_not_installed("aisdk")
  local_mocked_bindings(
    generate_object = mock_generate_object,
    .package = "aisdk"
  )

  df <- data.frame(
    gene = c("TP53", "BRCA1"),
    logFC = c(2.5, -1.8),
    p.adjust = c(0.001, 0.01),
    stringsAsFactors = FALSE
  )

  res <- interpret_table(df, top_n = 2)
  expect_s3_class(res, "tcm_ai_analysis")
  expect_equal(res$input$type, "table")
})

# ── Tests: draft function (mocked) ──────────────────────────────────────────

test_that("draft_result_paragraph returns tcm_ai_draft", {
  skip_if_not_installed("aisdk")
  local_mocked_bindings(
    generate_object = mock_generate_object_draft,
    .package = "aisdk"
  )

  df <- data.frame(
    ID = c("GO:0001"),
    p.adjust = c(0.001),
    GeneRatio = c("3/100"),
    geneID = c("TP53/BRCA1"),
    stringsAsFactors = FALSE
  )

  res <- draft_result_paragraph(df, type = "enrichment")
  expect_s3_class(res, "tcm_ai_draft")
  expect_true(!is.null(res$draft$paragraph))
})

# ── Tests: print methods ────────────────────────────────────────────────────

test_that("print.tcm_ai_analysis works", {
  obj <- .new_tcm_ai_analysis(
    input = list(type = "enrichment"),
    context = "test context",
    output = list(
      summary = "Test summary",
      key_findings = list("F1"),
      biological_interpretation = "Test interp",
      tcm_relevance = "TCM note",
      caveats = list("C1")
    ),
    metadata = list(
      model = "test-model", language = "en", audience = "researcher",
      input_class = "data.frame", generated_at = Sys.time(),
      prompt_version = "1.0"
    )
  )
  expect_output(print(obj), "TCM AI Analysis")
})

test_that("print.tcm_ai_draft works", {
  obj <- .new_tcm_ai_draft(
    input = list(type = "enrichment"),
    context = "test context",
    draft = list(
      paragraph = "Test paragraph.",
      figure_legend_hint = "Test legend."
    ),
    metadata = list(
      model = "test-model", language = "en", audience = "paper",
      input_class = "data.frame", generated_at = Sys.time(),
      prompt_version = "1.0"
    )
  )
  expect_output(print(obj), "TCM AI Draft")
})

# ── Regression: dispatch, text-mode params, provider whitelist ───────────────

test_that("tcm_interpret type= overrides class-based auto-detection", {
  skip_if_not_installed("aisdk")
  local_mocked_bindings(
    generate_object = mock_generate_object,
    .package = "aisdk"
  )
  # A plain data.frame would normally dispatch to interpret_table();
  # type = "enrichment" must force it to interpret_enrichment() instead.
  df <- data.frame(
    ID = c("GO:0001", "GO:0002"),
    p.adjust = c(0.001, 0.05),
    GeneRatio = c("3/100", "5/100"),
    geneID = c("TP53/BRCA1", "AKT1/MTOR"),
    stringsAsFactors = FALSE
  )
  res <- tcm_interpret(df, type = "enrichment")
  expect_s3_class(res, "tcm_ai_analysis")
  expect_equal(res$input$type, "enrichment")
})

test_that("tcm_interpret text mode passes audience and role to system prompt", {
  skip_if_not_installed("aisdk")
  captured <- list()
  fake_agent <- list(
    run = function(task, model) {
      captured$task <<- task
      list(text = "mock text")
    }
  )
  local_mocked_bindings(
    create_agent = function(name, description, system_prompt, ...) {
      captured$system <<- system_prompt
      fake_agent
    },
    get_model = function() "mock-model",
    .package = "aisdk"
  )
  out <- tcm_interpret(
    "test query",
    audience = "wetlab",
    role     = "You are a pharmacologist.",
    prompt   = "Focus on TCM:",
    verbose  = FALSE
  )
  expect_type(out, "character")
  expect_true(grepl("pharmacologist",  captured$system))
  expect_true(grepl("wet-lab",         captured$system))
  expect_true(grepl("Focus on TCM",    captured$task))
})

test_that("custom is absent from provider whitelist", {
  providers <- .available_providers()
  expect_false("custom" %in% providers)
  expect_true("openai"    %in% providers)
  expect_true("anthropic" %in% providers)
  expect_true("gemini"    %in% providers)
})

# ── Regression: artifacts, tools, routing ───────────────────────────────────

test_that("artifact registry keeps artifact metadata", {
  clear_tcm_artifacts()
  on.exit(clear_tcm_artifacts(), add = TRUE)

  handle <- save_tcm_artifact(
    object = data.frame(gene = "TP53", score = 1, stringsAsFactors = FALSE),
    artifact_type = "table_result",
    summary = "One-row test table."
  )

  art_df <- list_tcm_artifacts()
  row <- art_df[art_df$artifact_id == handle$artifact_id, , drop = FALSE]

  expect_equal(nrow(row), 1)
  expect_equal(row$artifact_type, "table_result")
  expect_equal(row$r_class, "data.frame")
})

test_that("create_tcm_tools exposes expanded tool modules", {
  skip_if_not_installed("aisdk")

  tools <- create_tcm_tools()
  tool_names <- vapply(tools, function(tool) tool$name, character(1))

  expect_true(all(c(
    "run_go_enrichment",
    "run_kegg_enrichment",
    "get_ppi_network",
    "prepare_ml_dataset",
    "run_ml_screening",
    "get_pubmed_evidence",
    "resolve_compound_cid",
    "plot_enrichment_result"
  ) %in% tool_names))
})

test_that("route_tcm_task recognizes new module categories", {
  pubmed_route <- route_tcm_task("Retrieve PubMed evidence for ginseng and diabetes")
  compound_route <- route_tcm_task("Resolve the PubChem CID for aspirin")
  visualization_route <- route_tcm_task("Plot a heatmap of stored metrics")

  expect_equal(pubmed_route$task_type, "pubmed")
  expect_equal(compound_route$task_type, "compound")
  expect_equal(visualization_route$task_type, "visualization")
})

test_that("get_ppi_network tool uses the local get_ppi wrapper", {
  skip_if_not_installed("aisdk")
  skip_if_not_installed("igraph")

  clear_tcm_artifacts()
  on.exit(clear_tcm_artifacts(), add = TRUE)

  mock_graph <- igraph::graph_from_data_frame(
    data.frame(from = "TP53", to = "BRCA1", score = 0.9),
    directed = FALSE
  )

  local_mocked_bindings(
    get_ppi = function(x, taxID = 9606, ...) mock_graph,
    .package = "TCMDATA"
  )

  result <- tool_get_ppi_network()$run(list(genes = c("TP53", "BRCA1"), tax_id = 9606L))

  expect_true(isTRUE(result$ok))
  expect_equal(result$artifact_type, "ppi_graph")
  expect_true(artifact_exists(result$artifact_id))
})

test_that("plot_ml_result tool supports upset plots", {
  skip_if_not_installed("aisdk")

  clear_tcm_artifacts()
  on.exit(clear_tcm_artifacts(), add = TRUE)

  mock_ml_1 <- .new_tcm_ml(
    method = "lasso",
    model = NULL,
    importance = data.frame(gene = c("TP53", "BRCA1"), importance = c(1, 0.8)),
    selected_features = c("TP53", "BRCA1"),
    cv_performance = list(auc = 0.91, sensitivity = 0.8, specificity = 0.9),
    ml_data = list()
  )
  mock_ml_1$genes <- mock_ml_1$selected_features

  mock_ml_2 <- .new_tcm_ml(
    method = "rf",
    model = NULL,
    importance = data.frame(gene = c("TP53", "AKT1"), importance = c(0.9, 0.7)),
    selected_features = c("TP53", "AKT1"),
    cv_performance = list(auc = 0.88, sensitivity = 0.75, specificity = 0.85),
    ml_data = list()
  )
  mock_ml_2$genes <- mock_ml_2$selected_features

  ml_list <- create_tcm_ml_list(lasso = mock_ml_1, rf = mock_ml_2)
  handle <- save_tcm_artifact(ml_list, artifact_type = "ml_result")
  fake_plot <- structure(list(), class = c("a_upset_plot", "list"))

  local_mocked_bindings(
    upsetplot = function(list, ...) fake_plot,
    .package = "TCMDATA"
  )

  result <- tool_plot_ml_result()$run(list(artifact_id = handle$artifact_id, plot_type = "upset"))

  expect_true(isTRUE(result$ok))
  expect_equal(result$artifact_type, "plot")
  expect_true(artifact_exists(result$artifact_id))
})
