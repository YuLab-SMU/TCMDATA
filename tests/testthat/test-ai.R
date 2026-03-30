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
