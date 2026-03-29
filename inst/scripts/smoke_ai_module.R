# tests/smoke_ai_module.R
# ─────────────────────────────────────────────────────────────────────────────
# Lightweight smoke test for the AI module.
# Does NOT require devtools::test() or a compiled package.
#
# Usage (from project root):
#   source("tests/smoke_ai_module.R")
#
# Useful when devtools::test() is blocked (e.g. OpenMP shared-memory errors
# in some environments). Mocks aisdk internals — no live API calls needed.
# ─────────────────────────────────────────────────────────────────────────────

stopifnot("Run from project root" = file.exists("DESCRIPTION"))

for (pkg in c("aisdk", "jsonlite", "igraph")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Need package: ", pkg)
}

message("Loading AI module source files...")
invisible(lapply(
  c("R/ai_core.R", "R/ai_schema.R", "R/ai_analysis.R",
    "R/ai_draft.R", "R/ai_agent.R",  "R/ai_print.R"),
  source
))

# ── Mocks ─────────────────────────────────────────────────────────────────────
mock_analysis_obj <- list(
  summary                   = "Immune signalling enriched.",
  key_findings              = list("IL6 and TNF upregulated"),
  biological_interpretation = "Inflammation pathway activated.",
  tcm_relevance             = "Quercetin targets IL6.",
  caveats                   = list()
)

assignInNamespace("generate_object",
  function(...) list(object = mock_analysis_obj),
  ns = "aisdk")

aisdk::set_model("mock-model")

# ── Helpers ───────────────────────────────────────────────────────────────────
pass <- function(msg) message("  PASS  ", msg)
fail <- function(msg) stop(  "  FAIL  ", msg, call. = FALSE)

# ── Test 1: compress_enrichment sorts by p.adjust before top_n ───────────────
{
  df <- data.frame(
    ID = c("A", "B", "C"),
    p.adjust = c(0.05, 0.001, 0.01),
    GeneRatio = c("2/100", "5/100", "3/100"),
    geneID = c("GX", "GY", "GZ"),
    stringsAsFactors = FALSE
  )
  ctx <- .compress_enrichment(df, top_n = 2)
  # After sorting by p.adjust asc, B (0.001) is row 1
  if (!grepl("^1\\. B", ctx))
    fail("compress_enrichment sort: expected B first (p.adjust=0.001)")
  pass("compress_enrichment sorts by p.adjust before top_n")
}

# ── Test 2: tcm_interpret type='enrichment' overrides data.frame dispatch ─────
{
  df <- data.frame(
    ID = c("GO:0001"),
    p.adjust = c(0.001),
    GeneRatio = c("3/100"),
    geneID = c("TP53/BRCA1"),
    stringsAsFactors = FALSE
  )
  res <- tcm_interpret(df, type = "enrichment")
  if (!inherits(res, "tcm_ai_analysis"))
    fail("dispatch: expected tcm_ai_analysis")
  if (res$input$type != "enrichment")
    fail(paste("dispatch: expected type='enrichment', got", res$input$type))
  pass("tcm_interpret type='enrichment' overrides data.frame class dispatch")
}

# ── Test 3: draft_result_paragraph consumes tcm_ai_analysis ──────────────────
{
  assignInNamespace("generate_object",
    function(...) list(object = list(
      paragraph           = "Mock publication paragraph.",
      figure_legend_hint  = "Mock legend."
    )),
    ns = "aisdk")

  ai_res <- .new_tcm_ai_analysis(
    input    = list(type = "enrichment", top_n = 10L,
                    max_genes = 5L, object_class = "data.frame"),
    context  = "1. GO:0006954 | p.adjust=1e-04",
    output   = mock_analysis_obj,
    metadata = .build_metadata("mock-model", "zh", "researcher", "data.frame")
  )

  draft <- draft_result_paragraph(ai_res)
  if (!inherits(draft, "tcm_ai_draft"))
    fail("draft: expected tcm_ai_draft class")
  if (!nzchar(draft$draft$paragraph))
    fail("draft: paragraph is empty")
  pass("draft_result_paragraph consumes tcm_ai_analysis")

  # Restore analysis mock for subsequent tests
  assignInNamespace("generate_object",
    function(...) list(object = mock_analysis_obj),
    ns = "aisdk")
}

# ── Test 4: text mode passes audience and role to system prompt ───────────────
{
  captured <- list()
  fake_agent <- list(
    run = function(task, model) {
      captured$task <<- task
      list(text = "mock text response")
    }
  )
  assignInNamespace("create_agent",
    function(name, description, system_prompt, ...) {
      captured$system <<- system_prompt
      fake_agent
    },
    ns = "aisdk")

  tcm_interpret(
    "test query",
    audience = "wetlab",
    role     = "You are a pharmacologist.",
    prompt   = "Focus on TCM relevance:",
    verbose  = FALSE
  )

  if (!grepl("pharmacologist", captured$system))
    fail("text mode: role not in system prompt")
  if (!grepl("wet-lab", captured$system))
    fail("text mode: audience 'wetlab' not reflected in system prompt")
  if (!grepl("Focus on TCM relevance", captured$task))
    fail("text mode: prompt prefix not prepended to task")
  pass("text mode passes role, audience, and prompt correctly")
}

# ── Test 5: custom not in provider whitelist ──────────────────────────────────
{
  pvs <- .available_providers()
  if ("custom" %in% pvs)
    fail("'custom' should not be in provider whitelist")
  for (p in c("openai", "anthropic", "gemini", "deepseek")) {
    if (!p %in% pvs)
      fail(paste("expected provider missing from whitelist:", p))
  }
  pass("provider whitelist is correct (custom absent, core providers present)")
}

message("\n5/5 smoke tests passed.")
