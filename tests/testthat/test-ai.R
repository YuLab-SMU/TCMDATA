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

test_that(".call_generate_object uses aisdk 1.5 tool-mode output", {
  captured <- list()
  local_mocked_bindings(
    generate_object = function(...) {
      captured$args <<- list(...)
      list(
        object = list(summary = "ok"),
        valid = TRUE,
        attempts = 1L
      )
    },
    .package = "aisdk"
  )
  old_options <- options(
    tcm.supports_native_tools = TRUE,
    tcm.force_json_schema = TRUE
  )
  on.exit(options(old_options), add = TRUE)

  result <- .call_generate_object(
    model = "mock-model",
    prompt = "prompt",
    schema = list(type = "object"),
    system = "system"
  )

  expect_identical(captured$args$mode, "tool")
  expect_identical(captured$args$schema_name, "tcm_result")
  expect_identical(captured$args$max_retries, 1L)
  expect_false("response_format" %in% names(captured$args))
  expect_identical(
    attr(result, "tcm_structured_output_mode", exact = TRUE),
    "tool"
  )
})

test_that(".call_generate_object retains JSON-only relay compatibility", {
  captured <- list()
  local_mocked_bindings(
    generate_object = function(...) {
      captured$args <<- list(...)
      list(object = list(summary = "ok"), valid = TRUE, attempts = 1L)
    },
    .package = "aisdk"
  )
  old_options <- options(
    tcm.supports_native_tools = FALSE,
    tcm.force_json_schema = TRUE
  )
  on.exit(options(old_options), add = TRUE)
  schema <- list(type = "object")

  result <- .call_generate_object(
    model = "mock-model",
    prompt = "prompt",
    schema = schema,
    system = "system"
  )

  expect_identical(captured$args$mode, "json")
  expect_identical(captured$args$response_format, schema)
  expect_identical(
    attr(result, "tcm_structured_output_mode", exact = TRUE),
    "json"
  )
})

test_that("invalid aisdk structured output is not accepted as valid", {
  result <- list(
    object = list(summary = "incomplete"),
    raw_text = "",
    valid = FALSE,
    attempts = 2L
  )

  extracted <- .extract_object(result, type = "analysis")

  expect_identical(extracted$output_mode, "fallback_text")
  expect_identical(extracted$output$summary, "")
})

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
  expect_match(ctx, "3 nodes total; showing top 2", fixed = TRUE)
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
  expect_match(ctx, "Table: 3 rows x 3 cols, showing top 3", fixed = TRUE)
})

test_that(".compress_table reports total rows before truncation", {
  df <- data.frame(gene = LETTERS[1:5], score = seq_len(5))

  ctx <- .compress_table(df, top_n = 2)

  expect_match(ctx, "Table: 5 rows x 2 cols, showing top 2", fixed = TRUE)
})

# ── Tests: workflow execution ───────────────────────────────────────────────

test_that("run_tcm_workflow uses the aisdk 1.5 tool API and retains trace", {
  fake_tools <- list(
    list(
      name = "first_tool",
      run = function(args) {
        list(
          ok = TRUE,
          artifact_id = "artifact_001",
          received = args
        )
      }
    ),
    list(
      name = "second_tool",
      run = function(args) {
        list(ok = TRUE, received = args, detail = list(value = 42L))
      }
    )
  )
  local_mocked_bindings(
    create_tcm_tools = function(...) fake_tools,
    .package = "TCMDATA"
  )

  workflow <- create_tcm_workflow(
    "trace_test",
    list(
      list(
        tool = "first_tool",
        params = list(genes = "{{genes}}")
      ),
      list(
        tool = "second_tool",
        params = list(
          artifact_id = "{{from_previous}}",
          label = "source-{{cohort}}"
        )
      )
    )
  )

  result <- run_tcm_workflow(
    workflow,
    genes = c("IL6", "CASP1"),
    cohort = "GSE154918",
    verbose = FALSE
  )
  trace <- attr(result, "tool_calls")

  expect_s3_class(result, "tcm_workflow_result")
  expect_length(result, 2L)
  expect_identical(result[[1L]]$received$genes, c("IL6", "CASP1"))
  expect_identical(result[[2L]]$received$artifact_id, "artifact_001")
  expect_identical(result[[2L]]$received$label, "source-GSE154918")
  expect_identical(attr(result, "status"), "completed")
  expect_length(trace, 2L)
  expect_identical(trace[[1L]]$result, result[[1L]])
  expect_identical(trace[[2L]]$arguments, result[[2L]]$received)
  expect_identical(trace[[2L]]$status, "completed")
  expect_true(is.numeric(trace[[2L]]$duration_seconds))
})

test_that("run_tcm_workflow retains the failed call and prior results", {
  fake_tools <- list(
    list(
      name = "successful_tool",
      run = function(args) list(ok = TRUE, artifact_id = "artifact_002")
    ),
    list(
      name = "failing_tool",
      run = function(args) stop("tool execution failed")
    )
  )
  local_mocked_bindings(
    create_tcm_tools = function(...) fake_tools,
    .package = "TCMDATA"
  )

  workflow <- create_tcm_workflow(
    "failure_test",
    list(
      list(tool = "successful_tool", params = list()),
      list(
        tool = "failing_tool",
        params = list(artifact_id = "{{from_previous}}")
      )
    )
  )

  expect_warning(
    result <- run_tcm_workflow(workflow, verbose = FALSE),
    "Step 2 failed"
  )
  trace <- attr(result, "tool_calls")

  expect_length(result, 2L)
  expect_true(result[[1L]]$ok)
  expect_false(result[[2L]]$ok)
  expect_match(result[[2L]]$error, "tool execution failed")
  expect_identical(attr(result, "status"), "failed")
  expect_identical(trace[[2L]]$arguments$artifact_id, "artifact_002")
  expect_identical(trace[[2L]]$status, "failed")
  expect_identical(trace[[2L]]$result, result[[2L]])
})

test_that("run_tcm_workflow executes a real aisdk 1.5 Tool object", {
  skip_if_not_installed("aisdk")

  real_tool <- aisdk::tool(
    name = "echo_tool",
    description = "Return the supplied value.",
    parameters = aisdk::z_object(
      value = aisdk::z_string(description = "Value to return")
    ),
    execute = function(args) list(ok = TRUE, value = args$value)
  )
  local_mocked_bindings(
    create_tcm_tools = function(...) list(real_tool),
    .package = "TCMDATA"
  )
  workflow <- create_tcm_workflow(
    "real_tool_test",
    list(list(tool = "echo_tool", params = list(value = "{{value}}")))
  )

  result <- run_tcm_workflow(workflow, value = "IL6", verbose = FALSE)

  expect_true(result[[1L]]$ok)
  expect_identical(result[[1L]]$value, "IL6")
  expect_identical(attr(result, "tool_calls")[[1L]]$status, "completed")
})

# ── Tests: analysis replay scripts ──────────────────────────────────────────

test_that("tool traces produce parseable replay scripts with dependencies", {
  generation_result <- list(
    all_tool_calls = list(
      list(
        id = "call_001",
        name = "search_herb_records",
        arguments = list(herb = "Huangqi", type = "Herb_pinyin_name")
      ),
      list(
        id = "call_002",
        name = "plot_herb_sankey",
        arguments = list(
          artifact_id = "search_001",
          axis_order = c("herb", "molecule", "target")
        )
      )
    ),
    all_tool_results = list(
      list(
        id = "call_001",
        raw_result = list(ok = TRUE, artifact_id = "search_001"),
        is_error = FALSE
      ),
      list(
        id = "call_002",
        raw_result = list(ok = TRUE, artifact_id = "plot_001"),
        is_error = FALSE
      )
    )
  )

  script <- .build_tcm_analysis_script(
    generation_result,
    task = "Search Huangqi and draw a Sankey plot",
    turn = 2L,
    model = "deepseek-chat"
  )
  code <- as.character(script)

  expect_s3_class(script, "tcm_analysis_script")
  expect_match(code, "step_01 <-", fixed = TRUE)
  expect_match(code, "artifact_id = step_01$artifact_id", fixed = TRUE)
  expect_match(code, "# Output artifact: plot_001", fixed = TRUE)
  expect_length(attr(script, "tool_calls"), 2L)
  expect_silent(parse(text = code))

  replay_tools <- list(
    list(
      name = "search_herb_records",
      run = function(args) list(
        ok = TRUE,
        artifact_id = "search_001",
        arguments = args
      )
    ),
    list(
      name = "plot_herb_sankey",
      run = function(args) list(
        ok = TRUE,
        artifact_id = "plot_001",
        arguments = args
      )
    )
  )
  replay_env <- new.env(parent = globalenv())
  replay_env$create_tcm_tools <- function(...) replay_tools

  expect_silent(eval(parse(text = code), envir = replay_env))
  expect_identical(replay_env$step_02$arguments$artifact_id, "search_001")
})

test_that("analysis scripts are exported to unique global-style names", {
  target_env <- new.env(parent = emptyenv())
  assign("tcm_script_001", "existing", envir = target_env)
  script <- structure(
    "library(TCMDATA)",
    class = c("tcm_analysis_script", "character")
  )

  script_name <- .export_tcm_analysis_script(script, envir = target_env)

  expect_identical(script_name, "tcm_script_002")
  expect_true(exists(script_name, envir = target_env, inherits = FALSE))
  exported <- get(script_name, envir = target_env, inherits = FALSE)
  expect_s3_class(exported, "tcm_analysis_script")
  expect_identical(attr(exported, "global_name"), script_name)
  expect_output(print(exported), "library\\(TCMDATA\\)")
})

test_that("analysis script generation ignores turns without tool calls", {
  expect_null(.build_tcm_analysis_script(
    list(text = "No analysis was requested."),
    task = "Hello"
  ))
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

test_that("create_tcm_task_agent disables skills with NULL for aisdk 1.5", {
  skip_if_not_installed("aisdk")

  captured <- list()
  local_mocked_bindings(
    create_agent = function(...) {
      captured$args <<- list(...)
      list(name = "mock-agent")
    },
    .package = "aisdk"
  )

  create_tcm_task_agent(
    tools = list(),
    system_prompt = "test prompt",
    skills = character(0)
  )

  expect_null(captured$args$skills)
})

test_that("the default agent prompt requires plot artifacts in .GlobalEnv", {
  prompt <- .default_tcm_system_prompt()

  expect_match(prompt, "PLOT OUTPUT RULE", fixed = TRUE)
  expect_match(prompt, ".GlobalEnv", fixed = TRUE)
  expect_match(prompt, "global_name", fixed = TRUE)
})

test_that("custom is absent from provider whitelist", {
  providers <- .available_providers()
  expect_false("custom" %in% providers)
  expect_true("openai"    %in% providers)
  expect_true("anthropic" %in% providers)
  expect_true("gemini"    %in% providers)
})

test_that("tcm_config persists relay protocol settings", {
  env_file <- tempfile(fileext = ".env")
  writeLines(c("OTHER_SETTING=keep", "TCM_API_FORMAT=chat_completions"), env_file)

  expect_message(
    tcm_config(
      provider = "deepseek",
      api_key = "test-key",
      model = "deepseek-chat",
      base_url = "https://relay.example.com/v1",
      path = env_file,
      api_format = "responses",
      supports_native_tools = FALSE,
      disable_stream_options = TRUE,
      responses_state_mode = "stateless"
    ),
    "saved"
  )

  config <- readLines(env_file)
  expect_true("OTHER_SETTING=keep" %in% config)
  expect_true("TCM_API_FORMAT=responses" %in% config)
  expect_true("TCM_SUPPORTS_NATIVE_TOOLS=false" %in% config)
  expect_true("TCM_DISABLE_STREAM_OPTIONS=true" %in% config)
  expect_true("TCM_RESPONSES_STATE_MODE=stateless" %in% config)
})

test_that("tcm_setup enables aisdk internet-check bypass by default", {
  skip_if_not_installed("aisdk")

  fake_provider <- list(
    language_model = function(model_id) {
      structure(
        list(model_id = model_id, provider = "openai"),
        class = "LanguageModelV1"
      )
    }
  )

  local_mocked_bindings(
    create_openai = function(...) fake_provider,
    set_model = function(model) invisible(NULL),
    .package = "aisdk"
  )

  old <- getOption("aisdk.skip_internet_check", NULL)
  on.exit(options(aisdk.skip_internet_check = old), add = TRUE)
  options(aisdk.skip_internet_check = FALSE)

  tcm_setup(
    provider = "openai",
    api_key = "test-key",
    model = "test-model",
    test = FALSE
  )

  expect_true(isTRUE(getOption("aisdk.skip_internet_check")))
})

test_that("tcm_setup test request omits temperature for proxy compatibility", {
  skip_if_not_installed("aisdk")

  fake_provider <- list(
    language_model = function(model_id) {
      structure(
        list(model_id = model_id, provider = "openai"),
        class = "LanguageModelV1"
      )
    }
  )
  captured <- list()

  local_mocked_bindings(
    create_openai = function(...) fake_provider,
    set_model = function(model) invisible(NULL),
    generate_text = function(...) {
      captured$args <<- list(...)
      list(text = "pong")
    },
    .package = "aisdk"
  )

  tcm_setup(
    provider = "openai",
    api_key = "test-key",
    model = "test-model",
    test = TRUE
  )

  expect_true("temperature" %in% names(captured$args))
  expect_null(captured$args$temperature)
})

test_that("tcm_setup resolves DeepSeek through the companion provider", {
  skip_if_not_installed("aisdk")

  captured <- list()
  fake_provider <- list(
    language_model = function(model_id) {
      structure(
        list(model_id = model_id, provider = "deepseek"),
        class = "LanguageModelV1"
      )
    }
  )

  local_mocked_bindings(
    .resolve_tcm_provider_factory = function(provider) {
      captured$provider <<- provider
      function(...) {
        captured$factory_args <<- list(...)
        fake_provider
      }
    },
    .package = "TCMDATA"
  )
  local_mocked_bindings(
    set_model = function(model) captured$model <<- model,
    .package = "aisdk"
  )

  model <- tcm_setup(
    provider = "deepseek",
    api_key = "test-key",
    model = "deepseek-chat",
    test = FALSE
  )

  expect_equal(captured$provider, "deepseek")
  expect_equal(captured$factory_args$api_key, "test-key")
  expect_equal(model$provider, "deepseek")
  expect_equal(captured$model$model_id, "deepseek-chat")
})

test_that("the installed companion package constructs a DeepSeek model", {
  skip_if_not_installed("aisdk.providers")

  create_deepseek <- .resolve_tcm_provider_factory("deepseek")
  expect_identical(
    environmentName(environment(create_deepseek)),
    "aisdk.providers"
  )

  provider <- create_deepseek(
    api_key = "test-key",
    base_url = "https://api.deepseek.com"
  )
  model <- provider$language_model("deepseek-chat")

  expect_true(inherits(model, "LanguageModelV1"))
  expect_equal(model$provider, "deepseek")
  expect_equal(model$model_id, "deepseek-chat")
})

test_that("tcm_setup configures OpenAI-compatible relay endpoints", {
  skip_if_not_installed("aisdk")

  captured <- list()
  fake_provider <- list(
    language_model = function(model_id) {
      structure(
        list(model_id = model_id, provider = "deepseek"),
        class = "LanguageModelV1"
      )
    }
  )

  local_mocked_bindings(
    create_custom_provider = function(...) {
      captured$provider_args <<- list(...)
      fake_provider
    },
    set_model = function(model) captured$model <<- model,
    .package = "aisdk"
  )

  tcm_setup(
    provider = "deepseek",
    api_key = "test-key",
    model = "deepseek-chat",
    base_url = "https://relay.example.com/v1",
    test = FALSE
  )

  expect_equal(captured$provider_args$provider_name, "deepseek")
  expect_equal(captured$provider_args$base_url, "https://relay.example.com/v1")
  expect_equal(captured$provider_args$api_format, "chat_completions")
  expect_true(captured$provider_args$supports_native_tools)
  expect_true(captured$provider_args$disable_stream_options)
  expect_equal(captured$provider_args$responses_state_mode, "stateless")
})

test_that("tcm_setup supports Responses-compatible relay endpoints", {
  skip_if_not_installed("aisdk")

  captured <- list()
  fake_provider <- list(
    language_model = function(model_id) {
      structure(
        list(model_id = model_id, provider = "openai"),
        class = "LanguageModelV1"
      )
    }
  )

  local_mocked_bindings(
    create_custom_provider = function(...) {
      captured$provider_args <<- list(...)
      fake_provider
    },
    set_model = function(model) invisible(NULL),
    .package = "aisdk"
  )

  tcm_setup(
    provider = "openai",
    api_key = "test-key",
    model = "gpt-test",
    base_url = "https://relay.example.com/v1",
    api_format = "responses",
    responses_state_mode = "auto",
    test = FALSE
  )

  expect_equal(captured$provider_args$api_format, "responses")
  expect_equal(captured$provider_args$responses_state_mode, "auto")
})

test_that("tcm_setup selects the native OpenAI Responses model explicitly", {
  skip_if_not_installed("aisdk")

  captured <- list()
  fake_provider <- list(
    language_model = function(model_id) {
      stop("Chat Completions model should not be selected")
    },
    responses_model = function(model_id) {
      captured$responses_model <<- model_id
      structure(
        list(model_id = model_id, provider = "openai"),
        class = "LanguageModelV1"
      )
    }
  )

  local_mocked_bindings(
    create_openai = function(...) {
      captured$provider_args <<- list(...)
      fake_provider
    },
    set_model = function(model) captured$model <<- model,
    .package = "aisdk"
  )

  tcm_setup(
    provider = "openai",
    api_key = "test-key",
    model = "gpt-test",
    api_format = "responses",
    responses_state_mode = "stateless",
    test = FALSE
  )

  expect_equal(captured$provider_args$api_format, "responses")
  expect_equal(captured$provider_args$responses_state_mode, "stateless")
  expect_equal(captured$responses_model, "gpt-test")
  expect_equal(captured$model$model_id, "gpt-test")
})

test_that("tcm_setup rejects a mismatched native provider protocol", {
  skip_if_not_installed("aisdk")

  expect_error(
    tcm_setup(
      provider = "gemini",
      api_key = "test-key",
      model = "gemini-test",
      api_format = "responses",
      test = FALSE
    ),
    "not supported by the native 'gemini' provider",
    fixed = TRUE
  )
})

test_that(".send_tcm_chat_turn falls back when streaming fails", {
  skip_if_not_installed("aisdk")

  failing_stream_session <- list(
    send_stream = function(prompt, callback, ...) {
      stop("stream failed")
    },
    get_last_response = function() ""
  )
  fallback_session <- list(
    send = function(prompt, ...) {
      captured$send_args <<- list(...)
      list(text = paste("batch", prompt), tool_calls = list())
    }
  )
  captured <- list()

  local_mocked_bindings(
    create_chat_session = function(agent, model) fallback_session,
    .package = "aisdk"
  )

  saw_stream_error <- FALSE
  out <- .send_tcm_chat_turn(
    chat_session = failing_stream_session,
    session_agent = list(),
    model = "mock-model",
    prompt = "hello",
    stream = TRUE,
    on_stream_error = function(e) {
      saw_stream_error <<- TRUE
    }
  )

  expect_true(saw_stream_error)
  expect_true(isTRUE(out$fallback))
  expect_equal(out$result$text, "batch hello")
  expect_identical(out$chat_session, fallback_session)
  expect_true("temperature" %in% names(captured$send_args))
  expect_null(captured$send_args$temperature)
})

test_that(".send_tcm_chat_turn omits temperature for batch chat", {
  captured <- list()
  batch_session <- list(
    send = function(prompt, ...) {
      captured$send_args <<- list(...)
      list(text = "batch ok", tool_calls = list())
    }
  )

  out <- .send_tcm_chat_turn(
    chat_session = batch_session,
    session_agent = list(),
    model = "mock-model",
    prompt = "hello",
    stream = FALSE
  )

  expect_false(isTRUE(out$fallback))
  expect_equal(out$result$text, "batch ok")
  expect_true("temperature" %in% names(captured$send_args))
  expect_null(captured$send_args$temperature)
})

test_that(".send_tcm_chat_turn preserves the complete streamed result", {
  complete_calls <- list(
    list(name = "first_tool", arguments = list(value = 1L)),
    list(name = "second_tool", arguments = list(value = 2L))
  )
  stream_result <- list(
    text = "stream ok",
    tool_calls = complete_calls[2L],
    all_tool_calls = complete_calls,
    steps = 3L
  )
  stream_session <- list(
    send_stream = function(prompt, callback, ...) {
      callback("stream ok", done = TRUE)
      stream_result
    },
    get_last_response = function() "stream ok"
  )

  out <- .send_tcm_chat_turn(
    chat_session = stream_session,
    session_agent = list(),
    model = "mock-model",
    prompt = "hello",
    stream = TRUE
  )

  expect_identical(out$result, stream_result)
  expect_length(.tcm_result_tool_calls(out$result), 2L)
  expect_identical(.tcm_result_tool_calls(out$result), complete_calls)
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

test_that("tool artifacts are immediately available to eval_r_code", {
  skip_if_not_installed("aisdk")

  clear_tcm_artifacts()
  on.exit(clear_tcm_artifacts(), add = TRUE)

  result <- .save_tool_artifact(
    object = data.frame(target = c("TP53", "AKT1"), stringsAsFactors = FALSE),
    artifact_type = "search_result",
    function_name = "unit_test"
  )
  on.exit(
    if (exists(result$artifact_id, envir = globalenv(), inherits = FALSE)) {
      rm(list = result$artifact_id, envir = globalenv())
    },
    add = TRUE
  )

  expect_true(exists(result$artifact_id, envir = globalenv(), inherits = FALSE))

  code <- sprintf(
    "df <- get('%s', envir = .GlobalEnv); cat(paste(unique(df$target), collapse = ','))",
    result$artifact_id
  )
  eval_result <- tool_eval_r_code()$run(list(code = code))

  expect_true(isTRUE(eval_result$ok))
  expect_true(grepl("TP53,AKT1", eval_result$output, fixed = TRUE))
})

test_that("plot artifacts report their .GlobalEnv object name", {
  clear_tcm_artifacts()
  on.exit(clear_tcm_artifacts(), add = TRUE)

  plot_obj <- structure(list(label = "test plot"), class = "test_plot")
  result <- .save_tool_artifact(
    object = plot_obj,
    artifact_type = "plot",
    function_name = "unit_test_plot"
  )
  on.exit(
    if (exists(result$artifact_id, envir = globalenv(), inherits = FALSE)) {
      rm(list = result$artifact_id, envir = globalenv())
    },
    add = TRUE
  )

  expect_equal(result$global_name, result$artifact_id)
  expect_equal(result$global_environment, ".GlobalEnv")
  expect_true(exists(result$global_name, envir = globalenv(), inherits = FALSE))
  expect_identical(get(result$global_name, envir = globalenv()), plot_obj)
})

test_that("create_tcm_tools exposes expanded tool modules", {
  skip_if_not_installed("aisdk")

  tools <- create_tcm_tools()
  tool_names <- vapply(tools, function(tool) tool$name, character(1))

  expect_true(all(c(
    "search_disease_targets",
    "run_go_enrichment",
    "run_kegg_enrichment",
    "get_ppi_network",
    "prepare_ml_dataset",
    "run_ml_screening",
    "get_pubmed_evidence",
    "resolve_compound_cid",
    "plot_enrichment_result"
  ) %in% tool_names))
  expect_false("search_gene_diseases" %in% tool_names)
})

test_that("route_tcm_task recognizes new module categories", {
  pubmed_route <- route_tcm_task("Retrieve PubMed evidence for ginseng and diabetes")
  compound_route <- route_tcm_task("Resolve the PubChem CID for aspirin")
  visualization_route <- route_tcm_task("Plot a heatmap of stored metrics")
  disease_route <- route_tcm_task("Retrieve sepsis-associated targets from Open Targets")

  expect_equal(pubmed_route$task_type, "pubmed")
  expect_equal(compound_route$task_type, "compound")
  expect_equal(visualization_route$task_type, "visualization")
  expect_equal(disease_route$task_type, "disease_lookup")
  expect_equal(disease_route$source_hint, "api")
  expect_true("search_disease_targets" %in% disease_route$tools)
  expect_false("search_gene_diseases" %in% disease_route$tools)
})

test_that("disease-target AI tool uses Open Targets wrapper", {
  skip_if_not_installed("aisdk")

  clear_tcm_artifacts()
  on.exit(clear_tcm_artifacts(), add = TRUE)

  local_mocked_bindings(
    query_disease_targets = function(disease_name, size, score_threshold) {
      data.frame(
        ensembl_id = "ENSG00000136244",
        gene_symbol = "IL6",
        gene_name = "interleukin 6",
        biotype = "protein_coding",
        score = 0.82,
        stringsAsFactors = FALSE
      )
    },
    .package = "TCMDATA"
  )

  result <- tool_query_disease_targets()$run(list(
    disease = "sepsis",
    size = 25,
    score_threshold = 0.2
  ))

  expect_true(isTRUE(result$ok))
  expect_match(result$summary, "Open Targets")
  expect_true(artifact_exists(result$artifact_id))
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
