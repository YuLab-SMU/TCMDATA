# End-to-end test: tcm_config -> tcm_setup -> tcm_interpret -> draft
# Run from project root: source("tests/test_ai_e2e.R")

devtools::load_all()

api_key <- Sys.getenv("TCM_API_KEY", "")
if (!nzchar(api_key)) {
  message("TCM_API_KEY not set — skipping live API tests.")
  message("Set via: Sys.setenv(TCM_API_KEY='sk-...')")
  quit(save = "no", status = 0)
}

cat("=== Step 1: tcm_config ===\n")
tcm_config(
  provider = Sys.getenv("TCM_PROVIDER", "openai"),
  api_key  = api_key,
  model    = Sys.getenv("TCM_MODEL",    "gpt-4o-mini"),
  base_url = Sys.getenv("TCM_BASE_URL", "")
)

cat("\n=== Step 2: tcm_setup ===\n")
tcm_setup()

cat("\n=== Step 3: mock enrichment data ===\n")
enrich_df <- data.frame(
  ID = c(
    "GO:0006954", "GO:0007165", "hsa04151",
    "GO:0042981", "GO:0008284"
  ),
  p.adjust = c(0.0001, 0.0008, 0.002, 0.008, 0.03),
  GeneRatio = c("12/200", "8/200", "6/200", "5/200", "4/200"),
  geneID = c(
    "IL6/TNF/IL1B/CXCL8/CCL2",
    "AKT1/MAPK1/PIK3CA/SRC/EGFR",
    "AKT1/MTOR/PIK3CA/PTEN/GSK3B",
    "BCL2/BAX/CASP3/CASP9/TP53",
    "CCND1/CDK4/MYC/PCNA/RB1"
  ),
  stringsAsFactors = FALSE
)
print(enrich_df)

cat("\n=== Step 4: tcm_interpret (auto-dispatch) ===\n")
ai_res <- tcm_interpret(enrich_df, type = "enrichment", language = "zh")
print(ai_res)

cat("\n=== Step 5: as.character() ===\n")
cat("Summary text:\n", as.character(ai_res), "\n")

cat("\n=== Step 6: draft_result_paragraph (pipe) ===\n")
draft <- draft_result_paragraph(ai_res, language = "zh")
print(draft)

cat("\n=== Step 7: as.character() on draft ===\n")
cat("Paragraph:\n", as.character(draft), "\n")

cat("\n=== All tests passed ===\n")
