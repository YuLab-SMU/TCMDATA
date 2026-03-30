# ============================================================
# TCMDATA AI Module — Complete Demo
# ============================================================
# This script demonstrates every AI capability in TCMDATA.
# Run line by line in RStudio for best experience.
# ============================================================

library(TCMDATA)

# ── 1. Configuration (one-time) ─────────────────────────────
# Write credentials to .env (only need to do this once per machine).
# Set environment variables before running:
#   Sys.setenv(TCM_API_KEY  = "sk-...")
#   Sys.setenv(TCM_PROVIDER = "openai")   # or anthropic, gemini, deepseek, ...
#   Sys.setenv(TCM_MODEL    = "gpt-4o-mini")
#   Sys.setenv(TCM_BASE_URL = "")         # optional proxy base URL

api_key <- Sys.getenv("TCM_API_KEY", "")
if (!nzchar(api_key)) stop("Set TCM_API_KEY before running this demo.")

tcm_config(
  provider = Sys.getenv("TCM_PROVIDER", "openai"),
  api_key  = api_key,
  model    = Sys.getenv("TCM_MODEL",    "gpt-4o-mini"),
  base_url = Sys.getenv("TCM_BASE_URL", "")
)

# ── 2. Setup (every session) ────────────────────────────────
tcm_setup()  # reads .env, initialises model

# ============================================================
# Layer 1: Structured Interpretation (interpret_* family)
# ============================================================

# ── 3a. Enrichment Analysis ─────────────────────────────────
enrich_df <- data.frame(
  ID = c("GO:0006954", "GO:0007165", "hsa04151",
         "GO:0042981", "GO:0008284"),
  p.adjust = c(1e-4, 8e-4, 2e-3, 8e-3, 3e-2),
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

# Default usage
res1 <- tcm_interpret(enrich_df, type = "enrichment")
print(res1)

# ── 3b. Table (e.g. DEG results) ────────────────────────────
deg_df <- data.frame(
  gene = c("IL6", "TNF", "AKT1", "TP53", "VEGFA",
           "EGFR", "MYC", "CASP3", "BCL2", "MTOR"),
  logFC = c(3.2, 2.8, -1.5, -2.1, 2.5,
            1.9, -1.8, 1.2, -2.3, -1.1),
  p.adjust = c(1e-6, 5e-6, 1e-4, 2e-4, 3e-4,
               5e-4, 8e-4, 1e-3, 2e-3, 5e-3),
  stringsAsFactors = FALSE
)

# Generic table interpretation — rows serialised as key=value pairs
res2 <- tcm_interpret(deg_df)
print(res2)

# ── 3c. Generic Table ───────────────────────────────────────
# Any data.frame works — rows are serialised as key=value pairs
summary_df <- data.frame(
  pathway = c("PI3K-Akt", "MAPK", "TNF signaling"),
  gene_count = c(15, 12, 8),
  avg_logFC = c(1.8, 2.1, 2.5),
  stringsAsFactors = FALSE
)

res3 <- tcm_interpret(summary_df)
print(res3)

# ============================================================
# Prompt Engineering: audience / prompt / system
# ============================================================

# ── 4a. audience — control tone with one word ───────────────
# Presets
tcm_interpret(enrich_df, type = "enrichment",
              audience = "wetlab")

# Free-text audience description
tcm_interpret(deg_df,
              audience = "临床医生，不懂生物信息学")

tcm_interpret(deg_df,
              audience = "本科毕业论文，需要详细解释每个术语")

# ── 4b. prompt — control what to focus on ────────────────────
tcm_interpret(enrich_df, type = "enrichment",
  prompt = "这是槲皮素靶点的富集结果，请重点关注抗炎通路：")

tcm_interpret(deg_df,
  prompt = "从这些差异基因中找出潜在的药物靶点：")

# ── 4c. system — full control over AI role ──────────────────
tcm_interpret(enrich_df, type = "enrichment",
  system = paste(
    "你是中医药理学专家，专注于草药-靶点-通路关系。",
    "用中文回答，重点分析中药活性成分的作用机制。",
    "必须返回 JSON 格式。"
  ))

# ── 4d. Combine audience + prompt ────────────────────────────
tcm_interpret(deg_df,
  audience = "中医药方向的审稿人",
  prompt = "请从网络药理学角度评估这些差异基因的药理学意义：")

# ============================================================
# Draft: Generate publication-ready paragraphs
# ============================================================

# ── 5. From interpret result → paper paragraph ──────────────
ai_res <- tcm_interpret(enrich_df, type = "enrichment")
draft  <- draft_result_paragraph(ai_res, language = "zh")
print(draft)

# Extract plain text
cat(as.character(draft))

# One-liner pipe
enrich_df |>
  tcm_interpret(type = "enrichment") |>
  draft_result_paragraph() |>
  as.character() |>
  cat()

# ============================================================
# Layer 2: Free-text Interpretation (agent-based)
# ============================================================

# ── 6a. Single query ─────────────────────────────────────────
tcm_interpret("IL6 logFC=3.2, TNF logFC=2.8, AKT1 logFC=-1.5")

# ── 6b. Batch queries ───────────────────────────────────────
tcm_interpret(c(
  "HALLMARK_INFLAMMATORY_RESPONSE enriched, p.adjust=1e-6",
  "AKT1 degree=42, betweenness=0.35 in PPI network",
  "quercetin targets: AKT1, TNF, IL6, VEGFA"
))

# ── 6c. Custom system for text ──────────────────────────────
tcm_interpret(
  "槲皮素的主要靶点是AKT1、TNF和IL6，这些靶点在PI3K-Akt通路中富集",
  system = "你是中医药理学教授，请用通俗易懂的语言解释给本科生听。"
)

# ============================================================
# Layer 3: Custom Agent (advanced)
# ============================================================

# ── 7a. Create a reusable agent ──────────────────────────────
tcm_agent <- create_tcm_agent(
  name = "TCMExpert",
  description = "TCM pharmacology interpreter",
  system_prompt = paste(
    "你是中医药网络药理学专家。",
    "回答要结合中医理论和现代药理学。",
    "用中文回答，3-5句话。"
  )
)

run_tcm_agent(tcm_agent,
  "黄芪中的黄芪甲苷(astragaloside IV)如何通过PI3K-Akt通路发挥心肌保护作用？")

# ── 7b. Batch with sapply ────────────────────────────────────
questions <- c(
  "槲皮素的抗炎机制是什么？",
  "丹参酮IIA如何抑制肿瘤血管生成？",
  "小檗碱降糖的分子靶点有哪些？"
)
answers <- sapply(questions, \(q) run_tcm_agent(tcm_agent, q, verbose = FALSE))

# ============================================================
# Switching providers (no code changes needed)
# ============================================================

# Switch to Anthropic Claude
# tcm_config("anthropic", "sk-ant-xxx", "claude-3-5-haiku-20241022")
# tcm_setup()

# Switch to Google Gemini
# tcm_config("gemini", "AIza-xxx", "gemini-2.0-flash")
# tcm_setup()

# Switch to DeepSeek
# tcm_config("deepseek", "sk-xxx", "deepseek-chat")
# tcm_setup()

# Temporary override (doesn't touch .env)
# tcm_setup("openai", api_key = "sk-xxx", model = "gpt-4o")
