---
name: analysis-preferences
description: >
  User analysis preferences and default parameters for TCM network pharmacology.
  This skill defines default behaviors, parameter choices, reporting style, and
  quality standards. The agent should consult these preferences before executing
  any analysis to ensure consistent, user-aligned results.
  Active on every turn as a background constraint layer.
---

# Analysis Preferences

## Purpose

This skill acts as a persistent "taste vector" — a set of user-defined defaults and stylistic preferences that shape how the agent conducts analyses. Rather than asking the user to specify parameters every time, the agent should apply these defaults unless the user explicitly overrides them.

## Default Parameters

### Search
- Herb name type: `Herb_cn_name` (Chinese name) preferred; fall back to `Herb_pinyin_name`
- Always confirm herb name spelling before searching (e.g., "黄芪" not "黄氏")

### Enrichment
- `pvalueCutoff = 0.05`
- `qvalueCutoff = 0.2`
- GO ontology: run all three (BP, CC, MF) by default, unless user specifies one
- KEGG: organism = "hsa" (human), readable = TRUE

### PPI Network
- `score_threshold = 400` (medium confidence) as default
- If user says "high confidence", use 700; if "highest", use 900
- Always compute topological metrics after PPI retrieval
- Rank by integrated hub score, report top 10 genes by default

### Machine Learning
- Default methods: `c("lasso", "rf", "svm", "xgboost")`
- Default consensus threshold: `min_methods = 2`
- Always set seed (default 42) for reproducibility
- Report AUC, sensitivity, specificity for each method

### Visualization
- Enrichment: show top 20 terms by default
- PPI: label top 10 hub genes
- Default plot dimensions: 10 x 8 inches at 300 DPI

## Reporting Style

### Language
- Detect user language and respond in kind (Chinese input → Chinese response)
- Keep gene symbols, pathway IDs, and compound names in their original English form
- Section headers in user's language, technical terms in English

### Structure
- For complete workflow results, organize as:
  1. Target retrieval summary (herb targets / disease targets / intersection count)
  2. PPI & hub gene analysis (top genes with metrics)
  3. Functional enrichment highlights (top GO terms + KEGG pathways)
  4. Key findings (2-3 bullet points of biological insight)
  5. Suggested next steps
- For single-tool results, keep it concise (3-5 sentences)

### Quantitative Standards
- Always cite exact numbers: gene counts, p-values, enrichment ratios
- Round p-values to 3 decimal places, percentages to 1 decimal
- When ranking, show the metric values alongside gene names

## Quality Standards

### Minimum Thresholds
- Intersection genes < 5: warn user, suggest broadening search
- Intersection genes > 500: warn user, suggest narrowing or filtering
- Enrichment with 0 significant terms: suggest parameter adjustment
- PPI with > 50% isolated nodes: suggest lowering score_threshold
- ML with < 2 methods: warn about weak consensus

### Grounding Rules
- Never claim a biological mechanism without tool output evidence
- Always attribute findings to specific artifacts (cite artifact_id)
- When interpreting results, distinguish between "statistically significant" and "biologically relevant"
- If literature validation returns 0 PubMed hits, explicitly note the lack of prior evidence

## Override Behavior

Users can override any preference by stating it explicitly:
- "Use p < 0.01" → override pvalueCutoff
- "Only run LASSO" → override default ML methods
- "Show top 30 terms" → override visualization defaults

When the user overrides a default, apply the override for the current turn only. Do not permanently change preferences unless the user says "always use..." or "from now on...".
