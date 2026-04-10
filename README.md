# TCMDATA: Traditional Chinese Medicine Data Analysis and Visualization R Package

<img src="man/figures/logo.png" height="200" align="right" />

[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen?logo=r)](https://github.com/Hinna0818/TCMDATA)
[![bookdown](https://github.com/Hinna0818/TCMDATA/actions/workflows/bookdown.yaml/badge.svg)](https://github.com/Hinna0818/TCMDATA/actions/workflows/bookdown.yaml)
[![R version](https://img.shields.io/badge/R%20≥-4.1.0-276DC3?logo=r)](https://cran.r-project.org/)
[![Platform](https://img.shields.io/badge/Platform-Linux%20|%20macOS%20|%20Windows-informational)](https://github.com/Hinna0818/TCMDATA)
[![Documentation](https://img.shields.io/badge/docs-bookdown-orange?logo=bookstack)](https://hinna0818.github.io/TCMDATA/)

- **TCMDATA** is a comprehensive R toolkit for **Traditional Chinese Medicine (TCM) network pharmacology** research.
- It provides an end-to-end computational workflow—from herb–compound–target data retrieval and pharmacological network construction, through enrichment analysis and PPI mining, to machine-learning-based biomarker screening and AI-powered result interpretation.
- Publication-ready visualizations including Sankey diagrams, docking heatmaps, lollipop plots, and more.

For details, please visit the [full documentation](https://hinna0818.github.io/TCMDATA/).

---

<img src="man/figures/workflow.png" width="100%" />

## Highlights

- 🌿 **Built-in TCM Database**: Manually curated herb–compound–target interaction data ready for analysis
- 🔬 **PubChem Integration**: Compound identification, property retrieval, and structure download
- 📊 **PPI Network Analysis**: 17+ centrality metrics with community detection and robustness evaluation
- 🤖 **Machine Learning-based Feature Selection**: 6 algorithms (LASSO, Elastic Net, Ridge, RF+Boruta, SVM-RFE, XGBoost) with consensus scoring
- 🧠 **AI Agent System**: LLM-powered agent with natural language task routing, skill-based workflow, and interactive chat (`tcm_agent()`, `tcm_chat()`)
- 📚 **Literature & Data Mining**: PubMed search and GEO dataset discovery for TCM–disease studies
- 🔗 **Seamless Integration**: Works with clusterProfiler, enrichplot, and other Bioconductor tools for enrichment analysis

## Feature Overview

| Module | Description | Key Functions |
|--------|-------------|---------------|
| **Data Retrieval** | Bidirectional query of herbs, compounds, and validated targets | `search_herb()`, `search_target()` |
| **Disease Targets** | DisGeNET-based disease–gene queries and GEO dataset discovery | `search_disease()`, `search_gene_disease()`, `search_geo_datasets()` |
| **Molecule Detection** | PubChem-based CID resolution, property annotation, similarity search | `resolve_cid()`, `getprops()`, `compound_similarity()` |
| **Network Construction** | Build herb–compound–target networks with topological metrics | `prepare_herb_graph()` |
| **Enrichment Analysis** | Herb-based over-representation analysis; GO/KEGG compatible | `herb_enricher()` |
| **PPI Analysis** | STRING retrieval, 17+ centrality metrics, community detection, robustness | `get_ppi()`, `compute_nodeinfo()`, `ppi_knock()` |
| **Clustering** | Louvain, MCL, and MCODE community detection | `run_louvain()`, `run_MCL()`, `run_mcode()` |
| **ML Screening** | 6 algorithms × 3 validation modes with consensus analysis | `run_ml_screening()`, `plot_ml_roc()` |
| **AI Agent** | Natural language task routing, skill-based workflow, interactive chat | `tcm_agent()`, `tcm_chat()`, `tcm_interpret()` |
| **AI Interpretation** | LLM-powered interpretation for enrichment, PPI, tables | `tcm_interpret()`, `draft_result_paragraph()` |
| **Visualization** | Sankey, docking heatmaps, lollipop, UpSet, radar charts | `tcm_sankey()`, `ggdock()`, `gglollipop()`, `upsetplot()` |

## Installation

```r
# install.packages("devtools")
options(timeout = 600)
devtools::install_github("Hinna0818/TCMDATA")
```

## Quick Start

```r
library(TCMDATA)

# Search by herb name (supports Chinese pinyin)
huangqi <- search_herb("huangqi", "Herb_pinyin_name")
head(huangqi)

# Reverse lookup: find herbs targeting a specific gene
il6_herbs <- search_target("IL6")
head(il6_herbs)
```

## AI Agent System

TCMDATA integrates an AI agent layer via [aisdk](https://github.com/YuLab-SMU/aisdk), featuring natural language task routing, a skill-based workflow engine, and 36 registered tools.

```r
# One-time setup
devtools::install_github("YuLab-SMU/aisdk")

tcm_setup(
  provider = "openai",
  api_key  = "sk-xxxx",
  model    = "gpt-5-mini",
  save     = TRUE
)

# One-shot agent: auto-routes task to appropriate tools
result <- tcm_agent("Query the targets of Astragalus and list the top 10 most frequently occurring ones.")

# Interactive multi-turn chat session
tcm_chat()

# AI-powered result interpretation
ai_res <- tcm_interpret(enrich_res, language = "en")
draft  <- draft_result_paragraph(ai_res, language = "en")
```

The agent supports 11 task types including herb/disease/target lookup, enrichment, PPI analysis, ML screening, PubMed search, GEO dataset discovery, and more. Skills can be customized via `tcm_init_skills()`.

## Documentation

Complete tutorials with worked examples can be found [here](<https://hinna0818.github.io/TCMDATA/>).

## Citation

If you use TCMDATA in your research, please cite:

```
DOI pending
```
