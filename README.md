# TCMDATA: Traditional Chinese Medicine Data Analysis and Visualization

<!-- badges: start -->
[![R-CMD-check](https://img.shields.io/badge/R--CMD--check-passing-brightgreen?logo=r)](https://github.com/Hinna0818/TCMDATA)
[![bookdown](https://github.com/Hinna0818/TCMDATA/actions/workflows/bookdown.yaml/badge.svg)](https://github.com/Hinna0818/TCMDATA/actions/workflows/bookdown.yaml)
[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic--2.0-blue.svg)](https://opensource.org/license/artistic-2-0)
[![R version](https://img.shields.io/badge/R%20≥-3.5.0-276DC3?logo=r)](https://cran.r-project.org/)
[![Platform](https://img.shields.io/badge/Platform-Linux%20|%20macOS%20|%20Windows-informational)](https://github.com/Hinna0818/TCMDATA)
[![Documentation](https://img.shields.io/badge/docs-bookdown-orange?logo=bookstack)](https://hinna0818.github.io/TCMDATA/)
<!-- badges: end -->

**TCMDATA** is an integrated R package for Traditional Chinese Medicine (TCM)
network pharmacology research. It provides a unified computational framework
for herb–molecule–target data retrieval, pharmacological network construction,
enrichment analysis, protein–protein interaction (PPI) analysis, and
publication-ready visualization.

## Features

| Module | Description | Key Functions |
|--------|-------------|---------------|
| **Data retrieval** | Bidirectional query of 500+ herbs, 10 000+ compounds, and validated targets | `search_herb()`, `search_target()` |
| **Molecule detection** | PubChem-based compound identifier resolution, property annotation, similarity search, structure download, and format conversion | `resolve_cid()`, `getprops()`, `compound_similarity()`, `download_ligand_structure()` |
| **Network construction** | Build herb–molecule–target networks with topological metrics | `prepare_herb_graph()` |
| **Enrichment analysis** | ORA using herbs as functional categories; GO/KEGG compatible | `herb_enricher()` |
| **PPI analysis** | Filtering, 15+ centrality metrics, community detection, and PPI robustness analysis | `ppi_subset()`, `compute_nodeinfo()`, `rank_ppi_nodes()` |
| **Clustering** | Louvain, MCL, and MCODE algorithms | `run_louvain()`, `run_MCL()`, `runMCODE()` |
| **ML screening** | LASSO, Elastic Net, Ridge, Random Forest + Boruta, SVM-RFE, XGBoost; three validation modes (A/B/C) with consensus analysis | `run_ml_screening()`, `plot_ml_roc()`, `plot_ml_venn()` |
| **Data interpretation** | AI-assisted interpretation for text, enrichment objects, PPI objects, tables, and result drafting | `tcm_interpret()`, `draft_result_paragraph()`, `tcm_interpret_schema()` |
| **Visualization** | Sankey, docking heatmap, PPI heatmap, network plots | `tcm_sankey()`, `ggdock()`, `plot_node_heatmap()` |
| **Other resources** | Supplementary datasets: gut microbiota–metabolite associations (gutMGene) and transcription factor–target regulation pairs | `gutMGene`, `tf_targets`, `dn_targets` |


## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
options(timeout = 600)
devtools::install_github("YuLab-SMU/TCMDATA")
```

## AI Quick Start

`TCMDATA` includes an AI interpretation module built on top of
[aisdk](https://github.com/YuLab-SMU/aisdk).

```r
# One-time install
devtools::install_github("YuLab-SMU/aisdk")

# Configure and initialise the model
tcm_setup(
  provider = "openai",
  api_key = "sk-xxxx",
  model = "gpt-5",
  base_url= "xxx", # if neccessary
  save = TRUE,
  test = TRUE
)
```

Object interpretation:

```r
library(clusterProfiler)

# Example: build a GO enrichment result from Lingzhi targets
lz_targets <- search_herb("lingzhi", "Herb_pinyin_name")$target
lz_targets <- unique(na.omit(sample(lz_targets, 100)))

enrich_res <- enrichGO(
  gene = lz_targets,
  OrgDb = 'org.Hs.eg.db',
  keyType = "SYMBOL",
  ont = "BP"
)

ai_res <- tcm_interpret(
  enrich_res,
  prompt = "This is a GO BP enrichment result derived from the potential targets of Lingzhi. Please briefly summarise its core biological significance.",
  language = "en"
)
print(ai_res)

draft <- draft_result_paragraph(ai_res, language = "en")
cat(as.character(draft))
```

Free-text interpretation:

```r
txt <- tcm_interpret(
  "Please introduce the major pharmacological functions of Huangqi (Astragalus membranaceus), and briefly explain its potential roles in immunoregulation and disease treatment.",
  verbose = FALSE,
  language = "en"
)
cat(txt)
```

Custom structured output:

```r
my_schema <- tcm_schema(
  summary = tcm_field_string("Brief summary"),
  key_targets = tcm_field_array("Top targets")
)

res <- tcm_interpret_schema(enrich_res, schema = my_schema, language = "en")
print(res)
```

## Documentation

Full documentation with worked examples is available at **[here](https://hinna0818.github.io/TCMDATA/)**.
