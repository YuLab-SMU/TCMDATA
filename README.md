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
| **Network construction** | Build herb–molecule–target networks with topological metrics | `prepare_herb_graph()` |
| **Enrichment analysis** | ORA using herbs as functional categories; GO/KEGG compatible | `herb_enricher()` |
| **PPI analysis** | Filtering, 15+ centrality metrics, community detection, and PPI robustness analysis | `ppi_subset()`, `compute_nodeinfo()`, `rank_ppi_nodes()` |
| **Clustering** | Louvain, MCL, and MCODE algorithms | `run_louvain()`, `run_MCL()`, `runMCODE()` |
| **ML screening** | LASSO, Elastic Net, Ridge, Random Forest + Boruta, SVM-RFE, XGBoost; three validation modes (A/B/C) with consensus analysis | `run_ml_screening()`, `plot_ml_roc()`, `plot_ml_venn()` |
| **Visualization** | Sankey, docking heatmap, PPI heatmap, network plots | `tcm_sankey()`, `ggdock()`, `plot_node_heatmap()` |
| **Other resources** | Supplementary datasets: gut microbiota–metabolite associations (gutMGene) and transcription factor–target regulation pairs | `gutMGene`, `tf_targets`, `dn_targets` |

## Installation

Install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("YuLab-SMU/TCMDATA")
```

## Documentation

Full documentation with worked examples is available at **[here](https://hinna0818.github.io/TCMDATA/)**.


