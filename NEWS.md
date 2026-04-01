# TCMDATA 0.1.0

*Released: 2026-04-01*

This is the first public release of **TCMDATA**, a comprehensive R toolkit for Traditional Chinese Medicine (TCM) network pharmacology research.

## 🌿 Data Retrieval

- `search_herb()`: Query herb–compound–target associations by herb name (Chinese, Pinyin, or Latin)
- `search_target()`: Reverse lookup—find herbs targeting specific genes

## 🔬 Molecule Detection (PubChem Integration)

- `resolve_cid()`: Resolve compound names to PubChem CIDs
- `getprops()`: Retrieve compound properties (MW, LogP, TPSA, etc.)
- `compound_similarity()`: Tanimoto similarity search against PubChem
- `download_ligand_structure()`: Download molecular structures (SDF/MOL2/PDB)
- `convert_structure()`: Convert between molecular file formats

## 📊 Network Analysis

- `prepare_herb_graph()`: Build herb–compound–target networks
- `ppi_subset()`: Filter PPI networks by confidence score
- `compute_nodeinfo()`: Calculate 15+ centrality metrics (degree, betweenness, closeness, MCC, DMNC, etc.)
- `rank_ppi_nodes()`: Rank nodes by multiple topological features
- `ppi_knock()`: Evaluate network robustness via drug-attack simulation

## 🧬 Clustering Algorithms

- `run_louvain()`: Louvain community detection
- `run_MCL()`: Markov Clustering (MCL)
- `run_mcode()`: MCODE algorithm for dense subgraph detection

## 📈 Enrichment Analysis

- `herb_enricher()`: Herb-based over-representation analysis (ORA)
- Compatible with GO/KEGG via clusterProfiler

## 🤖 Machine Learning Screening

- `run_ml_screening()`: Run 6 ML algorithms with 3 validation modes
  - LASSO, Elastic Net, Ridge

  - Random Forest + Boruta
  - SVM-RFE
  - XGBoost
- `get_ml_consensus()`: Consensus feature selection across methods
- `plot_ml_roc()`, `plot_ml_venn()`: Visualize ML results

## 💡 AI-Powered Interpretation

- `tcm_setup()`: Configure LLM provider (OpenAI, etc.)
- `tcm_interpret()`: Interpret enrichment results, PPI data, or free text
- `draft_result_paragraph()`: Generate manuscript-ready result paragraphs
- `tcm_interpret_schema()`: Structured output with custom schemas
- `create_tcm_agent()`, `make_tcm_function()`: Build reusable AI functions

## 📚 Literature Mining

- `get_pubmed_data()`: Search PubMed for TCM–disease literature
- `get_pubmed_table()`: Export results as publication tables

## 📊 Visualization

- `tcm_sankey()`: Herb–compound–target Sankey diagrams
- `ggdot_sankey()`: Combined dot-Sankey plots for enrichment
- `ggdock()`: Molecular docking affinity heatmaps
- `gglollipop()`: Lollipop plots for enrichment results
- `plot_node_heatmap()`: PPI centrality heatmaps
- `radar_plot()`: Radar charts for multi-dimensional comparisons
- `go_barplot()`, `gocircle_plot()`: GO enrichment visualizations
- `ggvenn_plot()`: Venn diagrams

## 📦 Bundled Datasets
 
- `gutMGene`: Gut microbiota–metabolite associations
- `tf_targets`: Transcription factor–target regulation pairs
- `dn_gcds`, `dn_otp`: Diabetic nephropathy reference datasets
- `demo_ppi`: Example PPI network for tutorials
