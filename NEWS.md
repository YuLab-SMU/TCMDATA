# TCMDATA 0.1.0

*Released: 2026-04-01*

This is the first public release of **TCMDATA**, a comprehensive R toolkit for Traditional Chinese Medicine (TCM) network pharmacology research.

## 🌿 Data Retrieval

- `search_herb()`: Query herb–compound–target associations by herb name (Chinese, Pinyin, or Latin)
- `search_target()`: Reverse lookup—find herbs targeting specific genes
- `search_disease()`: Query DisGeNET (via DOSE) for disease-associated genes; supports UMLS CUI, exact name, and fuzzy matching
- `search_gene_disease()`: Reverse lookup—find all diseases associated with given gene symbols or Entrez IDs
- `search_geo_datasets()`: Query NCBI GEO for expression datasets (GDS/GSE) by disease keyword; filters by organism and dataset type

## 🔬 Molecule Detection (PubChem Integration)

- `resolve_cid()`: Resolve compound names to PubChem CIDs
- `getprops()`: Retrieve compound properties (MW, LogP, TPSA, etc.)
- `compound_similarity()`: Tanimoto similarity search against PubChem
- `download_ligand_structure()`: Download molecular structures (SDF/MOL2/PDB)
- `convert_structure()`: Convert between molecular file formats

## 📊 Network Analysis

- `get_ppi()`: Retrieve STRING PPI networks (wrapper around `clusterProfiler::getPPI()`)
- `prepare_herb_graph()`: Build herb–compound–target networks
- `ppi_subset()`: Filter PPI networks by confidence score
- `compute_nodeinfo()`: Calculate 17+ centrality metrics (degree, betweenness, closeness, MCC, DMNC, BN, EPC, radiality, Stress, etc.)
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

- `tcm_setup()`: Configure LLM provider (OpenAI, Anthropic, Gemini, Deepseek, aihubmix, and 10+ others)
- `tcm_interpret()`: Interpret enrichment results, PPI data, or free text
- `draft_result_paragraph()`: Generate manuscript-ready result paragraphs
- `tcm_interpret_schema()`: Structured output with custom schemas
- `create_tcm_agent()`, `make_tcm_function()`: Build reusable AI functions

## 🧠 AI Agent System (NEW)

- `tcm_agent()`: One-shot natural language entry point with automatic task routing to appropriate tools
- `tcm_chat()`: Interactive multi-turn chat session with built-in commands (`/help`, `/artifacts`, `/history`, `/model`, `/stats`, `/clear`, `/quit`)
- `route_tcm_task()`: Rule-based task router supporting 11 task types (herb_lookup, disease_lookup, enrichment, ppi_analysis, ml_screening, pubmed, geo_search, visualization, interpretation, etc.) with Chinese and English keyword patterns
- `create_tcm_task_agent()`: Create a configured aisdk agent with tools and skills for custom workflows
- `create_tcm_tools()`: Generate 36 aisdk Tool objects wrapping TCMDATA functions; supports filtering by task_type or tool_names
- `tcm_init_skills()`: Copy bundled skills to local directory for customization
- `tcm_use_skills()` / `tcm_reset_skills()` / `tcm_skill_dir()`: Skill directory management

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
- `upsetplot()`: UpSet plot for gene set intersection visualization (wrapper around `aplotExtra::upset_plot()`)

## 📦 Bundled Datasets
 
- `gutMGene`: Gut microbiota–metabolite associations
- `tf_targets`: Transcription factor–target regulation pairs
- `dn_gcds`, `dn_otp`: Diabetic nephropathy reference datasets (GeneCards & Open Targets)
- `deg_earlydn`: DESeq2 results for early diabetic nephropathy (GSE142025)
- `demo_ppi`: Example PPI network for tutorials
- `covid19`: COVID-19 case study dataset
