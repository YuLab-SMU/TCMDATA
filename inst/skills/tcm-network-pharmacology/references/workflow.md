# TCM Network Pharmacology — Unified Workflow

Complete analysis pipeline from herb/disease input through multi-omics validation.
Covers standard (database-only), expression-enhanced (RNA-seq/DEGs), WGCNA, single-cell, and ML workflows.

> **Scope rule**: This document describes the FULL pipeline. Only execute phases the user actually asks for. For a simple herb target query, run Phase 1 Step 1 only. For a single enrichment, run only that step. Execute all phases only when the user explicitly requests a complete or systematic analysis.

---

## Phase 0: Data Reconnaissance (Optional but Recommended)

Before starting analysis, help the user discover relevant public datasets.

### Step 0a: Search GEO for RNA-seq / Microarray Datasets

```
geo_result <- search_geo_datasets(
  disease = [disease_name],
  organism = "Homo sapiens",
  max_results = 20
)
```

Present results as a table: GSE accession, title, platform, sample count, summary.
Let the user choose which dataset(s) to use for downstream analysis.

**When to suggest this step:**
- User asks about a specific disease but hasn't provided expression data
- User mentions "GEO", "dataset", "transcriptome", or "expression profile"
- At the start of any analysis, proactively mention: "I found N relevant GEO datasets for [disease]. Would you like to incorporate expression data?"

### Step 0b: Search GEO for Single-Cell Datasets

```
geo_sc_result <- search_geo_datasets(
  disease = [disease_name],
  organism = "Homo sapiens",
  dataset_type = "single-cell",
  max_results = 10
)
```

Single-cell datasets can be used for cell-type-level validation in Phase 4.

---

## Phase 1: Target Collection

### Step 1: Herb Target Retrieval

```
herb_result <- search_herb_records(herb = [herb_names], type = [auto-detect])
herb_targets <- unique(herb_result$target)
```

Auto-detect type:
- Contains Chinese characters → `"Herb_cn_name"`
- Capitalized English → `"Herb_pinyin_name"`
- Lowercase English → `"Herb_en_name"`

If multiple herbs (a formula), combine all targets:
```
herb_result <- search_herb_records(herb = c("黄芪", "当归", "甘草"), type = "Herb_cn_name")
herb_targets <- unique(herb_result$target)
```

Report: "[N] targets found for [herb(s)] from TCMDATA database (source: TCMSP/SymMap)"

### Step 2: Disease Target Retrieval

```
disease_result <- search_disease_targets(disease = [disease_name])
disease_genes <- unique(disease_result$symbol)
```

For ambiguous disease names, try the most specific term first. If results are too few (< 50 genes), try broader terms.

Report: "[N] disease-associated genes found for [disease] from DisGeNET"

### Step 3: Target Intersection

**Mode A — Standard (two-way intersection):**
```
intersection <- compute_target_intersection(
  gene_lists = list(herb_targets, disease_genes),
  set_names = c("[Herb] Targets", "[Disease] Targets")
)
common_targets <- intersection$intersection
```

**Mode B — Expression-enhanced (three-way intersection, with DEGs):**
```
intersection <- compute_target_intersection(
  gene_lists = list(herb_targets, disease_genes, deg_list),
  set_names = c("[Herb] Targets", "[Disease] Genes", "DEGs")
)
common_targets <- intersection$intersection
```

**Mode C — WGCNA-enhanced (three-way intersection, with WGCNA hub genes):**
```
intersection <- compute_target_intersection(
  gene_lists = list(herb_targets, disease_genes, wgcna_hub_genes),
  set_names = c("[Herb] Targets", "[Disease] Genes", "WGCNA Hub Genes")
)
common_targets <- intersection$intersection
```

**Mode D — Multi-source (four-way intersection):**
When both DEGs and WGCNA results are available:
```
intersection <- compute_target_intersection(
  gene_lists = list(herb_targets, disease_genes, deg_list, wgcna_hub_genes),
  set_names = c("[Herb] Targets", "[Disease] Genes", "DEGs", "WGCNA Hub Genes")
)
common_targets <- intersection$intersection
```

**Decision points:**
- `length(common_targets) < 5`: Warn user. Suggest relaxing parameters or using two-way intersection instead.
- `length(common_targets) > 500`: Suggest narrowing disease scope or selecting specific herbs.
- `5 ≤ length(common_targets) ≤ 500`: Proceed normally.

Report: "[N] common targets at the intersection"

---

## Phase 2: Network & Enrichment Analysis

### Step 4: PPI Network Construction

```
ppi <- get_ppi_network(genes = common_targets, score_threshold = 400)
metrics <- compute_ppi_metrics(artifact_id = ppi$artifact_id)
hub_genes <- rank_ppi_nodes(artifact_id = metrics$artifact_id, top_n = 10)
```

If PPI returns very few edges (< 10), try lowering score_threshold to 150.

Report: "PPI network: [N] nodes, [M] edges; Top hub genes: [list]"

### Step 5: Functional Enrichment

Run both GO-BP and KEGG on the intersection genes (NOT on herb targets or disease targets alone):

```
go_result <- run_go_enrichment(genes = common_targets, ontology = "BP")
kegg_result <- run_kegg_enrichment(genes = common_targets)
```

Report top 5 significant terms (p.adjust < 0.05) for each.

### Step 6: Visualization

Generate publication-ready figures in this order:
1. `plot_herb_sankey` — Herb-Compound-Target network (use herb artifact)
2. `plot_ppi_result(type = "heatmap")` — Hub gene topology metrics
3. `plot_enrichment_result(type = "lollipop")` — GO/KEGG lollipop plot

---

## Phase 3: Expression-Based Validation (When Expression Data Available)

Skip this phase entirely if user has no expression data / DEG list / expression matrix.

### Step 7a: WGCNA Module Analysis

**When to use:** User has expression matrix with ≥ 15 samples per group.

WGCNA identifies co-expression modules correlated with the disease phenotype. The hub genes of key modules serve as an independent gene set for intersection.

Recommended workflow (user runs in R with guidance):
```r
# Recommend the user to perform WGCNA using the WGCNA package:
# 1. Construct weighted adjacency matrix
# 2. Identify modules via dynamic tree cut
# 3. Correlate modules with clinical traits
# 4. Extract hub genes from key module(s) (module membership > 0.8, gene significance > 0.2)

# The agent should guide these parameter choices:
# - softPower: use pickSoftThreshold() to determine
# - minModuleSize: 30 (default)
# - mergeCutHeight: 0.25 (default)
# - Module selection: choose module(s) with highest |correlation| to disease trait

# After WGCNA, user provides the hub gene list:
wgcna_hub_genes <- [user_provided_list]
```

**Integration with network pharmacology:**
- Use `wgcna_hub_genes` in Phase 1 Step 3 (Mode C or D) for intersection
- WGCNA hub genes that overlap with herb targets are high-confidence therapeutic targets
- Report module color, module-trait correlation, and number of hub genes extracted

### Step 7b: ML Feature Selection

**When to use:** User has expression matrix with group labels (e.g. disease vs control).

```
ml_data <- prepare_ml_dataset(
  expr_matrix = [expression_matrix],
  group = [group_labels],
  gene_filter = common_targets  # restrict to intersection genes
)
ml_result <- run_ml_screening(artifact_id = ml_data$artifact_id)
consensus <- get_ml_consensus(artifact_id = ml_result$artifact_id, min_methods = 2)
```

The consensus genes are the final "key targets" with both network pharmacology supporting and ML validation.

Add ML-specific plots:
```
plot_ml_result(type = "venn")   # ML method overlap
plot_ml_result(type = "roc")    # ROC curves for top genes
```

### Step 7c: DEG Analysis Integration

**When to use:** User provides a DEG list (e.g. from limma/DESeq2 on GEO data).

- Use DEG list directly in Phase 1 Step 3 (Mode B) for three-way intersection
- DEG direction (up/down) can further refine biological interpretation
- If user has full DESeq2/limma results, suggest filtering: |log2FC| > 1, adj.p < 0.05

---

## Phase 4: Single-Cell Validation (When scRNA-seq Data Available)

Skip this phase entirely if no single-cell data is available.

### Step 8: Single-Cell Analysis for Target Validation

**Purpose:** Validate hub genes / key targets at single-cell resolution — confirm cell-type-specific expression patterns.

**Recommended packages** (guide user to use in their R session):

1. **Seurat** (most common):
```r
# Standard Seurat workflow:
library(Seurat)
obj <- CreateSeuratObject(counts = counts_matrix)
obj <- NormalizeData(obj) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA()
obj <- FindNeighbors(obj) %>% FindClusters() %>% RunUMAP()

# Validate hub genes at cell-type level:
FeaturePlot(obj, features = hub_genes[1:4])
DotPlot(obj, features = hub_genes, group.by = "cell_type")
VlnPlot(obj, features = hub_genes[1:4], group.by = "cell_type")

# Differential expression in specific cell types:
markers <- FindMarkers(obj, ident.1 = "disease", ident.2 = "control",
                       subset.ident = "target_cell_type")
```

2. **SingleCellExperiment / Bioconductor ecosystem**:
```r
library(SingleCellExperiment)
library(scater)
library(scran)
sce <- SingleCellExperiment(assays = list(counts = counts_matrix))
sce <- logNormCounts(sce)
sce <- runPCA(sce) %>% runUMAP()

# Validate targets:
plotExpression(sce, features = hub_genes, x = "cell_type")
```

3. **sclet** (lightweight alternative):
Details can be found at: https://github.com/YuLab-SMU/sclet
```r
library(sclet)
# Follow sclet workflow for QC, normalization, clustering
# Validate hub gene expression across cell types
```

**Integration with network pharmacology results:**
- Check which cell types express the hub genes most strongly
- Verify that key pathway genes (from KEGG enrichment) are active in disease-relevant cell types
- Cell-type-specific expression adds biological depth to the network pharmacology findings
- Report: "Hub gene [X] is highly expressed in [cell_type], consistent with [pathway] enrichment"

---

## Phase 5: Validation & Reporting

### Step 9: Literature Validation

```
pubmed <- get_pubmed_evidence(herb = [herb_pinyin], disease = [disease_en], max_results = 20)
```

Summarize: how many papers, publication trend, key findings from abstracts.

### Step 10: External Cross-Validation

```
urls <- generate_verification_urls(herbs = [herb_names], targets = hub_genes, disease = [disease])
```

Provide HERB database and ShenNong Alpha links for user to cross-check key findings.

### Step 11: Compound-Target Molecular Docking (Optional)

If user requests docking or molecular validation:
```
cid <- resolve_compound_cid(identifier = [compound_name])
props <- get_compound_properties(cid = cid)
```

Report on druglikeness (Lipinski's Rule of Five from properties).

### Step 12: Interpretation & Report

Call `interpret_artifact` on:
1. The enrichment result (GO)
2. The PPI hub gene ranking
3. The KEGG result

Synthesize into the report template defined in SKILL.md.

---

## Workflow Decision Tree

```
User Input → Herb + Disease
  │
  ├─ Phase 0: Search GEO for relevant datasets (recommend to user)
  │
  ├─ Phase 1: Collect targets → Intersection
  │   ├─ No expression data → Two-way intersection (Mode A)
  │   ├─ Has DEGs → Three-way intersection (Mode B)
  │   ├─ Has WGCNA hub genes → Three-way intersection (Mode C)
  │   └─ Has DEGs + WGCNA → Four-way intersection (Mode D)
  │
  ├─ Phase 2: PPI + Enrichment (always)
  │
  ├─ Phase 3: Expression validation (if data available)
  │   ├─ WGCNA → Module hub genes
  │   ├─ ML screening → Consensus biomarkers
  │   └─ DEG integration → Direction-aware filtering
  │
  ├─ Phase 4: Single-cell validation (if scRNA-seq available)
  │   └─ Cell-type-specific expression of hub genes
  │
  └─ Phase 5: Literature + Cross-validation + Report
```

## Phase Selection Guide

| User has | Phases to run | Intersection mode |
|----------|--------------|-------------------|
| Herb + Disease only | 0→1→2→5 | Mode A (two-way) |
| + DEG list | 0→1→2→3c→5 | Mode B (three-way) |
| + Expression matrix | 0→1→2→3a→3b→5 | Mode C or D |
| + scRNA-seq data | 0→1→2→(3)→4→5 | Mode A/B/C/D |
| Full multi-omics | 0→1→2→3→4→5 | Mode D (four-way) |
