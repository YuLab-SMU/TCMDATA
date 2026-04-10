---
name: tcm-network-pharmacology
description: >
  Domain knowledge and tool usage reference for TCM network pharmacology.
  Consult this skill for parameter defaults, tool combinations, and quality
  thresholds. This is a REFERENCE — not an execution plan.
  Execute ONLY what the user asks for.
---

# TCM Network Pharmacology — Tool Reference

This document is a **look-up reference**, not an execution plan.
Read the section relevant to the user's current request; ignore the rest.

## Core Principle

**Do exactly what the user asks — nothing more, nothing less.**

- "查黄芪靶点" → call `search_herb_records`, return results. Done.
- "做个GO富集" → call `run_go_enrichment`. Done.
- "找出黄芪靶点，然后做GO富集" → do BOTH: search + enrichment in one turn.
- "做黄芪治疗脓毒症的网络药理学分析" → chain multiple tools as a full pipeline.

Keyword signals for multi-step: 然后, 接着, 并且, and then, comma-separated actions.
After completing ALL requested steps, you may mention 1–2 possible next steps but do NOT auto-execute them.

## Tool Quick Reference

| Task | Tool(s) | Key defaults |
|------|---------|-------------|
| Herb targets | `search_herb_records` | auto-detect type by input language |
| Disease targets | `search_disease_targets` | — |
| Target intersection | `compute_target_intersection` | supports 2/3/4-way |
| PPI network | `get_ppi_network` → `compute_ppi_metrics` → `rank_ppi_nodes` | score_threshold=400, top_n=10 |
| GO enrichment | `run_go_enrichment` | ontology="BP", p<0.05 |
| KEGG enrichment | `run_kegg_enrichment` | p<0.05 |
| Herb enrichment | `run_herb_enrichment` | — |
| ML screening | `prepare_ml_dataset` → `run_ml_screening` → `get_ml_consensus` | min_methods=2 |
| PubMed evidence | `get_pubmed_evidence` | max_results=20 |
| Compound info | `resolve_compound_cid` → `get_compound_properties` | — |
| GEO search | `search_geo_datasets` | organism="Homo sapiens" |
| Visualization | `plot_enrichment_result`, `plot_ppi_result`, `plot_herb_sankey`, `plot_ml_result`, `plot_pubmed_result`, `plot_docking_heatmap` | — |
| Interpretation | `interpret_artifact` | pass artifact_id |

## Parameter Defaults

- **Herb name type**: Chinese characters → `"Herb_cn_name"`, capitalized → `"Herb_pinyin_name"`, lowercase → `"Herb_en_name"`
- **PPI score**: 400 (medium). Use 700 for stringent.
- **Enrichment cutoffs**: pvalueCutoff=0.05, qvalueCutoff=0.2
- **Hub gene ranking**: composite hub_score from `rank_ppi_nodes`, top 10

## Quality Thresholds

- Intersection < 5 genes → warn user, suggest relaxing parameters
- Intersection > 500 genes → suggest narrowing disease scope
- PPI < 10 edges → try lowering score_threshold to 150
- Enrichment 0 significant terms → suggest relaxing cutoffs
- PubMed 0 results → report honestly, do NOT fabricate evidence

## Data Integrity Rules

- NEVER fabricate gene names, pathway names, or p-values
- ALWAYS cite exact numbers from tool outputs
- Distinguish evidence levels: clinical (HERB 2.0) > experimental > computational (TCMDATA) > single-database
- Use `generate_verification_urls` for cross-validation links when reporting key findings

## Multi-Step Pipeline (only when explicitly requested)

When the user asks for a "complete analysis", "network pharmacology", or "full pipeline":

1. **Targets**: herb targets + disease targets → intersection
2. **Network**: PPI → metrics → hub genes
3. **Enrichment**: GO-BP + KEGG on intersection genes
4. **Visualization**: sankey + PPI heatmap + enrichment lollipop
5. **Validation** (optional): PubMed evidence + cross-DB links
6. **Expression** (if data provided): WGCNA / ML / DEG integration
7. **Report**: structured summary with numbers from each step

Key rule for pipelines: always compute intersection BEFORE PPI/enrichment.
Use intersection genes (not herb or disease genes alone) for downstream analysis.

## Report Template (for full pipeline only)

```
## [Herb] treats [Disease] — Network Pharmacology

### Target Retrieval
- Herb targets: [N] (TCMDATA) | Disease targets: [N] (DisGeNET)
- Intersection: [N] common targets

### PPI & Hub Genes
- Network: [N] nodes, [M] edges (STRING ≥ [threshold])
- Top hub genes: [list]

### Enrichment
- GO-BP: [top 3] | KEGG: [top 3]

### Key Findings
[Grounded interpretation]
```
