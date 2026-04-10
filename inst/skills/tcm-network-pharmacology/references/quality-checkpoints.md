# Quality Checkpoints

Thresholds and sanity checks at each analysis stage. Warn or halt when values fall outside expected ranges.

## Stage 1: Target Retrieval

| Check | Expected | Action if Failed |
|-------|----------|-----------------|
| Herb found in TCMDATA | herb name matches ≥ 1 record | Suggest alternative spelling; try Pinyin/English name |
| Herb target count | 20-2000 | < 20: herb may have limited data; > 2000: normal for common herbs |
| Disease found in DisGeNET | ≥ 1 disease match | Try broader disease term, or use UMLS CUI |
| Disease gene count | 50-5000 | < 50: very specific disease; > 5000: very broad term |

## Stage 2: Intersection

| Check | Expected | Action if Failed |
|-------|----------|-----------------|
| Intersection size | 5-500 | < 5: check input, try broader terms; > 500: narrow disease |
| Intersection ratio | > 1% of smaller set | < 1%: weak herb-disease association, warn user |

## Stage 3: PPI Network

| Check | Expected | Action if Failed |
|-------|----------|-----------------|
| Connected nodes | > 50% of input genes | < 50%: many genes have no interactions at threshold |
| Edge count | > node count | < node count: sparse network, lower threshold to 150 |
| Largest component | > 60% of nodes | < 60%: fragmented, report largest component separately |

## Stage 4: Enrichment

| Check | Expected | Action if Failed |
|-------|----------|-----------------|
| Significant GO-BP terms | ≥ 5 at p.adjust < 0.05 | If 0: try p < 0.1 or check gene ID conversion |
| Significant KEGG pathways | ≥ 3 at p.adjust < 0.05 | If 0: try p < 0.1 |
| Background gene count | > 10000 | Much lower suggests gene ID mapping issue |

## Stage 5: Hub Genes

| Check | Expected | Action if Failed |
|-------|----------|-----------------|
| Top 10 hub genes selected | 10 genes (or all if < 10) | Report available count |
| Hub genes in herb targets | ≥ 50% overlap | Low overlap suggests indirect mechanism |
| Hub genes in disease targets | ≥ 50% overlap | Expected by construction |

## Stage 6: Literature

| Check | Expected | Action if Failed |
|-------|----------|-----------------|
| PubMed results | ≥ 1 | If 0: report "no direct literature evidence found", NOT a failure |
| Recent papers (last 5 years) | ≥ 30% of total | Low ratio: research is declining in this area |
