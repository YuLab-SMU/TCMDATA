# Data Sources & Grounding Rules

This document describes the authoritative data sources used by TCMDATA tools. All analysis claims MUST be grounded in these sources. Never fabricate data.

## TCMDATA Herb-Compound-Target Database

- **Content**: 1,808 herbs (Chinese name), 2,166 compounds, 15,889 targets, 208,693 records
- **Original sources**: CTD, BATMAN-TCM, TM-MC, TCMIO, NPASS, ITCM
- **Schema**: `Herb_cn_name | Herb_en_name | Herb_pinyin_name | molecule | target | entrez_id`
- **Coverage**: Covers major herbs in Chinese Pharmacopoeia (中国药典)
- **Limitations**:
  - Not all herb-target interactions are experimentally validated; some are predicted by network inference
  - Compound names may differ from PubChem canonical names
  - Target gene symbols follow HGNC nomenclature
- **Grounding rule**: When reporting "X herb has Y targets", cite the exact count from `search_herb_records` output

## DisGeNET (via DOSE package)

- **Content**: 30,170 diseases, 21,671 genes, 1,134,942 disease-gene associations
- **Source**: DisGeNET v7 integrated from UNIPROT, CTD, ClinGen, CGI, GWAS Catalog, and curated literature
- **ID system**: UMLS CUI for diseases, Entrez Gene ID for genes
- **Score**: GDA score (Gene-Disease Association) ranges 0-1; higher = stronger evidence
- **Limitations**:
  - Associations include both causal and correlational evidence
  - Common diseases (e.g. diabetes, cancer) have thousands of associated genes — be specific
  - Rare diseases may have very few associations
- **Grounding rule**: Always report the exact gene count from `search_disease_targets` output

## STRING PPI Database

- **Content**: Protein-protein interactions from multiple evidence channels
- **Score**: Combined score 0-1000; recommend ≥ 400 (medium confidence) or ≥ 700 (high confidence)
- **Evidence channels**: Experimental, co-expression, text-mining, databases, co-occurrence, gene fusion, neighborhood
- **API**: REST API v11.5, species 9606 (Homo sapiens)
- **Limitations**:
  - Text-mining channel may include indirect associations
  - High threshold (700+) strongly reduces false positives but may lose real interactions
  - API rate limit: 1 request/second
- **Grounding rule**: Always report the score threshold used and the resulting node/edge counts

## Gene Ontology (GO)

- **Content**: Three ontologies — Biological Process (BP), Cellular Component (CC), Molecular Function (MF)
- **Method**: Over-representation analysis (ORA) via clusterProfiler::enrichGO()
- **Background**: org.Hs.eg.db (Homo sapiens, NCBI gene annotations)
- **Multiple testing**: Benjamini-Hochberg (p.adjust / qvalue)
- **Grounding rule**: Report p.adjust values from tool output, not raw p-values

## KEGG Pathways

- **Content**: ~350 human pathways covering metabolism, signaling, disease
- **Method**: ORA via clusterProfiler::enrichKEGG()
- **ID system**: hsa##### (e.g. hsa04151 = PI3K-Akt signaling pathway)
- **Limitations**:
  - KEGG pathway membership is manually curated and may lag behind latest research
  - Disease-specific pathways (hsa05***) may overlap with query disease — note this
- **Grounding rule**: Always use the pathway name from KEGG, not a paraphrased version

## PubMed / NCBI

- **Content**: > 36 million biomedical citations
- **API**: NCBI E-utilities (esearch + efetch)
- **Search**: Herb Pinyin name + disease English name as keywords
- **Limitations**:
  - Herb names in Pinyin may not match all publications (some use Latin botanical names)
  - Abstracts only; full text not available via API
- **Grounding rule**: Report exact publication count and date range from tool output

## PubChem

- **Content**: > 116 million compounds
- **Properties available**: Molecular weight, SMILES, InChI, XLogP, TPSA, HBond donors/acceptors
- **Limitations**:
  - TCM compound names may not resolve in PubChem; try CAS number or SMILES as fallback
  - 2D fingerprint similarity (Tanimoto) ≥ 0.8 is typically considered structurally similar

## External Cross-Validation Databases

These databases serve as independent grounding sources. When making claims about herb efficacy, compound-target relationships, or traditional uses, the agent SHOULD generate verification links and suggest cross-checking.

### HERB Database (http://47.92.70.12/)

- **Full name**: High-throughput Experiment- and Reference-guided database of TCM
- **Version**: HERB 2.0 (Nucleic Acids Res. 2025; original: Nucleic Acids Res. 2021, PMID: 33264402)
- **Content**: Integrates pharmacotranscriptomics (high-throughput gene expression) + clinical trial evidence + meta-analysis evidence for TCM
- **Coverage**: ~7,200+ herbs, ~49,000+ ingredients, ~12,000+ targets, ~28,000+ diseases
- **Key features**:
  - Experiment-based herb-target connections (gene expression profiles after herb treatment)
  - Reference-guided connections (literature-mined herb-ingredient-target-disease relationships)
  - Clinical evidence integration (HERB 2.0): clinical trials and meta-analyses for TCM formulas
- **Cross-validation use**:
  - Verify herb-target relationships found by TCMDATA against HERB's experiment-based evidence
  - Check if herb-disease associations have clinical trial support in HERB 2.0
  - Search URL pattern: `http://47.92.70.12/#/Browse?keyword={herb_name}`
- **Grounding rule**: When TCMDATA identifies herb targets, note that these are database-derived associations. Suggest the user verify key targets in HERB database, especially for experimental (gene expression) evidence

### ShenNong Alpha (https://shennongalpha.westlake.edu.cn/)

- **Full name**: ShenNong — TCM Knowledge Base by Westlake University
- **Content**: Comprehensive TCM knowledge graph integrating classical texts, pharmacological data, and modern research
- **Coverage**: Herbs (single herbs + formulas), compounds, targets, diseases, TCM syndromes (zheng), meridian tropism, properties (四气五味)
- **Key features**:
  - Classical text annotations (source citations from 本草纲目, 伤寒论, etc.)
  - TCM syndrome-disease-herb mapping (辨证论治 framework)
  - Bilingual (Chinese/English) knowledge base
  - Links traditional TCM theory (性味归经/四气五味) with modern molecular targets
- **Cross-validation use**:
  - Verify TCM traditional uses and classical indications
  - Cross-check herb properties (性味归经) and traditional disease associations
  - Validate formula composition when users ask about compound prescriptions (复方)
  - Search URL pattern: `https://shennongalpha.westlake.edu.cn/#/herb?keyword={herb_name}`
- **Grounding rule**: For traditional TCM knowledge claims (e.g. "黄芪补气固表"), cite ShenNong as a verification source. Never fabricate classical text attributions

### Cross-Validation Protocol

When reporting herb-target or herb-disease findings, follow this protocol:

1. **Primary evidence**: Always cite TCMDATA tool output with exact counts
2. **Experimental grounding**: Suggest verifying key findings in HERB database (gene expression evidence level)
3. **Clinical grounding**: For disease-treatment claims, check HERB 2.0 for clinical trial/meta-analysis support
4. **Traditional grounding**: For classical TCM claims (properties, indications, formulas), reference ShenNong
5. **Generate verification links**: Use the `generate_verification_urls` tool to produce clickable links for the user

Evidence hierarchy for herb-disease claims:
- **Strongest**: Clinical trial evidence (HERB 2.0) + PubMed literature
- **Strong**: Experimental evidence (HERB gene expression) + computational prediction (TCMDATA)
- **Moderate**: Database association (TCMDATA) + text-mining (PubMed)
- **Weakest**: Single database association alone — always flag as "requires further validation"

## Important Anti-Hallucination Rules

1. **Never invent gene-disease associations** not present in DisGeNET output
2. **Never fabricate enrichment terms** — only report terms from GO/KEGG tool output with their exact p.adjust values
3. **Never claim a herb "is known to treat" a disease** without PubMed evidence from the tool
4. **Never make up PPI interaction counts** — always from STRING tool output
5. **If a tool returns an error or empty result**, report it honestly. Absence of evidence is NOT evidence of absence
6. **Distinguish between**: "X gene is a hub gene in the PPI network" (topology metric) vs "X gene is a validated therapeutic target" (experimental evidence). The former is computational, the latter requires literature support
7. **Never fabricate classical TCM attributions** — do not claim "according to 本草纲目..." unless verified via ShenNong or explicit user input
8. **Never claim clinical efficacy** without evidence — distinguish "computational prediction" from "clinically validated treatment"
9. **Always provide verification links** — when reporting key herb-target or herb-disease findings, generate HERB/ShenNong URLs so users can cross-check
