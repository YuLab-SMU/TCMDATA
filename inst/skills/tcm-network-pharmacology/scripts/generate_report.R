#!/usr/bin/env Rscript
# generate_report.R
# Collects all artifacts from the current session and generates a structured
# network pharmacology analysis report in Markdown format.
#
# Usage (by agent via execute_skill_script):
#   Rscript generate_report.R [output_file]
#
# Expects TCMDATA to be loaded and artifacts to be available in the session.

args <- commandArgs(trailingOnly = TRUE)
output_file <- if (length(args) >= 1) args[1] else "network_pharmacology_report.md"

if (!requireNamespace("TCMDATA", quietly = TRUE)) {
  stop("TCMDATA package is required.")
}

artifacts <- TCMDATA::list_tcm_artifacts()
if (is.null(artifacts) || nrow(artifacts) == 0) {
  stop("No artifacts found in current session. Run the analysis first.")
}

# Helper: find artifact by type or function name
find_artifact <- function(type = NULL, func = NULL) {
  if (!is.null(type)) {
    idx <- which(artifacts$artifact_type == type)
  } else if (!is.null(func)) {
    idx <- which(artifacts$function_name == func)
  } else {
    return(NULL)
  }
  if (length(idx) == 0) return(NULL)
  # Return the latest one
  TCMDATA::list_tcm_artifacts()$artifact_id[idx[length(idx)]]
}

# Collect sections
sections <- list()

sections$header <- paste0(
  "# Network Pharmacology Analysis Report\n\n",
  "Generated: ", Sys.time(), "\n\n",
  "---\n"
)

# Herb search
herb_id <- find_artifact(func = "search_herb")
if (!is.null(herb_id)) {
  herb_art <- get(herb_id, envir = TCMDATA:::.artifact_env)
  sections$herb <- sprintf(
    "## 1. Herb Target Retrieval\n\n%s\n\n- Targets: %d unique genes\n",
    herb_art$summary %||% "",
    length(unique(herb_art$object$target))
  )
}

# Disease search
disease_id <- find_artifact(func = "search_disease")
if (!is.null(disease_id)) {
  disease_art <- get(disease_id, envir = TCMDATA:::.artifact_env)
  sections$disease <- sprintf(
    "## 2. Disease Target Retrieval\n\n%s\n\n- Genes: %d\n",
    disease_art$summary %||% "",
    length(unique(disease_art$object$gene_id))
  )
}

# Intersection
intersect_id <- find_artifact(type = "intersection_result")
if (!is.null(intersect_id)) {
  int_art <- get(intersect_id, envir = TCMDATA:::.artifact_env)
  sections$intersection <- sprintf(
    "## 3. Target Intersection\n\n%s\n\n- Common targets: %d\n",
    int_art$summary %||% "",
    length(int_art$object$intersection)
  )
}

# PPI
ppi_id <- find_artifact(func = "rank_ppi_nodes")
if (is.null(ppi_id)) ppi_id <- find_artifact(type = "ppi_network")
if (!is.null(ppi_id)) {
  ppi_art <- get(ppi_id, envir = TCMDATA:::.artifact_env)
  sections$ppi <- sprintf(
    "## 4. PPI Network & Hub Genes\n\n%s\n",
    ppi_art$summary %||% ""
  )
}

# Enrichment
go_id <- find_artifact(func = "run_go_enrichment")
kegg_id <- find_artifact(func = "run_kegg_enrichment")
enrich_text <- "## 5. Functional Enrichment\n\n"
if (!is.null(go_id)) {
  go_art <- get(go_id, envir = TCMDATA:::.artifact_env)
  enrich_text <- paste0(enrich_text, "### GO-BP\n\n", go_art$summary %||% "", "\n\n")
}
if (!is.null(kegg_id)) {
  kegg_art <- get(kegg_id, envir = TCMDATA:::.artifact_env)
  enrich_text <- paste0(enrich_text, "### KEGG\n\n", kegg_art$summary %||% "", "\n\n")
}
sections$enrichment <- enrich_text

# PubMed
pubmed_id <- find_artifact(func = "get_pubmed_data")
if (!is.null(pubmed_id)) {
  pub_art <- get(pubmed_id, envir = TCMDATA:::.artifact_env)
  sections$pubmed <- sprintf(
    "## 6. Literature Validation\n\n%s\n",
    pub_art$summary %||% ""
  )
}

# Write report
report <- paste(unlist(sections), collapse = "\n\n")
writeLines(report, output_file)
cat(sprintf("Report written to: %s (%d bytes)\n", output_file, nchar(report)))
