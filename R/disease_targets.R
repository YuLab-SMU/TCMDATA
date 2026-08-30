# Open Targets Platform GraphQL API utilities
# API endpoint: https://api.platform.opentargets.org/api/v4/graphql

.OT_GRAPHQL_URL <- "https://api.platform.opentargets.org/api/v4/graphql"

#' @keywords internal
#' @noRd
.check_httr <- function() {
  if (!requireNamespace("httr", quietly = TRUE))
    stop("Package 'httr' is required. Install with: install.packages('httr')", call. = FALSE)
}

#' @keywords internal
#' @noRd
.ot_post <- function(query, variables = list()) {
  .check_httr()
  resp <- httr::POST(
    .OT_GRAPHQL_URL,
    httr::content_type_json(),
    body = jsonlite::toJSON(list(query = query, variables = variables), auto_unbox = TRUE)
  )
  httr::stop_for_status(resp)
  httr::content(resp, as = "parsed", simplifyVector = TRUE)
}


#' Search Disease EFO IDs from Open Targets Platform
#'
#' Queries the Open Targets Platform GraphQL API to search for diseases
#' matching a given name and returns their EFO/MONDO/Orphanet IDs.
#'
#' @param disease_name A character string. Disease name to search for
#'   (English, e.g. \code{"breast cancer"}).
#' @param size An integer. Maximum number of results to return (default \code{10}).
#'
#' @return A \code{data.frame} with columns \code{id} (EFO/MONDO/Orphanet ID)
#'   and \code{name} (disease name). Returns an empty \code{data.frame} if no
#'   match is found.
#'
#' @examples
#' \dontrun{
#' search_disease_efo("breast cancer")
#' search_disease_efo("diabetes", size = 5)
#' }
#'
#' @export
#'
search_disease_efo <- function(disease_name, size = 10) {
  query <- '
  query SearchDisease($query: String!, $size: Int!) {
    search(queryString: $query, entityNames: ["disease"], page: {index: 0, size: $size}) {
      hits {
        id
        name
        entity
      }
    }
  }'

  parsed <- .ot_post(query, list(query = disease_name, size = as.integer(size)))
  hits <- parsed$data$search$hits

  if (is.null(hits) || length(hits) == 0)
    return(data.frame(id = character(), name = character(), stringsAsFactors = FALSE))

  result <- hits[hits$entity == "disease", c("id", "name"), drop = FALSE]
  rownames(result) <- NULL
  result
}


#' Retrieve Disease-Associated Targets from Open Targets Platform
#'
#' Queries the Open Targets Platform GraphQL API for targets associated with
#' a disease, identified by its EFO/MONDO/Orphanet ID.
#'
#' @param efo_id A character string. Disease ontology ID (e.g.
#'   \code{"EFO_0000305"} for breast carcinoma, \code{"EFO_0000400"} for
#'   diabetes mellitus).
#' @param size An integer. Maximum number of targets to return (default \code{200}).
#' @param score_threshold A numeric value between 0 and 1. Only targets with an
#'   overall association score at or above this threshold are returned
#'   (default \code{0}).
#'
#' @return A \code{data.frame} with columns:
#'   \describe{
#'     \item{ensembl_id}{Ensembl gene ID.}
#'     \item{gene_symbol}{HGNC-approved gene symbol.}
#'     \item{gene_name}{Full gene name.}
#'     \item{biotype}{Gene biotype (e.g. \code{"protein_coding"}).}
#'     \item{score}{Overall association score (0-1) from Open Targets.}
#'   }
#'   Rows are sorted by \code{score} in descending order.
#'
#' @examples
#' \dontrun{
#' # Breast carcinoma targets with score >= 0.5
#' get_disease_targets("EFO_0000305", size = 100, score_threshold = 0.5)
#'
#' # All diabetes mellitus targets
#' get_disease_targets("EFO_0000400", size = 500)
#' }
#'
#' @export
#'
get_disease_targets <- function(efo_id, size = 200, score_threshold = 0) {
  query <- '
  query DiseaseTargets($efoId: String!, $size: Int!) {
    disease(efoId: $efoId) {
      id
      name
      associatedTargets(page: {index: 0, size: $size}) {
        count
        rows {
          target {
            id
            approvedSymbol
            approvedName
            biotype
          }
          score
        }
      }
    }
  }'

  parsed <- .ot_post(query, list(efoId = efo_id, size = as.integer(size)))
  disease_data <- parsed$data$disease

  if (is.null(disease_data))
    stop("No disease found for EFO ID: ", efo_id, call. = FALSE)

  rows <- disease_data$associatedTargets$rows
  if (is.null(rows) || nrow(rows) == 0)
    return(data.frame(
      ensembl_id = character(), gene_symbol = character(),
      gene_name  = character(), biotype = character(),
      score      = numeric(),
      stringsAsFactors = FALSE
    ))

  result <- data.frame(
    ensembl_id  = rows$target$id,
    gene_symbol = rows$target$approvedSymbol,
    gene_name   = rows$target$approvedName,
    biotype     = rows$target$biotype,
    score       = rows$score,
    stringsAsFactors = FALSE
  )

  result <- result[result$score >= score_threshold, , drop = FALSE]
  result <- result[order(-result$score), , drop = FALSE]
  rownames(result) <- NULL
  result
}


#' Query Disease Targets by Disease Name from Open Targets Platform
#'
#' A convenience wrapper that combines disease name search and target retrieval.
#' Searches for a disease by name, selects a match, and returns its associated
#' targets with overall association scores.
#'
#' @param disease_name A character string. Disease name to query
#'   (English, e.g. \code{"breast cancer"}, \code{"diabetes"}).
#' @param size An integer. Maximum number of targets to return (default \code{200}).
#' @param score_threshold A numeric value between 0 and 1. Minimum association
#'   score for included targets (default \code{0}).
#' @param efo_index An integer. When the search returns multiple disease matches,
#'   specifies which result to use (default \code{1}, i.e. the top hit). Use
#'   \code{search_disease_efo()} to inspect all candidates first.
#'
#' @return A \code{data.frame} with columns \code{ensembl_id}, \code{gene_symbol},
#'   \code{gene_name}, \code{biotype}, and \code{score}, sorted by \code{score}
#'   descending. Returns \code{NULL} if the disease name yields no results.
#'
#' @seealso \code{\link{search_disease_efo}}, \code{\link{get_disease_targets}}
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' targets <- query_disease_targets("breast cancer")
#' head(targets)
#'
#' # Stricter filter, more results
#' targets <- query_disease_targets("diabetes", size = 500, score_threshold = 0.1)
#'
#' # Inspect candidates, then select the correct one
#' search_disease_efo("diabetes")
#' targets <- query_disease_targets("diabetes", efo_index = 2)
#'
#' # Extract gene symbol vector for downstream analysis
#' gene_set <- targets$gene_symbol
#' }
#'
#' @export
#'
query_disease_targets <- function(disease_name, size = 200, score_threshold = 0, efo_index = 1) {
  diseases <- search_disease_efo(disease_name)

  if (nrow(diseases) == 0) {
    message("No matching disease found for: ", disease_name)
    return(NULL)
  }

  if (efo_index > nrow(diseases))
    stop("efo_index (", efo_index, ") exceeds number of matches found (",
         nrow(diseases), ").", call. = FALSE)

  message(sprintf("Using disease: %s (%s)  [match %d of %d]",
    diseases$name[efo_index], diseases$id[efo_index],
    efo_index, nrow(diseases)))

  get_disease_targets(diseases$id[efo_index], size = size, score_threshold = score_threshold)
}


# DisGeNET access was removed after its licensing model changed. Keep the
# exported functions below as compatibility stubs so existing scripts fail with
# an actionable message instead of a missing-function error.
.disgenet_unavailable <- function() {
  stop(
    paste(
      "DisGeNET is now a commercial resource.",
      "The DisGeNET-derived data are unavailable in TCMDATA due to licensing and copyright restrictions.",
      "Use query_disease_targets() to retrieve disease-associated targets from the Open Targets Platform."
    ),
    call. = FALSE
  )
}

#' DisGeNET disease-to-gene lookup (unavailable)
#'
#' This compatibility interface is retained for existing scripts, but its
#' underlying DisGeNET-derived data are no longer available because of
#' licensing and copyright restrictions. Use \code{query_disease_targets()}
#' for disease-associated targets from the Open Targets Platform.
#'
#' @param disease Character. Disease name or UMLS CUI (e.g. "sepsis" or
#'   "C0243026"). Supports a vector for multiple diseases.
#' @param readable Logical. Convert Entrez IDs to gene symbols (default TRUE).
#'
#' @return This function always raises an error explaining that the data are
#'   unavailable.
#' @seealso \code{\link{query_disease_targets}}
#' @export
search_disease <- function(disease, readable = TRUE) {
  .disgenet_unavailable()
}


#' DisGeNET gene-to-disease lookup (unavailable)
#'
#' This compatibility interface is retained for existing scripts, but its
#' underlying DisGeNET-derived data are no longer available because of
#' licensing and copyright restrictions. There is currently no reverse
#' gene-to-disease replacement in TCMDATA.
#'
#' @param gene Character. Gene symbols (e.g. "TNF") or Entrez IDs.
#'   Supports a vector for multiple genes.
#' @param readable Logical. Attach gene symbol column (default TRUE).
#'
#' @return This function always raises an error explaining that the data are
#'   unavailable.
#' @seealso \code{\link{query_disease_targets}}
#' @export
search_gene_disease <- function(gene, readable = TRUE) {
  .disgenet_unavailable()
}
