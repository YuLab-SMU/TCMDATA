# geo_search.R
# Search NCBI GEO for relevant datasets using E-utilities.

#' Search GEO for datasets related to a disease or keyword
#'
#' Queries the NCBI GEO database (via E-utilities) to find relevant
#' expression datasets (GDS/GSE) for a given search term. Useful for
#' discovering RNA-seq, microarray, or single-cell datasets before analysis.
#'
#' @param query Character. Search query (e.g. disease name, tissue + disease).
#' @param organism Character. Organism filter. Default \code{"Homo sapiens"}.
#' @param dataset_type Character. Optional type filter: \code{"expression profiling by array"},
#'   \code{"expression profiling by high throughput sequencing"}, or \code{"single-cell"}.
#'   Pass \code{NULL} to search all types.
#' @param email Character. Email address required by NCBI E-utilities.
#' @param retmax Integer. Maximum number of results. Default 20.
#'
#' @return A data.frame with columns: gse, title, summary, platform, n_samples, type, organism.
#'   Returns \code{NULL} if no results found.
#'
#' @examples
#' \dontrun{
#'   res <- search_geo_datasets("sepsis", email = "user@example.com")
#'   head(res)
#' }
#' @export
search_geo_datasets <- function(query,
                                organism = "Homo sapiens",
                                dataset_type = NULL,
                                email = NULL,
                                retmax = 20L) {
  if (!requireNamespace("httr", quietly = TRUE)) {
    stop("Package 'httr' is required. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("xml2", quietly = TRUE)) {
    stop("Package 'xml2' is required. Please install it.", call. = FALSE)
  }
  if (is.null(email) || !nzchar(email)) {
    stop("NCBI requires an email address to use E-utils. Please provide one.", call. = FALSE)
  }

  # Build query
  search_term <- paste0(query, "[Title]")
  if (!is.null(organism) && nzchar(organism)) {
    search_term <- paste0(search_term, " AND ", organism, "[Organism]")
  }

  # Map dataset_type to GEO type terms
  if (!is.null(dataset_type) && nzchar(dataset_type)) {
    type_map <- list(
      "single-cell" = "expression profiling by high throughput sequencing",
      "rnaseq"      = "expression profiling by high throughput sequencing",
      "microarray"  = "expression profiling by array"
    )
    mapped <- type_map[[tolower(dataset_type)]]
    if (!is.null(mapped)) {
      search_term <- paste0(search_term, " AND ", mapped, "[DataSet Type]")
    } else {
      search_term <- paste0(search_term, " AND ", dataset_type, "[DataSet Type]")
    }
  }

  # If single-cell, also add scRNA keywords to title/summary
  if (!is.null(dataset_type) && tolower(dataset_type) == "single-cell") {
    search_term <- paste0(
      "(", query, " AND (single-cell OR scRNA-seq OR single cell RNA))[Title] AND ",
      organism, "[Organism]"
    )
  }

  message("GEO search query: ", search_term)

  # Step 1: ESearch to get GDS IDs
  search_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
  search_res <- httr::GET(search_url, query = list(
    db      = "gds",
    term    = search_term,
    retmax  = as.integer(retmax),
    retmode = "xml",
    email   = email,
    sort    = "relevance"
  ))
  httr::stop_for_status(search_res)

  search_xml <- xml2::read_xml(httr::content(search_res, "text", encoding = "UTF-8"))
  ids <- xml2::xml_text(xml2::xml_find_all(search_xml, "//IdList/Id"))
  total_count <- as.numeric(xml2::xml_text(xml2::xml_find_first(search_xml, "//Count")))

  message("GEO total matches: ", total_count, "; retrieving top ", min(length(ids), retmax))


  if (length(ids) == 0L) {
    message("No GEO datasets found for: ", query)
    return(NULL)
  }

  # Step 2: ESummary to get details
  summary_url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
  id_string <- paste(ids, collapse = ",")
  summary_res <- httr::GET(summary_url, query = list(
    db      = "gds",
    id      = id_string,
    retmode = "xml",
    email   = email
  ))
  httr::stop_for_status(summary_res)

  summary_xml <- xml2::read_xml(httr::content(summary_res, "text", encoding = "UTF-8"))
  doc_sums <- xml2::xml_find_all(summary_xml, "//DocSum")

  if (length(doc_sums) == 0L) {
    message("No GEO dataset summaries returned.")
    return(NULL)
  }

  # Parse each DocSum
  records <- lapply(doc_sums, function(ds) {
    .extract_item <- function(node, item_name) {
      item <- xml2::xml_find_first(node, paste0(".//Item[@Name='", item_name, "']"))
      if (is.na(item)) return(NA_character_)
      xml2::xml_text(item)
    }

    accession  <- .extract_item(ds, "Accession")
    title      <- .extract_item(ds, "title")
    summary    <- .extract_item(ds, "summary")
    taxon      <- .extract_item(ds, "taxon")
    gds_type   <- .extract_item(ds, "gdsType")
    platform   <- .extract_item(ds, "GPL")
    n_samples  <- .extract_item(ds, "n_samples")

    # Truncate summary if too long
    if (!is.na(summary) && nchar(summary) > 300) {
      summary <- paste0(substr(summary, 1, 297), "...")
    }

    data.frame(
      gse       = accession,
      title     = title,
      summary   = summary,
      platform  = ifelse(is.na(platform), NA_character_, paste0("GPL", platform)),
      n_samples = as.integer(ifelse(is.na(n_samples), NA_integer_, n_samples)),
      type      = gds_type,
      organism  = taxon,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, records)

  # Filter to only GSE entries (exclude GDS if present)
  if (nrow(result) > 0L) {
    result <- result[grepl("^GSE", result$gse, ignore.case = TRUE) |
                     grepl("^GDS", result$gse, ignore.case = TRUE), , drop = FALSE]
  }

  if (nrow(result) == 0L) {
    message("No valid GSE/GDS records after filtering.")
    return(NULL)
  }

  rownames(result) <- NULL
  result
}
