#' Retrieve Herbs and Retrieve Associated Molecules and Targets
#'
#' This function retrieves herb, molecule, and target information from the internal dataset, \code{tcm_data}, 
#'    based on specified herb names and their corresponding name types (Chinese, Pinyin, or English).
#'
#' @param herb A character vector containing the names of herbs to be queried.
#' @param type A string indicating the type of herb names provided. Must be one of \code{"Herb_cn_name"}, 
#'    \code{"Herb_pinyin_name"}, or \code{"Herb_en_name"}.
#'
#' @return A \code{data.frame} with three columns: \code{herb}, \code{molecule}, 
#'    and \code{target}, containing the corresponding information for the specified herbs. 
#'    Rows with \code{NA} values are excluded.
#'
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr distinct
#' @importFrom dplyr %>%
#' @importFrom tidyr drop_na
#' @importFrom rlang .data
#' @examples
#' \dontrun{
#' # Example usage with Chinese herb names
#' herbs <- c("灵芝")
#' lz <- search_herb(herb = herbs, type = "Herb_cn_name")
#' print(lz)
#'
#' # Example usage with pinyin herb names
#' herbs <- c("Ginseng", "Licorice")
#' comp <- search_herb(herb = herbs, type = "Herb_pinyin_name")
#' print(comp)
#' }
#'
#' @export
#' 
search_herb <- function(herb, type){
  # Validate the 'type' parameter
  type <- match.arg(type, c("Herb_cn_name", "Herb_pinyin_name", "Herb_en_name"))
  
  # Check existence of herbs based on the specified type
  if (type == "Herb_cn_name"){
    if (all(herb %in% unique(tcm_data$Herb_cn_name)) == FALSE){
      herb_not_exist <- setdiff(herb, unique(tcm_data$Herb_cn_name))
      warning(paste0(paste(herb_not_exist, collapse=", "), " doesn't/don't exist in our dataset."))
      herb <- herb[-match(herb_not_exist, herb)]
    }
    
    result <- tcm_data %>%
      dplyr::filter(.data$Herb_cn_name %in% herb) %>%
      dplyr::select(c("Herb_pinyin_name", "molecule", "target")) %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      tidyr::drop_na()
    colnames(result) <- c("herb", "molecule", "target")
    rownames(result) <- NULL
    return(result)
  }
  
  if (type == "Herb_pinyin_name"){
    if (all(herb %in% unique(tcm_data$Herb_pinyin_name)) == FALSE){
      herb_not_exist <- setdiff(herb, unique(tcm_data$Herb_pinyin_name))
      warning(paste0(paste(herb_not_exist, collapse=", "), " doesn't/don't exist in our dataset."))
      herb <- herb[-match(herb_not_exist, herb)]
    }
    
    result <- tcm_data %>%
      dplyr::filter(.data$Herb_pinyin_name %in% herb) %>%
      dplyr::select(c("Herb_pinyin_name", "molecule", "target")) %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      tidyr::drop_na()
    colnames(result) <- c("herb", "molecule", "target")
    rownames(result) <- NULL
    return(result)
  }
  
  #if (type == "Herb_en_name"){
  if (all(herb %in% unique(tcm_data$Herb_en_name)) == FALSE){
    herb_not_exist <- setdiff(herb, unique(tcm_data$Herb_en_name))
    warning(paste0(paste(herb_not_exist, collapse=", "), " doesn't/don't exist in our dataset."))
    herb <- herb[-match(herb_not_exist, herb)]
  }
    
  result <- tcm_data %>%
    dplyr::filter(.data$Herb_en_name %in% herb) %>%
    dplyr::select("Herb_pinyin_name", "molecule", "target") %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    tidyr::drop_na()
  colnames(result) <- c("herb", "molecule", "target")
  rownames(result) <- NULL

  return(result)
  #}
}


#' Retrieve Herbs, Molecules, and Targets Based on Gene List
#'
#' This function retrieves herb, molecule, and target information from the internal dataset, \code{tcm_data},
#'    based on a provided list of target genes.
#'
#' @param gene_list A character vector containing gene symbols which can be determined freely by users.
#'
#' @return A \code{data.frame} with three columns: \code{herb}, \code{molecule}, and \code{target}, 
#'    containing information related to the specified genes. Rows with \code{NA} values are excluded.
#'
#' @examples
#' \dontrun{
#' # Example usage with a list of gene symbols
#' genes <- c("TP53", "EGFR", "BRCA1")
#' herbs_targets <- search_target(genes)
#' print(herbs_targets)
#' }
#'
#' @export
#' 
search_target <- function(gene_list){  
  # Check existence of genes in the dataset
  if (!all(gene_list %in% unique(tcm_data$target))){
    gene_diff <- setdiff(gene_list, unique(tcm_data$target))
    warning(paste0(paste(gene_diff, collapse=", "), " doesn't/don't exist in the datasets."))
    gene_list <- gene_list[-match(gene_diff, gene_list)]
  }
  
  # Retrieve relevant data
  herbs_data <- tcm_data %>%
    dplyr::filter(.data$target %in% gene_list) %>%
    dplyr::select(c("Herb_pinyin_name", "molecule", "target")) %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    tidyr::drop_na()
  
  colnames(herbs_data) <- c("herb", "molecule", "target")
  rownames(herbs_data) <- NULL
  return(herbs_data)
}



# Disease-target search (DisGeNET via DOSE)
# Package-level cache
.disease_env <- new.env(parent = emptyenv())

#' @keywords internal
#' @noRd
.check_dose <- function() {
  if (!requireNamespace("DOSE", quietly = TRUE))
    stop("Package 'DOSE' is required. Install: BiocManager::install('DOSE')", call. = FALSE)
}

#' @keywords internal
#' @noRd
.get_disease_data <- function() {
  if (!is.null(.disease_env$df)) return(.disease_env$df)
  .check_dose()

  env <- new.env(parent = emptyenv())
  utils::data("DGN_PATHID2EXTID", package = "DOSE", envir = env)
  utils::data("DGN_PATHID2NAME", package = "DOSE", envir = env)

  id2gene <- env$DGN_PATHID2EXTID
  id2name <- env$DGN_PATHID2NAME

  df <- data.frame(
    disease_id = rep(names(id2gene), lengths(id2gene)),
    gene_id = as.character(unlist(id2gene, use.names = FALSE)),
    stringsAsFactors = FALSE
  )
  df$disease_name <- id2name[df$disease_id]

  .disease_env$df <- df
  .disease_env$id2name <- id2name
  df
}

#' @keywords internal
#' @noRd
.add_gene_symbol <- function(df) {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("Package 'org.Hs.eg.db' required. Install: BiocManager::install('org.Hs.eg.db')", call. = FALSE)
  mapping <- suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                          keys = unique(df$gene_id),
                          columns = "SYMBOL", keytype = "ENTREZID")
  )
  sym_map <- stats::setNames(mapping$SYMBOL, mapping$ENTREZID)
  df$symbol <- sym_map[df$gene_id]
  df
}

#' @keywords internal
#' @noRd
.symbol_to_entrez <- function(symbols) {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("Package 'org.Hs.eg.db' required. Install: BiocManager::install('org.Hs.eg.db')", call. = FALSE)
  suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                          keys = symbols, columns = "ENTREZID", keytype = "SYMBOL")
  )
}

#' Search disease targets (disease -> genes)
#'
#' Query DisGeNET (via DOSE) to find genes associated with a disease.
#' Supports UMLS CUI IDs, exact name match, or fuzzy name search.
#'
#' @param disease Character. Disease name or UMLS CUI (e.g. "sepsis" or
#'   "C0243026"). Supports a vector for multiple diseases.
#' @param readable Logical. Convert Entrez IDs to gene symbols (default TRUE).
#'
#' @return A \code{data.frame} with columns \code{disease_id}, \code{disease_name},
#'   \code{gene_id}, and optionally \code{symbol}. Returns NULL if no match.
#'
#' @examples
#' \dontrun{
#'   search_disease("sepsis")
#'   search_disease(c("sepsis", "asthma"))
#'   search_disease("C0243026")
#' }
#' @export
search_disease <- function(disease, readable = TRUE) {
  df <- .get_disease_data()
  id2name <- .disease_env$id2name

  # Match: exact ID -> exact name (case-insensitive) -> grep
  all_matched_ids <- character(0)
  for (q in disease) {
    if (q %in% names(id2name)) {
      all_matched_ids <- c(all_matched_ids, q)
    } else {
      exact <- names(id2name)[tolower(id2name) == tolower(q)]
      if (length(exact) > 0) {
        all_matched_ids <- c(all_matched_ids, exact)
      } else {
        fuzzy <- names(id2name)[grepl(q, id2name, ignore.case = TRUE)]
        if (length(fuzzy) > 0) {
          all_matched_ids <- c(all_matched_ids, fuzzy)
        } else {
          message("No disease found for: ", q)
        }
      }
    }
  }

  if (length(all_matched_ids) == 0) return(NULL)

  res <- df[df$disease_id %in% all_matched_ids, , drop = FALSE]
  if (readable) res <- .add_gene_symbol(res)
  res <- res[order(res$disease_name, res$gene_id), , drop = FALSE]
  rownames(res) <- NULL
  res
}


#' Search gene-associated diseases (gene -> diseases)
#'
#' Reverse lookup: given gene symbols or Entrez IDs, find associated diseases
#' from DisGeNET (via DOSE).
#'
#' @param gene Character. Gene symbols (e.g. "TNF") or Entrez IDs.
#'   Supports a vector for multiple genes.
#' @param readable Logical. Attach gene symbol column (default TRUE).
#'
#' @return A \code{data.frame} with columns \code{disease_id}, \code{disease_name},
#'   \code{gene_id}, and optionally \code{symbol}. Returns NULL if no match.
#'
#' @examples
#' \dontrun{
#'   search_gene_disease("TNF")
#'   search_gene_disease(c("IL6", "TNF", "PPARG"))
#'   search_gene_disease("7124")
#' }
#' @export
search_gene_disease <- function(gene, readable = TRUE) {
  df <- .get_disease_data()

  # Try as Entrez ID first
  hits <- df[df$gene_id %in% gene, , drop = FALSE]

  # If no hits, try Symbol -> Entrez conversion
  if (nrow(hits) == 0) {
    gene_map <- .symbol_to_entrez(gene)
    entrez_ids <- gene_map$ENTREZID[!is.na(gene_map$ENTREZID)]
    if (length(entrez_ids) > 0) {
      hits <- df[df$gene_id %in% entrez_ids, , drop = FALSE]
    }
  }

  if (nrow(hits) == 0) {
    message("No diseases found for: ", paste(gene, collapse = ", "))
    return(NULL)
  }

  if (readable) hits <- .add_gene_symbol(hits)
  hits <- hits[order(hits$gene_id, hits$disease_name), , drop = FALSE]
  rownames(hits) <- NULL
  hits
}

