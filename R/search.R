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
search_herb <- function(herb, type) {
  type <- match.arg(type, c("Herb_cn_name", "Herb_pinyin_name", "Herb_en_name"))

  available <- unique(tcm_data[[type]])
  herb_not_exist <- setdiff(herb, available)
  if (length(herb_not_exist) > 0L) {
    warning(paste0(
      paste(herb_not_exist, collapse = ", "),
      " doesn't/don't exist in our dataset."
    ))
  }
  herb <- intersect(herb, available)

  result <- tcm_data |>
    dplyr::filter(.data[[type]] %in% herb) |>
    dplyr::select("Herb_pinyin_name", "molecule", "target") |>
    dplyr::distinct(.keep_all = TRUE) |>
    tidyr::drop_na()
  names(result) <- c("herb", "molecule", "target")
  rownames(result) <- NULL
  result
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
  herbs_data <- tcm_data |>
    dplyr::filter(.data$target %in% gene_list) |>
    dplyr::select(c("Herb_pinyin_name", "molecule", "target")) |>
    dplyr::distinct(.keep_all = TRUE) |>
    tidyr::drop_na()
  
  colnames(herbs_data) <- c("herb", "molecule", "target")
  rownames(herbs_data) <- NULL
  return(herbs_data)
}
