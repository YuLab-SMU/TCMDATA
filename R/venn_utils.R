#' Construct a Logical Matrix for Venn Diagram Visualization
#'
#' This function takes 2 to 4 named character vectors as input, removes duplicates,
#' and returns a data frame in logical format that is compatible with `ggvenn` for Venn diagram plotting.
#'
#' @param ... 2 to 4 character vectors representing different gene sets or element sets.
#' @param set_names Optional character vector to assign custom names to the sets.
#'        If not provided, names will be automatically assigned as "Set1", "Set2", etc.
#'
#' @return A \code{data.frame} with:
#' \itemize{
#'   \item \code{Element}: all unique elements across the sets.
#'   \item One logical column per input set indicating presence (TRUE/FALSE) of each element.
#' }
#' @examples
#' # Usage with 3 sets
#' targets_db1 <- c("TP53", "EGFR", "KRAS", "MYC", "AKT1")
#' targets_db2 <- c("TP53", "KRAS", "PTEN", "BRCA1")
#' targets_db3 <- c("EGFR", "MYC", "BRAF", "PTEN", "TP53")
#'
#' df_venn3 <- getvenndata(targets_db1, targets_db2, targets_db3,
#'                         set_names = c("Database_1", "Database_2", "Database_3"))
#'
#' head(df_venn3)
#' @export
getvenndata <- function(..., set_names = NULL) {

  input_sets <- list(...)

  # check the number pf input vectors
  if (length(input_sets) < 2 || length(input_sets) > 4) {
    stop("please put 2-4 vectors instead.")
  }

  # set names
  if (is.null(set_names)) {
    if (is.null(names(input_sets)) || any(names(input_sets) == "")) {
      names(input_sets) <- paste0("Set", seq_along(input_sets))
    }
  } else {
    if (length(set_names) != length(input_sets)) {
      stop("the length of 'set_names' should be the same as the input sets")
    }
    names(input_sets) <- set_names
  }

  input_sets <- lapply(input_sets, function(x) unique(as.character(x)))
  all_elements <- unique(unlist(input_sets))

  # construct logical df for ggvenn
  venn_df <- data.frame(
    Element = all_elements,
    stringsAsFactors = FALSE
  )

  for (name in names(input_sets)) {
    venn_df[[name]] <- all_elements %in% input_sets[[name]]
  }

  return(venn_df)
}


#' Get all Venn set intersections from getvenndata-style logical matrix
#'
#' @param venn_df A data.frame from getvenndata(), with logical columns
#' @param col_names Character vector of columns to use (default: all except 1st)
#' @param drop_empty Logical; whether to drop combinations with 0 genes
#'
#' @importFrom dplyr filter
#' @examples
#' # Usage with 3 sets
#' targets_db1 <- c("TP53", "EGFR", "KRAS", "MYC", "AKT1")
#' targets_db2 <- c("TP53", "KRAS", "PTEN", "BRCA1")
#' targets_db3 <- c("EGFR", "MYC", "BRAF", "PTEN", "TP53")
#'
#' df_venn3 <- getvenndata(targets_db1, targets_db2, targets_db3,
#'                         set_names = c("Database_1", "Database_2", "Database_3"))
#'
#' venn_res <- getvennresult(df_venn3)
#' print(venn_res)
#' @return A data.frame with intersection combinations and their gene lists
#' @export
getvennresult <- function(venn_df,
                          col_names = NULL,
                          drop_empty = TRUE) {

  if (is.null(col_names)) {
    col_names <- colnames(venn_df)[-1]
  }

  if (!all(sapply(venn_df[col_names], is.logical))) {
    stop("All selected columns must be logical (TRUE/FALSE).")
  }

  gene_col <- venn_df[[1]]

  combos <- expand.grid(rep(list(c(TRUE, FALSE)), length(col_names)))
  colnames(combos) <- col_names

  results <- apply(combos, 1, function(mask) {
    filter <- Reduce(`&`, Map(function(col, val) venn_df[[col]] == val, col_names, mask))
    matched_genes <- gene_col[filter]
    set_combo <- paste(col_names[which(mask)], collapse = "&")
    if (set_combo == "") set_combo <- "None"
    data.frame(
      Set_Combination = set_combo,
      Gene_Count = length(matched_genes),
      Genes = paste(matched_genes, collapse = ", ")
    )
  })

  df <- do.call(rbind, results)
  rownames(df) <- NULL

  if (drop_empty) {
    df <- dplyr::filter(df, .data$Gene_Count > 0)
  }

  return(df)
}

