#' Get all Venn set intersections from getvenndata-style logical matrix
#'
#' @param venn_df A data.frame from getvenndata(), with logical columns
#' @param col_names Character vector of columns to use (default: all except 1st)
#' @param drop_empty Logical; whether to drop combinations with 0 genes
#'
#' @importFrom dplyr filter
#' @return A data.frame with intersection combinations and their gene lists
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
