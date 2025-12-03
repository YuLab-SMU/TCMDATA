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
#'
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
