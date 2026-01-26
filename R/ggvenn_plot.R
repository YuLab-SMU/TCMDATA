#' Plot a Venn Diagram for 2–4 Sets Using ggvenn
#'
#' This is a wrapper around \code{ggvenn()} to draw customizable Venn diagrams 
#' for 2–4 sets with automatic color matching and error handling.
#'
#' @param venn_df A data.frame where the first column is element ID and the remaining 2–4 columns are logical (TRUE/FALSE) indicators for set membership.
#' @param col_names Character vector of column names to be used as sets. If NULL, all columns except the first are used.
#' @param set.color Fill colors for Venn sets.
#' @param set.name.color Color for set name labels (ignored if \code{use.color.as.text = TRUE}).
#' @param use.color.as.text Logical; if TRUE, use fill color as set name text color.
#' @param name.size Font size for set names.
#' @param text.size Font size for intersection counts.
#' @param stroke.color Circle border color. Default is 'black'.
#' @param stroke.size Circle border thickness. Default is 0.6.
#' @param show.percentage Logical; whether to show percentages instead of raw counts.
#' @param show.elements Logical; whether to display individual elements in each region.
#' @param digits Number of decimal places for percentage display.
#' @param expand_ratio Plot margin expansion ratio for aesthetics.
#'
#' @return A \code{ggplot} object representing the Venn diagram.
#' 
#' @importFrom ggvenn ggvenn
#' @importFrom ggplot2 scale_x_continuous expansion
#' @importFrom dplyr setdiff
#' @export
#' 
ggvenn_plot <- function(venn_df,
                          col_names = NULL,
                          set.color = c("#E41A1C", "#1E90FF", 
                                        "#FF8C00", "#4DAF4A", "#75cbdc"),
                          set.name.color = "black",
                          use.color.as.text = TRUE,
                          name.size = 5,
                          text.size = 4,
                          stroke.color = "black",
                          stroke.size = 0.6,
                          show.percentage = FALSE,
                          show.elements = FALSE,
                          digits = 1,
                          expand_ratio = 0.2) {
  
  if (is.null(col_names)) {
    col_names <- colnames(venn_df)[-1]
  }
  
  if (length(col_names) < 2 || length(col_names) > 4) {
    stop("ggvenn_plot only supports 2 to 4 sets. Please use ggupset instead.")
  }
  
  if (!all(col_names %in% colnames(venn_df))) {
    stop("The following columns are not in venn_df: ", paste(dplyr::setdiff(col_names, colnames(venn_df)), collapse = ", "))
  }
  
  # check whether all cols selected are logical
  if (!all(sapply(venn_df[col_names], is.logical))) {
    stop("All selected cols must be logical (TRUE/FALSE).")
  }
  
  set_name_color_final <- if (use.color.as.text) {
    set.color[seq_along(col_names)]
  } else {
    set.name.color
  }
  
  p <- ggvenn::ggvenn(
    data = venn_df,
    columns = col_names,
    show_percentage = show.percentage,
    show_elements = show.elements,
    digits = digits,
    fill_color = set.color[seq_along(col_names)],
    stroke_color = stroke.color,
    stroke_size = stroke.size,
    set_name_color = set_name_color_final,
    set_name_size = name.size,
    text_size = text.size
  ) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = expand_ratio))
  
  return(p)
}
