#' Prepare an igraph object from a Herb–Molecule–Target data frame
#'
#' @param df A \code{data.frame} with three columns representing herb–molecule–target relationships.
#' @param n Integer. Optional. The number of rows to use from df. Default is 'nrow(df)'
#' @param herb_col Character. The column name corresponding to herbs. Default is \code{"herb"}.
#' @param molecule_col Character. The column name corresponding to molecules/compounds. Default is \code{"molecule"}.
#' @param target_col Character. The column name corresponding to targets. Default is \code{"target"}.
#' @param compute_metrics Logical. Whether to compute and attach node-level metrics. Default is TRUE.
#'
#' @return A directed \code{igraph}
#' @importFrom dplyr select bind_rows all_of case_when
#' @importFrom igraph graph_from_data_frame degree eigen_centrality V V<-
#' @importFrom utils head
#'
#' @examples
#' library(igraph)
#' library(ggplot2)
#' herbs <- c("灵芝")
#' lz <- search_herb(herb = herbs, type = "Herb_cn_name")[1:40, ]
#' res <- prepare_herb_graph(lz)
#'
#' # visualization
#' p <- ggtangle::ggplot(res, layout = "kk") +
#'  ggtangle::geom_edge(alpha = 0.4) +
#'  ggrepel::geom_text_repel(aes(label = name), size = 3.5, max.overlaps = 20) +
#'  geom_point(aes(color = type, size = centrality)) +
#'  scale_color_manual(values = c(
#'    "Herb" = "#E41A1C", "Molecule" = "#377EB8", "Target" = "#4DAF4A")) +
#'  theme_void()
#'  print(p)
#'
#' @export
prepare_herb_graph <- function(
    df,
    n = NULL,
    herb_col = "herb",
    molecule_col = "molecule",
    target_col = "target",
    compute_metrics = TRUE) {

  if (is.null(n)){
    n <- nrow(df)
  }
  else{
    n <- min(n, nrow(df))
  }

  stopifnot(is.data.frame(df))
  required_cols <- c(herb_col, molecule_col, target_col)
  if (!all(required_cols %in% colnames(df))) {
    stop("Input data frame must contain columns: ", paste(required_cols, collapse = ", "))
  }

  df <- df %>% utils::head(n)

  edges1 <- df %>% select(from = all_of(herb_col), to = all_of(molecule_col))
  edges2 <- df %>% select(from = all_of(molecule_col), to = all_of(target_col))
  edges <- bind_rows(edges1, edges2)

  g <- igraph::graph_from_data_frame(edges, directed = TRUE)

  if (compute_metrics) {
    V(g)$degree <- igraph::degree(g, mode = "all")
    V(g)$centrality <- igraph::eigen_centrality(g)$vector
    V(g)$betweenness <- igraph::betweenness(g, directed = TRUE)
    V(g)$closeness <- igraph::closeness(g, mode = "all")
    V(g)$pagerank <- igraph::page_rank(g)$vector
  }

  node_names <- V(g)$name

  V(g)$type <- case_when(
    node_names %in% df[[herb_col]] ~ "Herb",
    node_names %in% df[[molecule_col]] ~ "Molecule",
    node_names %in% df[[target_col]] ~ "Target",
    TRUE ~ "Other"
  )

  return(g)
}
