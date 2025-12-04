#' Compute MCL Clustering for PPI Network
#'
#' This function performs Markov Cluster Algorithm (MCL) on an igraph object
#' and assigns the cluster labels to the nodes.
#'
#' @param g An \code{igraph} object.
#' @param inflation Numeric. The inflation parameter controls the granularity of clusters.
#'   Common values for PPI: 2.0 - 3.0.
#'   Larger values (e.g., 3.0) produce smaller, tighter clusters.
#'   Smaller values (e.g., 1.5) produce larger, coarser clusters.
#'   Default is 2.5.
#' @param addLoops Logical; Self-loops with weight 1 are added to each vertex of g when TRUE; Default is TRUE (necessary).
#' @param ... Additional parameters in function `mcl()`.
#' @importFrom igraph as_adjacency_matrix V
#' @importFrom MCL mcl
#'
#' @return The input \code{igraph} object with a new vertex attribute \code{mcl_cluster}.
#' @examples
#' data(demo_ppi)
#' library(igraph)
#' ppi <- run_MCL(demo_ppi, inflation = 2.5)
#' head(V(ppi)$mcl_cluster)
#'
#' @export
run_MCL <- function(g, inflation = 2.5, addLoops = TRUE, ...) {

  stopifnot(inherits(g, "igraph"))
  message(sprintf("Running MCL with inflation = %.1f ...", inflation))

  adj_mat <- igraph::as_adjacency_matrix(g, sparse = TRUE)

  ## MCL running
  mcl_res <- MCL::mcl(x = adj_mat, addLoops = addLoops, inflation = inflation, ...)

  igraph::V(g)$mcl_cluster <- as.factor(mcl_res$Cluster)
  message(sprintf("Done! Identified %d modules (Iterations: %d).",
                  mcl_res$K, mcl_res$n.iterations))

  return(g)
}


#' Compute Louvain Clustering for PPI Network
#'
#' This function performs Louvain clustering (multi-level modularity optimization)
#' on an igraph object. It is suitable for detecting larger, macro-scale functional
#' modules in PPI networks.
#'
#' @param g An \code{igraph} object.
#' @param resolution Numeric. The resolution parameter controls the granularity of clusters. Default is 1.0 (standard modularity).
#' @param weights Numeric vector or NULL. Edge weights to use for clustering.
#'   If NULL (default), the function attempts to use the 'weight' or 'score' edge attribute.
#'   Set to NA to perform unweighted clustering.
#'
#' @return The input \code{igraph} object with a new vertex attribute \code{louvain_cluster}.
#' @importFrom igraph cluster_louvain V edge_attr_names E membership
#'
#' @examples
#' data(demo_ppi)
#' library(igraph)
#' ppi <- run_louvain(demo_ppi, resolution = 1)
#' head(V(ppi)$louvain_cluster)
#'
#' @export
run_louvain <- function(g, resolution = 1.0, weights = NULL) {

  stopifnot(inherits(g, "igraph"))

  # PPI networks often have 'score' or 'weight'. Using them improves accuracy.
  if (is.null(weights)) {
    edge_attrs <- igraph::edge_attr_names(g)
    if ("weight" %in% edge_attrs) {
      weights <- igraph::E(g)$weight
      message("Using edge attribute 'weight' for clustering.")
    } else if ("score" %in% edge_attrs) {
      weights <- igraph::E(g)$score
      message("Using edge attribute 'score' for clustering.")
    } else {
      weights <- NULL # Unweighted
      message("No weights found. Running unweighted clustering.")
    }
  } else if (identical(weights, NA)) {
    weights <- NULL # Explicitly unweighted
  }

  message(sprintf("Running Louvain (Resolution: %.1f)...", resolution))

  # Run Louvain Algorithm
  louvain_res <- igraph::cluster_louvain(g, weights = weights, resolution = resolution)

  # Assign Cluster Labels
  igraph::V(g)$louvain_cluster <- as.factor(igraph::membership(louvain_res))
  num_clusters <- length(unique(igraph::membership(louvain_res)))
  modularity_score <- max(louvain_res$modularity, na.rm = TRUE)

  message(sprintf("Done! Identified %d modules (Modularity: %.3f).",
                  num_clusters, modularity_score))

  return(g)
}


#' Score and Rank Network Clusters
#'
#' This function evaluates the clusters/modules identified in the graph.
#' It calculates a score for each cluster based on density and size (MCODE style),
#' and returns a ranked data frame.
#'
#' @param g An \code{igraph} object. The graph must have a vertex attribute containing cluster labels.
#' @param cluster_attr Character. The name of the vertex attribute that stores cluster labels. Default is louvain_cluster.
#' @param min_size Integer. Clusters smaller than this size will be ignored. Default is 3.
#'
#' @return A data frame containing cluster statistics, ranked by Score.
#' @importFrom igraph vertex_attr vertex_attr_names induced_subgraph vcount edge_density V ecount
#' @importFrom utils head
#' @examples
#' data(demo_ppi)
#' ppi <- run_louvain(demo_ppi, resolution = 1)
#' louvain_score <- Addclusterscore(ppi, cluster_attr = "louvain_cluster", min_size = 3)
#' head(louvain_score)
#'
#' @export
Addclusterscore <- function(g, cluster_attr = "louvain_cluster", min_size = 3) {

  stopifnot(inherits(g, "igraph"))

  if (!cluster_attr %in% igraph::vertex_attr_names(g)) {
    stop(paste0("Attribute '", cluster_attr, "' not found in graph. Please run clustering first."))
  }

  labels <- igraph::vertex_attr(g, cluster_attr)
  unique_clusters <- unique(labels)
  unique_clusters <- unique_clusters[!is.na(unique_clusters)]
  message(paste("Evaluating", length(unique_clusters), "clusters based on attribute:", cluster_attr))

  stats_list <- lapply(unique_clusters, function(cid) {
    nodes_in_mod <- igraph::V(g)$name[labels == cid]
    n_nodes <- length(nodes_in_mod)

    if (n_nodes < min_size) return(NULL)
    subg <- igraph::induced_subgraph(g, nodes_in_mod)

    # Density = E / (N * (N-1) / 2)
    dens <- igraph::edge_density(subg)

    # Score = Density * N
    score <- dens * n_nodes

    data.frame(
      Cluster_ID = as.character(cid),
      Score      = round(score, 3),
      Nodes      = n_nodes,
      Edges      = igraph::ecount(subg),
      Density    = round(dens, 3),
      Gene_List  = paste(utils::head(nodes_in_mod, 5), collapse = ", "),
      Full_Genes = paste(nodes_in_mod, collapse = ","),
      stringsAsFactors = FALSE
    )
  })

  stats_list <- stats_list[!sapply(stats_list, is.null)]

  if (length(stats_list) == 0) {
    warning("No clusters passed the size filter.")
    return(data.frame())
  }

  df_res <- do.call(rbind, stats_list)
  df_res <- df_res[order(df_res$Score, decreasing = TRUE), ]
  rownames(df_res) <- NULL

  return(df_res)
}

