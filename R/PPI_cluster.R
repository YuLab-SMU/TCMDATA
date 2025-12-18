#' Perform Markov Clustering (MCL) on a Graph
#'
#' @description
#' This function implements the Markov Clustering (MCL) algorithm for detecting
#' communities (clusters) in a graph. MCL simulates random walks within the graph
#' by alternating between two operations: expansion and inflation. It is particularly
#' efficient for biological networks.
#' @param g An \code{igraph} object. The graph to be clustered. It can be directed or undirected.
#' @param inflation Numeric. Controls cluster granularity. Higher values (e.g., > 2) yield smaller, tighter clusters; lower values yield larger clusters. Default is 2.5.
#' @param max_iter Integer. The maximum number of iterations to perform if convergence is not reached. Default is 100.
#' @param pruning Numeric. A threshold for pruning small values in the matrix to zero.
#'   This preserves the sparsity of the matrix and significantly speeds up computation
#'   while saving memory. Default is 1e-5.
#' @param allow1 Logical. If TRUE, clusters with only 1 node are kept as unique clusters. If FALSE,
#' cluster of size 1 are interpreted as background noise and grouped in one cluster. Default is FALSE.
#' @importFrom igraph as_adjacency_matrix is_igraph V
#' @importFrom Matrix Diagonal colSums drop0
#' @importFrom stats setNames
#' @return An igraph object containing MCL clustering labels.
#'
#' @examples
#' library(igraph)
#' g <- make_graph("Zachary")
#' g <- run_MCL(g, inflation = 2.5)
#' print(head(V(g)$MCL_cluster))
#'
#' # Visualize
#' plot(g,
#' vertex.color = V(g)$MCL_cluster,
#' vertex.size = 15,
#' vertex.label = V(g)$name)
#'
#' @export
run_MCL <- function(g,
                    inflation = 2.5,
                    max_iter = 100,
                    pruning = 1e-5,
                    allow1 = FALSE){

  stopifnot(igraph::is_igraph(g))

  ## create adjacency matrix
  adj <- igraph::as_adjacency_matrix(g, sparse = TRUE)

  ## add self-loops
  M <- adj + Matrix::Diagonal(nrow(adj))

  ## scale initially (Column Normalization)
  col_sum <- Matrix::colSums(M)
  col_sum[col_sum == 0] <- 1 # Avoid division by zero for isolated nodes
  M <- Matrix::t(Matrix::t(M) / col_sum)

  ## MCL iteration
  for (i in 1:max_iter){
    M_prev <- M

    # 1. expansion (Matrix Multiplication)
    M <- M %*% M

    # 2. inflation (Element-wise power)
    M <- M ^ inflation

    # 3. pruning (Keep matrix sparse)
    M[M < pruning] <- 0
    M <- Matrix::drop0(M) # Physically remove zeros from storage

    # 4. re-scale (Re-normalize columns)
    col_sum <- Matrix::colSums(M)
    col_sum[col_sum == 0] <- 1
    M <- Matrix::t(Matrix::t(M) / col_sum)

    # 5. check convergence
    diff <- sum((M - M_prev)^2)

    if (diff < 1e-5) {
      cat(sprintf("Converged at iteration %d (Diff: %.2e)\n", i, diff))
      break
    }
  }

  ## get results
  raw_clusters <- apply(M, 2, which.max)

  if (!allow1) {
    cluster_counts <- table(raw_clusters)
    singleton_ids <- as.integer(names(cluster_counts)[cluster_counts == 1])

    if (length(singleton_ids) > 0) {
      raw_clusters[raw_clusters %in% singleton_ids] <- 0
      message(sprintf("Note: Merged %d singleton nodes into background cluster (0).", length(singleton_ids)))
    }
  }

  ## Renumber clusters (Formatting)
  final_clusters <- raw_clusters
  valid_ids <- sort(unique(raw_clusters[raw_clusters != 0]))

  if (length(valid_ids) > 0) {
    mapping <- setNames(seq_along(valid_ids), valid_ids)
    mask_non_zero <- raw_clusters != 0
    final_clusters[mask_non_zero] <- mapping[as.character(raw_clusters[mask_non_zero])]
  }

  n_clusters <- length(valid_ids)
  message(paste("Result: Identified", n_clusters, "valid clusters (excluding noise)."))

  # Assign to graph
  igraph::V(g)$MCL_cluster <- as.integer(final_clusters)

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

