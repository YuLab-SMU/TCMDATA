#' Subset PPI network by edge score and top-n node degree
#'
#' This function filters a PPI network in two steps:
#' 1. Removes edges with a score below \code{score_cutoff}.
#' 2. (Optional) Keeps only the top \code{n} nodes with the highest degree.
#'
#' @param ppi_obj An igraph object of PPI network.
#' @param n Integer. Number of top-degree nodes to keep. If NULL (default), no degree filtering is performed.
#' @param score_cutoff Numeric. Minimum edge score to keep. Default is 0.7.
#' @param edge_attr Character. The name of the edge attribute representing confidence score. Default is "score" (keep the same as STRING).
#' @param rm_isolates Logical. Whether to remove isolated nodes (degree = 0) after edge filtering. Default is TRUE.
#'
#' @return A subgraph of the original PPI network.
#' @importFrom igraph degree E V subgraph.edges induced_subgraph vcount edge_attr_names
#' @importFrom utils head
#' @examples
#' data(demo_ppi)
#' library(igraph)
#' ppi.subset <- ppi_subset(demo_ppi, score_cutoff = 0.7)
#' str(ppi.subset)
#'
#' @export
ppi_subset <- function(ppi_obj,
                       n = NULL,
                       score_cutoff = 0.7,
                       edge_attr = "score",
                       rm_isolates = TRUE) {

  stopifnot(inherits(ppi_obj, "igraph"))

  # 1. Check edge attribute
  if (!edge_attr %in% igraph::edge_attr_names(ppi_obj)) {
    stop(sprintf("Edge attribute '%s' not found in the graph.", edge_attr))
  }

  # 2. Filter by Edge Score
  # Get the score vector
  scores <- igraph::edge_attr(ppi_obj, edge_attr)

  # Identify edges to keep
  eids_to_keep <- igraph::E(ppi_obj)[scores >= score_cutoff]

  if (length(eids_to_keep) == 0) {
    warning("No edges passed the score cutoff. Returning an empty graph.")
    # Return empty graph with same attributes
    return(igraph::delete_vertices(ppi_obj, igraph::V(ppi_obj)))
  }

  # Create subgraph based on edges
  # delete.vertices = FALSE here to handle isolates manually/explicitly later
  ppi_filtered <- igraph::subgraph.edges(ppi_obj, eids = eids_to_keep, delete.vertices = FALSE)

  # 3. Remove Isolates (Optional but Recommended)
  if (rm_isolates) {
    deg_all <- igraph::degree(ppi_filtered, mode = "all")
    ppi_filtered <- igraph::induced_subgraph(ppi_filtered, vids = which(deg_all > 0))
  }

  if (igraph::vcount(ppi_filtered) == 0) {
    return(ppi_filtered)
  }

  # 4. Filter by Top-N Degree
  if (!is.null(n)) {
    deg <- igraph::degree(ppi_filtered, mode = "all")
    deg_sorted <- sort(deg, decreasing = TRUE)
    top_nodes <- names(deg_sorted)[1:min(n, length(deg_sorted))]

    ppi_filtered <- igraph::induced_subgraph(ppi_filtered, vids = top_nodes)
  }

  return(ppi_filtered)
}


#' Compute some useful metrics for PPI nodes
#' @param g An igraph object containing PPI information.
#' @param weight_attr Character; The attribute of weight in this PPI object.
#' @param seed Numeric; Random seed for function `compute_EPC`. Default is 42.
#' @importFrom igraph E V edge_attr_names edge_attr
#' @importFrom igraph degree strength betweenness closeness
#' @importFrom igraph eigen_centrality page_rank coreness
#' @importFrom igraph transitivity eccentricity articulation_points
#' @return The input \code{igraph} object with additional vertex
#' @examples
#' data(demo_ppi)
#' library(igraph)
#' ppi <- compute_nodeinfo(demo_ppi, weight_attr = "score")
#' str(ppi)
#'
#' @export
compute_nodeinfo <- function(g, weight_attr = "score", seed = 42) {
  stopifnot(inherits(g, "igraph"))

  ## PPI score
  w <- NULL
  if (!is.null(weight_attr) && weight_attr %in% igraph::edge_attr_names(g)) {
    w <- igraph::edge_attr(g, weight_attr)
  }

  w_dist <- NULL
  if (!is.null(w)) {
    w_dist <- 1 / pmax(w, .Machine$double.eps)
  }

  ## compute node info
  message("Calculating degree / strength ...")
  V(g)$degree <- igraph::degree(g, mode = "all")
  V(g)$strength <- if (!is.null(w)) igraph::strength(g, mode = "all", weights = w) else NA_real_

  message("Calculating betweenness (unweighted) ...")
  V(g)$betweenness <- igraph::betweenness(g, directed = FALSE, normalized = TRUE)

  if (!is.null(w_dist)) {
    message("Calculating betweenness (weighted) ...")
    V(g)$betweenness_w <- igraph::betweenness(g, directed = FALSE,
                                              weights = w_dist,
                                              normalized = TRUE)
  } else {
    V(g)$betweenness_w <- NA_real_
  }

  message("Calculating closeness (unweighted) ...")
  V(g)$closeness <- igraph::closeness(g, normalized = TRUE)

  if (!is.null(w_dist)) {
    message("Calculating closeness (weighted) ...")
    V(g)$closeness_w <- igraph::closeness(g, normalized = TRUE,
                                          weights = w_dist)
  } else {
    V(g)$closeness_w <- NA_real_
  }

  message("Calculating eigenvector centrality ...")
  V(g)$eigen_centrality <- igraph::eigen_centrality(g, weights = w)$vector

  message("Calculating PageRank ...")
  V(g)$pagerank <- igraph::page_rank(g, weights = w)$vector

  message("Calculating coreness (k-core) ...")
  V(g)$coreness <- igraph::coreness(g, mode = "all")

  message("Calculating local clustering coefficient ...")
  V(g)$clustering_coef <- igraph::transitivity(g, type = "local", isolates = "zero")

  message("Calculating eccentricity ...")
  V(g)$eccentricity <- igraph::eccentricity(g)

  message("Detecting articulation points ...")
  art <- igraph::articulation_points(g)
  V(g)$is_articulation <- FALSE
  V(g)$is_articulation[art] <- TRUE

  ## compute MCC
  message("Calculating MCC ...")
  g <- compute_MCC(g)

  ## MNC
  message("Calculating MNC ...")
  g <- compute_MNC(g)

  ## DMNC
  message("Calculating DMNC ...")
  g <- compute_DMNC(g)

  ## BottleNeck (BN)
  message("Calculating BottleNeck (BN) ...")
  g <- compute_BN(g)

  ## Radiality
  message("Calculating radiality ...")
  g <- compute_radiality(g)

  ## Stress
  message("Calculating Stress ...")
  g <- compute_Stress(g)

  ## EPC
  message("Calculating EPC ...")
  g <- compute_EPC(g, seed = seed)

  return(g)
}


#' Rank PPI nodes by integrated network centrality
#'
#' @param g An igraph object that has already been processed by \code{compute_nodeinfo()} or have added node metrics manually.
#' @param metrics Character vector; which vertex attributes to use for scoring. Defaults are the same as Cytohubba.
#' @param weights Numeric vector of the same length as metrics; relative weights for each metric. If NULL (default), all metrics are equally weighted.
#' @param use_weight Logical; whether use weighted metrics(beweenness and closeness) instead. Default is TRUE.
#' @param na_rm Logical; if TRUE, NAs are ignored in normalization (set to 0.5).
#'
#' @importFrom igraph vertex_attr V
#' @return A list with:
#'   \item{graph}{igraph object with added vertex attributes \code{Score_network} and \code{Rank_network}.}
#'   \item{table}{data.frame with node-level metrics and scores, sorted by \code{Rank_network}.}
#'
#' @examples
#' data(demo_ppi)
#' library(igraph)
#' ppi <- compute_nodeinfo(demo_ppi)
#' rank_res <- rank_ppi_nodes(ppi)
#'
#' # update igraph object and extract rank info
#' ppi <- rank_res[[1]]
#' rank_df <- rank_res[[2]]
#' print(rank_df)
#'
#' @export
rank_ppi_nodes <- function(g,
                           metrics = c(
                             "degree",
                             "betweenness",
                             "closeness",
                             "eccentricity",
                             "radiality",
                             "Stress",
                             "MCC",
                             "MNC",
                             "DMNC",
                             "BN",
                             "EPC"),
                           weights = NULL,
                           use_weight = TRUE,
                           na_rm = TRUE) {
  stopifnot(inherits(g, "igraph"))

  available_metrics <- names(igraph::vertex_attr(g))
  message("Available metrics in graph: ", paste(available_metrics, collapse=", "))

  if (use_weight){
    metrics <- c("degree",
                 "betweenness_w",
                 "closeness_w",
                 "eccentricity",
                 "radiality",
                 "Stress",
                 "MCC",
                 "MNC",
                 "DMNC",
                 "BN",
                 "EPC")
  }

  vdat <- igraph::vertex_attr(g)
  df <- as.data.frame(vdat, stringsAsFactors = FALSE)

  metrics <- intersect(metrics, names(df))
  if (length(metrics) == 0L) {
    stop("None of the requested metrics are present in vertex attributes.")
  }
  message("Using the following metrics for ranking: ",
          paste(metrics, collapse = ", "))

  # default: each metrics are equally weighted
  if (is.null(weights)) {
    weights <- rep(1, length(metrics))
  }
  if (length(weights) != length(metrics)) {
    stop("Length of 'weights' must match length of 'metrics'.")
  }
  # normalize weights
  weights <- weights / sum(weights)

  M_raw  <- df[ , metrics, drop = FALSE]
  M_norm <- as.data.frame(lapply(M_raw, norm01))

  # calculate total scores
  Score_network <- as.numeric(as.matrix(M_norm) %*% weights)

  # rank the nodes
  Rank_network <- rank(-Score_network, ties.method = "min")
  df$Score_network <- Score_network
  df$Rank_network  <- Rank_network

  igraph::V(g)$Score_network <- Score_network
  igraph::V(g)$Rank_network  <- Rank_network

  df_out <- df[order(df$Rank_network), ]

  return(list(graph = g, table = df_out))
}


#' Normalize a numeric vector to the range \eqn{[0, 1]}.
#'
#' A simple and robust 0–1 scaler. Handles NA, Inf, or constant values by
#' returning 0.5 for those entries.
#'
#' @param x Numeric vector.
#' @param na_rm Logical; whether to ignore non-finite values when computing
#'   the range. Defaults to \code{TRUE}.
#'
#' @return Numeric vector which has been scaled.
#'
#' @keywords internal
norm01 <- function(x, na_rm = TRUE) {
  x <- as.numeric(x)
  if (na_rm) {
    finite_idx <- is.finite(x)
    if (!any(finite_idx)) return(rep(0.5, length(x)))
    rng <- range(x[finite_idx], na.rm = TRUE)
  } else {
    rng <- range(x, na.rm = TRUE)
  }
  if (diff(rng) == 0) {
    return(rep(0.5, length(x)))
  } else {
    out <- (x - rng[1]) / diff(rng)
    out[!is.finite(out)] <- 0.5
    return(out)
  }
}

#' compute MCC (maximum clique centrality) of PPI network.
#' @references
#' Chin, C.H., Chen, S.H., Wu, H.H., Ho, C.W., Ko, M.T., & Lin, C.Y. (2014).
#' cytoHubba: identifying hub objects and sub-networks from complex interactome.
#' \emph{BMC Systems Biology}, 8(Suppl 4), S11.
#' https://doi.org/10.1186/1752-0509-8-S4-S11
#' @param graph A igraph object from PPI network.
#' @importFrom igraph max_cliques vcount V V<-
#' @importFrom stats setNames
#'
#' @return An igraph object containing MCC
#' @keywords internal
compute_MCC <- function(graph) {

  max_cliques_list <- igraph::max_cliques(graph)
  mcc <- stats::setNames(rep(0, igraph::vcount(graph)), igraph::V(graph)$name)

  for (clq in max_cliques_list) {
    size <- length(clq)
    weight <- factorial(size - 1)
    mcc[igraph::V(graph)[clq]$name] <- mcc[igraph::V(graph)[clq]$name] + weight
  }

  igraph::V(graph)$MCC <- mcc

  return(graph)
}


#' calculate MNC (Maximum Neighborhood Component) of PPI network
#' @param g An igraph object of PPI network
#' @importFrom igraph V neighbors induced_subgraph components
#' @return An igraph object containing MNC score
#' @references
#' Chin, C.H., Chen, S.H., Wu, H.H., Ho, C.W., Ko, M.T., & Lin, C.Y. (2014).
#' cytoHubba: identifying hub objects and sub-networks from complex interactome.
#' \emph{BMC Systems Biology}, 8(Suppl 4), S11.
#' https://doi.org/10.1186/1752-0509-8-S4-S11
#' @keywords internal
compute_MNC <- function(g) {
  stopifnot(inherits(g, "igraph"))
  igraph::V(g)$MNC <- 0

  for (v in igraph::V(g)) {
    ## find all the neighbors of v firstly
    nbr <- igraph::neighbors(g, v)
    if (length(nbr) < 2) {
      igraph::V(g)$MNC[v] <- length(nbr)
      next
    }
    ## use this neighbors to construct subgraph
    subg <- igraph::induced_subgraph(g, nbr)
    ## find the max connected components of the subgraph
    comp <- igraph::components(subg)$csize
    igraph::V(g)$MNC[v] <- max(comp)
  }

  return(g)
}


#' calculate DMNC (Density of Maximum Neighborhood Component) of PPI network
#' @param g An igraph object of PPI network
#' @param alpha Numeric exponent used to penalize large neighborhood components. The original cytoHubba implementation uses \code{alpha = 1.7}.
#' @importFrom igraph V neighbors induced_subgraph components vcount ecount
#' @return An igraph object containing DMNC score
#' @references
#' Chin, C.H., Chen, S.H., Wu, H.H., Ho, C.W., Ko, M.T., & Lin, C.Y. (2014).
#' cytoHubba: identifying hub objects and sub-networks from complex interactome.
#' \emph{BMC Systems Biology}, 8(Suppl 4), S11.
#' https://doi.org/10.1186/1752-0509-8-S4-S11
#' @keywords internal
compute_DMNC <- function(g, alpha = 1.7) {
  igraph::V(g)$DMNC <- 0

  for (v in igraph::V(g)) {
    nbr <- igraph::neighbors(g, v)
    if (length(nbr) < 2) {
      V(g)$DMNC[v] <- 0
      next
    }

    subg <- igraph::induced_subgraph(g, nbr)
    comp <- igraph::components(subg)
    largest_comp_vertices <- which(comp$membership == which.max(comp$csize))

    sub2 <- igraph::induced_subgraph(subg, largest_comp_vertices)

    nV <- igraph::vcount(sub2)
    nE <- igraph::ecount(sub2)

    igraph::V(g)$DMNC[v] <- nE / (nV ^ alpha)
  }

  return(g)
}


#' Calculate BottleNeck (BN) score for PPI network
#' @description
#' BottleNeck (BN) identifies vertices that frequently serve as bottlenecks
#' in shortest-path trees rooted at all vertices, following the definition
#' in cytoHubba (Chin et al., 2014).
#'
#' @param g An \code{igraph} object (typically an undirected, unweighted PPI network).
#'
#' @return The input graph \code{g} with an additional vertex attribute
#'   \code{BN} storing the BottleNeck score for each node.
#'
#' @references
#' Chin, C.H., Chen, S.H., Wu, H.H., Ho, C.W., Ko, M.T., & Lin, C.Y. (2014).
#' cytoHubba: identifying hub objects and sub-networks from complex interactome.
#' \emph{BMC Systems Biology}, 8(Suppl 4), S11.
#'
#' @importFrom igraph vcount bfs V
#' @keywords internal
compute_BN <- function(g) {
  stopifnot(inherits(g, "igraph"))

  n <- igraph::vcount(g)
  BN <- numeric(n)

  for (s in seq_len(n)) {

    bfs_res <- igraph::bfs(
      graph       = g,
      root        = s,
      mode        = "all",
      unreachable = FALSE,
      father      = TRUE,
      dist        = TRUE
    )

    dist   <- as.numeric(bfs_res$dist)
    father <- as.numeric(bfs_res$father)

    reachable_indices <- which(!is.na(dist))
    nT <- length(reachable_indices)

    if (nT < 2) next

    subtree_size <- integer(n)
    subtree_size[reachable_indices] <- 1L

    ord <- reachable_indices[order(dist[reachable_indices], decreasing = TRUE)]

    for (v in ord) {
      p <- father[v]
      if (!is.na(p) && p > 0 && p != v) {
        subtree_size[p] <- subtree_size[p] + subtree_size[v]
      }
    }
    thr <- nT / 4
    candidates <- reachable_indices[subtree_size[reachable_indices] > thr]

    if (length(candidates) > 0) {
      BN[candidates] <- BN[candidates] + 1
    }
  }

  igraph::V(g)$BN <- BN
  return(g)
}

#' Calculate Radiality Centrality for PPI network
#'
#' This function calculates the radiality centrality for each node in a graph.
#' It implements the approach often used in tools like CytoHubba to handle
#' disconnected graphs by penalizing unreachable nodes based on the graph diameter.
#'
#' @param g An \code{igraph} object.
#' @return The input \code{igraph} object with a new vertex attribute \code{radiality}.
#' @references
#' Chin, C.H., Chen, S.H., Wu, H.H., Ho, C.W., Ko, M.T., & Lin, C.Y. (2014).
#' cytoHubba: identifying hub objects and sub-networks from complex interactome.
#' \emph{BMC Systems Biology}, 8(Suppl 4), S11.
#' @importFrom igraph vcount distances V
#' @keywords internal
compute_radiality <- function(g) {
  stopifnot(inherits(g, "igraph"))

  n <- igraph::vcount(g)
  if (n < 2) {
    igraph::V(g)$radiality <- ifelse(n == 1, 1, 0)
    return(g)
  }

  dist_mat <- igraph::distances(g, mode = "all")
  finite_dists <- dist_mat[is.finite(dist_mat)]

  if (length(finite_dists) == 0) {
    D <- 0
  } else {
    D <- max(finite_dists)
  }

  dist_mat[!is.finite(dist_mat)] <- D + 1

  # Calculate Radiality
  # Numerator: Sum of (D + 1 - distance) for all other nodes
  # We divide by (n - 1) to normalize the value (average closeness gain).
  # This normalization does not change the ranking but makes values comparable.
  numerator <- rowSums((D + 1) - dist_mat)
  radiality <- numerator / (n - 1)

  # Assign to graph
  igraph::V(g)$radiality <- radiality

  return(g)
}


#' Calculate Stress Centrality for PPI network
#'
#' This function calculates the Stress Centrality for each vertex in the graph.
#' Stress centrality of a node v is defined as the number of shortest paths
#' passing through v.
#'
#' @details
#' The Stress Centrality \eqn{C_{str}(v)} is calculated as:
#' \deqn{C_{str}(v) = \sum_{s \neq v \neq t} \sigma_{st}(v)}
#' where \eqn{\sigma_{st}(v)} is the number of shortest paths from node \eqn{s}
#' to node \eqn{t} which use the node \eqn{v}.
#'
#' This implementation uses a variation of Brandes' algorithm. Note that for
#' large graphs, this computation can be time-consuming as it essentially
#' involves all-pairs shortest path calculations with complexity \eqn{O(VE)}.
#'
#' @param g An \code{igraph} object.
#' @importFrom igraph as_adj_list vcount V
#' @references
#' Chin, C.H., Chen, S.H., Wu, H.H., Ho, C.W., Ko, M.T., & Lin, C.Y. (2014).
#' cytoHubba: identifying hub objects and sub-networks from complex interactome.
#' \emph{BMC Systems Biology}, 8(Suppl 4), S11.
#' @return The input \code{igraph} object with a new vertex attribute \code{Stress}.
#' @keywords internal
compute_Stress <- function(g) {
  stopifnot(inherits(g, "igraph"))

  adj_list <- igraph::as_adj_list(g, mode = "all")
  n <- igraph::vcount(g)
  stress_score <- numeric(n)

  for (s in seq_len(n)) {
    S <- integer()
    P <- vector("list", n)
    sigma <- numeric(n)
    dist <- rep(-1, n)

    sigma[s] <- 1
    dist[s] <- 0

    Q <- c(s)

    while (length(Q) > 0) {
      v <- Q[1]
      Q <- Q[-1]
      S <- c(S, v)

      neighbors <- adj_list[[v]]
      for (w in neighbors) {
        if (dist[w] < 0) {
          Q <- c(Q, w)
          dist[w] <- dist[v] + 1
        }
        if (dist[w] == dist[v] + 1) {
          sigma[w] <- sigma[w] + sigma[v]
          P[[w]] <- c(P[[w]], v)
        }
      }
    }

    gamma <- numeric(n)

    for (i in length(S):1) {
      w <- S[i]
      total_paths_from_w <- 1 + gamma[w]
      for (v in P[[w]]) {
        gamma[v] <- gamma[v] + total_paths_from_w
      }
    }

    valid_nodes <- which(dist >= 0)
    valid_nodes <- valid_nodes[valid_nodes != s]

    if (length(valid_nodes) > 0) {
      stress_score[valid_nodes] <- stress_score[valid_nodes] + (sigma[valid_nodes] * gamma[valid_nodes])
    }
  }

  igraph::V(g)$Stress <- stress_score
  return(g)
}

#' Calculate Edge Percolated Component (EPC) for a PPI network
#'
#' EPC centrality is computed by simulating random edge percolation
#' (random removal of edges) and measuring how well each node remains
#' connected across multiple reduced networks, following the definition
#' in cytoHubba.
#'
#' @param g An \code{igraph} object representing a PPI network.
#' @param threshold Numeric in \code{[0, 1]}. Edges with a random value
#'   below this threshold are removed in each simulation. Default is \code{0.5}.
#' @param n_iter Integer. Number of percolation simulations to run.
#'   Default is \code{1000}.
#' @param seed Integer. Random seed for reproducibility. Default is \code{42}.
#'
#' @return The input \code{igraph} object with a new vertex attribute
#'   \code{EPC} storing the EPC score for each node.
#'
#' @importFrom igraph vcount ecount subgraph_from_edges components V
#' @importFrom stats runif
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @references
#' Chin, C.H., Chen, S.H., Wu, H.H., Ho, C.W., Ko, M.T., & Lin, C.Y. (2014).
#' cytoHubba: identifying hub objects and sub-networks from complex interactome.
#' \emph{BMC Systems Biology}, 8(Suppl 4), S11.
#' @keywords internal
compute_EPC <- function(g, threshold = 0.5, n_iter = 1000, seed = 42) {
  stopifnot(inherits(g, "igraph"))

  set.seed(seed)

  n_vertices <- igraph::vcount(g)
  n_edges    <- igraph::ecount(g)

  # accumulate component sizes for each vertex across simulations
  sum_comp_sizes <- numeric(n_vertices)

  pb <- utils::txtProgressBar(min = 0, max = n_iter, style = 3)

  for (k in seq_len(n_iter)) {
    random_vals   <- stats::runif(n_edges)
    edges_to_keep <- which(random_vals >= threshold)

    g_reduced <- igraph::subgraph_from_edges(
      g,
      eids = edges_to_keep,
      delete.vertices = FALSE
    )

    comps <- igraph::components(g_reduced)

    current_sizes   <- comps$csize[comps$membership]
    sum_comp_sizes  <- sum_comp_sizes + current_sizes

    utils::setTxtProgressBar(pb, k)
  }
  close(pb)

  epc_scores <- sum_comp_sizes / n_iter

  igraph::V(g)$EPC <- epc_scores
  return(g)
}
