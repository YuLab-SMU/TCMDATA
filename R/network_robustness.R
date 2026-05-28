#' Evaluate PPI network robustness via drug attack model
#'
#' @description
#' Simulates a targeted node-knockout attack on a weighted PPI network and
#' evaluates the statistical significance of network disruption using a
#' permutation test.  Four topological metrics are tracked (ASPL, AD, DC, CC),
#' exactly following the methodology of Xi et al. (2022).
#'
#' @param g An undirected \code{igraph} object representing the PPI network.
#' @param targets Character vector of node names to knock out (drug targets).
#' @param n_perm Integer. Number of permutation iterations for the null
#'   distribution.  Default is \code{100}, matching the original paper.
#' @param metrics Character vector of whole-network metrics to track. Defaults
#'   to the original \code{ppi_knock()} metrics \code{"ASPL"}, \code{"AD"},
#'   \code{"DC"}, and \code{"CC"}. Additional supported exploratory metrics are
#'   \code{"density"}, \code{"transitivity"}, \code{"diameter"},
#'   \code{"components"}, and \code{"largest_component_fraction"}. The
#'   integrated \code{Total_Score} is only computed when all original four
#'   metrics are present.
#' @param weight_attr Character. Edge attribute name for confidence weights.
#'   Default is \code{"score"}.
#' @param rewire_niter Integer. Multiplier for edge-swap attempts per rewiring
#'   step (\code{ecount(g) * rewire_niter} swaps).  Default is \code{10}.
#' @param seed Integer. Random seed for reproducibility. Default is \code{42}.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{Summary}{A data.frame with one row per metric (\code{ASPL},
#'     \code{AD}, \code{DC}, \code{CC}) and columns \code{Baseline},
#'     \code{Post_KO}, \code{Raw_RI}, \code{Mu_Random}, \code{Sd_Random},
#'     \code{Normalized_RI} (Z-score), and \code{Pvalue}
#'     (two-sided, normal approximation).}
#'   \item{Total_Score}{Numeric. The integrated disruption score
#'     \eqn{Z_{ASPL} - Z_{AD} - Z_{DC} - Z_{CC}}.}
#'   \item{Total_Pvalue}{Numeric. Two-sided p-value for \code{Total_Score},
#'     derived by comparing the real combined RI against the permutation
#'     null distribution of combined RIs.}
#' }
#'
#' @references
#' Xi, Y., et al. (2022). Exploration of the Specific Pathology of HXMM Tablet
#' Against Retinal Injury Based on Drug Attack Model to Network Robustness.
#' \emph{Frontiers in Pharmacology}, 13, 826535.
#' \doi{10.3389/fphar.2022.826535}
#'
#' @importFrom igraph V vcount ecount edge_attr set_edge_attr delete_vertices
#'   components mean_distance rewire keeping_degseq degree centr_degree
#'   closeness edge_density transitivity diameter
#' @importFrom stats sd pnorm
#'
#' @examples
#' \dontrun{
#' data(demo_ppi)
#' targets <- "IL6"
#' res <- ppi_knock(demo_ppi, targets, n_perm = 100)
#' print(res$Summary)
#' cat("Total Score:", res$Total_Score, "\n")
#' cat("Total P-value:", res$Total_Pvalue, "\n")
#' }
#'
#' @export
ppi_knock <- function(g, targets,
                      n_perm = 100L,
                      metrics = c("ASPL", "AD", "DC", "CC"),
                      weight_attr = "score",
                      rewire_niter = 10L,
                      seed = 42L) {

  stopifnot(inherits(g, "igraph"))

  set.seed(seed)

  node_names <- .get_node_names(g)
  missing <- setdiff(targets, node_names)
  if (length(missing) > 0)
    stop("Targets not found in network: ", paste(missing, collapse = ", "))

  n_targets <- length(targets)
  n_edges   <- igraph::ecount(g)
  weights   <- igraph::edge_attr(g, weight_attr)

  if (is.null(weights))
    stop(sprintf("Edge attribute '%s' not found.", weight_attr))

  metric_names <- .validate_network_robustness_metrics(metrics)
  n_metrics <- length(metric_names)

  ## Stage 1: Baseline metrics
  m_base <- .network_metrics(g, weight_attr, metrics = metric_names)

  ## Stage 2: Real attack — remove targets, compute RI (relative change)
  g_ko   <- igraph::delete_vertices(g, targets)
  m_ko   <- .network_metrics(g_ko, weight_attr, metrics = metric_names)
  raw_ri <- .safe_ri(m_ko, m_base)

  ## Stage 3: Permutation null distribution
  rand_ri <- vapply(seq_len(n_perm), function(i) {
    g_rand <- igraph::rewire(g, igraph::keeping_degseq(niter = n_edges * rewire_niter))
    g_rand <- igraph::set_edge_attr(g_rand, weight_attr,
                                    value = sample(weights))
    m_rand_base  <- .network_metrics(g_rand, weight_attr, metrics = metric_names)
    rand_targets <- sample(node_names, n_targets)
    g_rand_ko    <- igraph::delete_vertices(g_rand, rand_targets)
    m_rand_ko    <- .network_metrics(g_rand_ko, weight_attr, metrics = metric_names)
    .safe_ri(m_rand_ko, m_rand_base)
  }, numeric(n_metrics))
  rownames(rand_ri) <- metric_names

  mu_rand <- rowMeans(rand_ri)
  sd_rand <- apply(rand_ri, 1, stats::sd)

  ## Stage 4: Z-score normalisation + p-value (normal approximation)
  sd_safe       <- ifelse(sd_rand == 0, 1, sd_rand)
  normalized_ri <- (raw_ri - mu_rand) / sd_safe

  # Two-sided p-value via Z-score (paper's 95% CI approach)
  p_zscore <- 2 * stats::pnorm(-abs(normalized_ri))

  ## Total Score = Z_ASPL - Z_AD - Z_DC - Z_CC
  if (all(c("ASPL", "AD", "DC", "CC") %in% metric_names)) {
    total_score <- normalized_ri[["ASPL"]] - normalized_ri[["AD"]] -
      normalized_ri[["DC"]] - normalized_ri[["CC"]]

    ## Total Score p-value from permutation null distribution
    ## Combine per-permutation RIs using the same sign convention
    null_combined <- rand_ri["ASPL", ] - rand_ri["AD", ] - rand_ri["DC", ] - rand_ri["CC", ]
    real_combined <- raw_ri[["ASPL"]] - raw_ri[["AD"]] - raw_ri[["DC"]] - raw_ri[["CC"]]
    mu_comb  <- mean(null_combined)
    sd_comb  <- stats::sd(null_combined)
    sd_comb  <- ifelse(sd_comb == 0, 1, sd_comb)
    z_total  <- (real_combined - mu_comb) / sd_comb
    p_total  <- 2 * stats::pnorm(-abs(z_total))
  } else {
    total_score <- NA_real_
    p_total <- NA_real_
  }

  ## Summary table
  summary_df <- data.frame(
    Metric        = metric_names,
    Baseline      = m_base,
    Post_KO       = m_ko,
    Raw_RI        = raw_ri,
    Mu_Random     = mu_rand,
    Sd_Random     = sd_rand,
    Normalized_RI = normalized_ri,
    Pvalue        = p_zscore,
    row.names     = NULL
  )

  res <- list(Summary = summary_df, Total_Score = unname(total_score), Total_Pvalue = unname(p_total))
  return(res)
}

#' Estimate node-level PPI perturbation after virtual protein knockout
#'
#' @description
#' Extends \code{\link{ppi_knock}} from whole-network robustness to node-level
#' perturbation. After removing one or more protein nodes, the function compares
#' selected topological metrics for every remaining protein before and after
#' knockout, then evaluates whether the observed metric changes are larger than
#' expected from the same random-attack background used by \code{ppi_knock}:
#' degree-sequence rewiring, edge-weight shuffling, and random removal of the
#' same number of nodes.
#'
#' This analysis should be interpreted as PPI topological sensitivity, not as
#' direct evidence of expression change or regulatory direction.
#'
#' @param g An undirected \code{igraph} object representing the PPI network.
#' @param targets Character vector of node names to knock out.
#' @param metrics Character vector of node-level topological metrics to track.
#'   Defaults to the same metric family used by \code{\link{ppi_knock}}:
#'   \code{"ASPL"}, \code{"AD"}, \code{"DC"}, and \code{"CC"}. At node level,
#'   \code{"ASPL"} is the weighted average shortest path length from a protein
#'   to all other reachable proteins, \code{"AD"} is node degree, \code{"DC"}
#'   is normalized node degree centrality, and \code{"CC"} is normalized
#'   closeness centrality. Additional supported exploratory metrics include
#'   \code{"degree"}, \code{"strength"}, \code{"betweenness"},
#'   \code{"betweenness_w"}, \code{"closeness"}, \code{"closeness_w"},
#'   \code{"pagerank"}, \code{"eigen_centrality"}, \code{"coreness"},
#'   \code{"clustering_coef"}, \code{"eccentricity"}, \code{"radiality"},
#'   \code{"Stress"}, \code{"MCC"}, \code{"MNC"}, \code{"DMNC"},
#'   \code{"BN"}, and \code{"EPC"}.
#' @param n_perm Integer. Number of random perturbations for the null
#'   distribution.
#' @param weight_attr Character. Edge attribute name for confidence weights.
#' @param rewire_niter Integer. Multiplier for edge-swap attempts per rewiring
#'   step, as in \code{\link{ppi_knock}}.
#' @param seed Integer. Random seed for reproducibility.
#' @param p_adjust_method P-value adjustment method passed to
#'   \code{\link[stats]{p.adjust}}.
#' @param alpha Significance cutoff used to label affected proteins in the
#'   protein-level summary.
#'
#' @return A \code{tcm_ppi_knock_impact} object with:
#' \describe{
#'   \item{impact}{A long-format data frame with one row per affected protein
#'     and metric, containing before/after values, deltas, random background,
#'     Z-scores, empirical and adjusted P-values.}
#'   \item{protein_summary}{A protein-level summary aggregating across metrics.}
#'   \item{params}{Analysis parameters.}
#' }
#'
#' @importFrom igraph V vcount ecount edge_attr edge_attr_names set_edge_attr
#'   delete_vertices rewire keeping_degseq degree strength betweenness closeness
#'   page_rank eigen_centrality coreness transitivity eccentricity distances
#' @importFrom stats pnorm p.adjust sd
#'
#' @examples
#' \dontrun{
#' data(demo_ppi)
#' res <- ppi_knock_impact(
#'   demo_ppi,
#'   targets = "IL6",
#'   n_perm = 100
#' )
#' head(res$impact)
#' head(res$protein_summary)
#' }
#'
#' @export
ppi_knock_impact <- function(g,
                             targets,
                             metrics = c("ASPL", "AD", "DC", "CC"),
                             n_perm = 100L,
                             weight_attr = "score",
                             rewire_niter = 10L,
                             seed = 42L,
                             p_adjust_method = "BH",
                             alpha = 0.05) {
  stopifnot(inherits(g, "igraph"))
  set.seed(seed)

  node_names <- .get_node_names(g)
  targets <- unique(as.character(targets))
  missing <- setdiff(targets, node_names)
  if (length(missing) > 0L) {
    stop("Targets not found in network: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  if (length(targets) >= length(node_names)) {
    stop("At least one non-target node must remain after knockout.", call. = FALSE)
  }

  metrics <- .validate_node_impact_metrics(metrics)
  n_perm <- as.integer(n_perm)
  if (n_perm < 1L) stop("n_perm must be at least 1.", call. = FALSE)

  weights <- igraph::edge_attr(g, weight_attr)
  if (is.null(weights)) {
    stop(sprintf("Edge attribute '%s' not found.", weight_attr), call. = FALSE)
  }

  n_targets <- length(targets)
  n_edges <- igraph::ecount(g)
  remaining <- setdiff(node_names, targets)

  base_df <- .node_metrics_selected(g, metrics = metrics, weight_attr = weight_attr)
  g_ko <- igraph::delete_vertices(g, targets)
  after_df <- .node_metrics_selected(g_ko, metrics = metrics, weight_attr = weight_attr)

  base_mat <- as.matrix(base_df[match(remaining, base_df$name), metrics, drop = FALSE])
  after_mat <- as.matrix(after_df[match(remaining, after_df$name), metrics, drop = FALSE])
  rownames(base_mat) <- remaining
  rownames(after_mat) <- remaining

  observed_delta <- after_mat - base_mat
  observed_ri <- .safe_ri(after_mat, base_mat)
  dimnames(observed_ri) <- dimnames(base_mat)
  metric_sign <- .node_impact_metric_sign(metrics)
  observed_impact <- sweep(observed_ri, 2, metric_sign, `*`)
  dimnames(observed_impact) <- dimnames(base_mat)

  random_impact <- array(
    NA_real_,
    dim = c(length(remaining), length(metrics), n_perm),
    dimnames = list(remaining, metrics, paste0("perm", seq_len(n_perm)))
  )

  for (i in seq_len(n_perm)) {
    g_rand <- igraph::rewire(g, igraph::keeping_degseq(niter = n_edges * rewire_niter))
    g_rand <- igraph::set_edge_attr(g_rand, weight_attr, value = sample(weights))

    rand_base_df <- .node_metrics_selected(g_rand, metrics = metrics, weight_attr = weight_attr)
    rand_targets <- sample(node_names, n_targets)
    g_rand_ko <- igraph::delete_vertices(g_rand, rand_targets)
    rand_after_df <- .node_metrics_selected(g_rand_ko, metrics = metrics, weight_attr = weight_attr)

    valid_nodes <- intersect(remaining, setdiff(node_names, rand_targets))
    if (length(valid_nodes) == 0L) next

    rand_base_mat <- as.matrix(rand_base_df[match(valid_nodes, rand_base_df$name), metrics, drop = FALSE])
    rand_after_mat <- as.matrix(rand_after_df[match(valid_nodes, rand_after_df$name), metrics, drop = FALSE])
    rand_ri <- .safe_ri(rand_after_mat, rand_base_mat)
    rand_impact <- sweep(rand_ri, 2, metric_sign, `*`)
    random_impact[valid_nodes, metrics, i] <- rand_impact
  }

  random_mean <- apply(random_impact, c(1, 2), function(x) mean(x, na.rm = TRUE))
  random_sd <- apply(random_impact, c(1, 2), stats::sd, na.rm = TRUE)
  random_n <- apply(random_impact, c(1, 2), function(x) sum(is.finite(x)))
  sd_safe <- ifelse(!is.finite(random_sd) | random_sd == 0, 1, random_sd)
  z_delta <- (observed_impact - random_mean) / sd_safe
  z_delta[!is.finite(z_delta)] <- 0
  p_value <- 2 * stats::pnorm(-abs(z_delta))
  p_empirical <- .empirical_node_impact_p(observed_impact, random_impact, random_mean)

  dist_info <- .knockout_distance_info(g, remaining = remaining, targets = targets)
  impact <- .node_impact_long_table(
    targets = targets,
    remaining = remaining,
    metrics = metrics,
    before = base_mat,
    after = after_mat,
    delta = observed_delta,
    relative_change = observed_ri,
    impact_change = observed_impact,
    metric_sign = metric_sign,
    random_mean = random_mean,
    random_sd = random_sd,
    random_n = random_n,
    z_delta = z_delta,
    p_value = p_value,
    p_empirical = p_empirical,
    dist_info = dist_info
  )
  impact$P_adjust <- stats::p.adjust(impact$Pvalue, method = p_adjust_method)
  impact$P_adjust_metric <- ave(
    impact$Pvalue,
    impact$metric,
    FUN = function(x) stats::p.adjust(x, method = p_adjust_method)
  )
  impact$direction <- ifelse(
    impact$impact_change > 0,
    "disruption_like",
    ifelse(impact$impact_change < 0, "opposite", "unchanged")
  )
  impact <- impact[order(impact$P_adjust, -abs(impact$z_delta), impact$affected_protein, impact$metric), , drop = FALSE]
  rownames(impact) <- NULL

  protein_summary <- .node_impact_protein_summary(impact, alpha = alpha)

  out <- list(
    impact = impact,
    protein_summary = protein_summary,
    params = list(
      targets = targets,
      metrics = metrics,
      n_perm = n_perm,
      weight_attr = weight_attr,
      rewire_niter = rewire_niter,
      seed = seed,
      p_adjust_method = p_adjust_method,
      alpha = alpha
    )
  )
  class(out) <- c("tcm_ppi_knock_impact", "list")
  return(out)
}

#' @export
print.tcm_ppi_knock_impact <- function(x, ...) {
  cat("PPI node-level knockout impact\n")
  cat("Knocked targets:", paste(x$params$targets, collapse = ", "), "\n")
  cat("Metrics:", paste(x$params$metrics, collapse = ", "), "\n")
  cat("Random perturbations:", x$params$n_perm, "\n")
  cat("Affected proteins tested:", nrow(x$protein_summary), "\n")
  print(utils::head(x$protein_summary, 10), row.names = FALSE)
  invisible(x)
}


#' Compute four network-level topological metrics
#'
#' Returns ASPL (weighted), AD (average degree), DC (Freeman's degree
#' centralization), and CC (mean closeness centrality).
#' The closeness and degree calculations follow the same convention as
#' \code{\link{compute_nodeinfo}} in \file{PPI_utils.R}.
#'
#' @param g An igraph object.
#' @param weight_attr Edge attribute name for weights.
#' @return Named numeric vector of length 4.
#' @keywords internal
#' @noRd
.network_metrics <- function(g, weight_attr, metrics = c("ASPL", "AD", "DC", "CC")) {
  n <- igraph::vcount(g)
  metrics <- .validate_network_robustness_metrics(metrics)
  if (n == 0) {
    out <- setNames(rep(0, length(metrics)), metrics)
    out["ASPL"] <- if ("ASPL" %in% metrics) Inf else out["ASPL"]
    return(out)
  }

  w <- igraph::edge_attr(g, weight_attr)
  dist_weights <- if (!is.null(w) && length(w) > 0) {
    1 / pmax(w, .Machine$double.eps)
  } else {
    NULL
  }

  out <- lapply(metrics, function(metric) {
    switch(
      metric,
      ASPL = {
        aspl <- igraph::mean_distance(g, weights = dist_weights, directed = FALSE)
        if (is.nan(aspl) || is.na(aspl)) Inf else aspl
      },
      AD = mean(igraph::degree(g, mode = "all")),
      DC = igraph::centr_degree(g, mode = "all", normalized = TRUE)$centralization,
      CC = {
        cc_vals <- igraph::closeness(g, mode = "all", normalized = TRUE)
        mean(cc_vals, na.rm = TRUE)
      },
      density = igraph::edge_density(g, loops = FALSE),
      transitivity = {
        tr <- igraph::transitivity(g, type = "global")
        if (is.nan(tr) || is.na(tr)) 0 else tr
      },
      diameter = {
        dia <- suppressWarnings(igraph::diameter(g, directed = FALSE, weights = dist_weights, unconnected = TRUE))
        if (!is.finite(dia)) 0 else dia
      },
      components = igraph::components(g)$no,
      largest_component_fraction = {
        comp <- igraph::components(g)
        if (length(comp$csize) == 0L) 0 else max(comp$csize) / igraph::vcount(g)
      }
    )
  })
  stats::setNames(as.numeric(out), metrics)
}


#' Compute relative Robustness Index (RI)
#'
#' \eqn{RI = (M_{after} - M_{before}) / M_{before}}.
#' Falls back to the absolute difference when the baseline is zero.
#'
#' @param m_after Numeric vector of post-attack metrics.
#' @param m_before Numeric vector of baseline metrics.
#' @return Numeric vector of the same length.
#' @keywords internal
#' @noRd
.safe_ri <- function(m_after, m_before) {
  ri <- ifelse(m_before == 0,
               m_after - m_before,
               (m_after - m_before) / abs(m_before))
  ri[!is.finite(ri)] <- 0
  return(ri)
}

.validate_network_robustness_metrics <- function(metrics) {
  supported <- c(
    "ASPL", "AD", "DC", "CC",
    "density", "transitivity", "diameter",
    "components", "largest_component_fraction"
  )
  metrics_raw <- unique(as.character(metrics))
  idx <- match(tolower(metrics_raw), tolower(supported))
  missing <- metrics_raw[is.na(idx)]
  if (length(missing) > 0L) {
    stop(
      "Unsupported whole-network metric(s): ",
      paste(missing, collapse = ", "),
      ". Supported metrics are: ",
      paste(supported, collapse = ", "),
      call. = FALSE
    )
  }
  unique(supported[idx])
}

.validate_node_impact_metrics <- function(metrics) {
  supported <- c(
    "ASPL", "AD", "DC", "CC",
    "degree", "strength",
    "betweenness", "betweenness_w",
    "closeness", "closeness_w",
    "pagerank", "eigen_centrality",
    "coreness", "clustering_coef", "eccentricity",
    "radiality", "Stress", "MCC", "MNC", "DMNC", "BN", "EPC"
  )
  metrics_raw <- unique(as.character(metrics))
  idx <- match(tolower(metrics_raw), tolower(supported))
  missing <- metrics_raw[is.na(idx)]
  metrics <- supported[idx[!is.na(idx)]]
  metrics <- unique(metrics)
  if (length(missing) > 0L) {
    stop(
      "Unsupported node-level metric(s): ",
      paste(missing, collapse = ", "),
      ". Supported metrics are: ",
      paste(supported, collapse = ", "),
      call. = FALSE
    )
  }
  metrics
}

.node_impact_metric_sign <- function(metrics) {
  signs <- c(
    ASPL = 1,
    AD = -1,
    DC = -1,
    CC = -1,
    degree = -1,
    strength = -1,
    betweenness = -1,
    betweenness_w = -1,
    closeness = -1,
    closeness_w = -1,
    pagerank = -1,
    eigen_centrality = -1,
    coreness = -1,
    clustering_coef = -1,
    eccentricity = 1,
    radiality = -1,
    Stress = -1,
    MCC = -1,
    MNC = -1,
    DMNC = -1,
    BN = -1,
    EPC = -1
  )
  signs[metrics]
}

.node_metrics_selected <- function(g, metrics, weight_attr = "score") {
  node_names <- .get_node_names(g)
  out <- data.frame(name = node_names, stringsAsFactors = FALSE)
  edge_attrs <- igraph::edge_attr_names(g)
  w <- if (!is.null(weight_attr) && weight_attr %in% edge_attrs) {
    suppressWarnings(as.numeric(igraph::edge_attr(g, weight_attr)))
  } else {
    NULL
  }
  if (!is.null(w)) {
    w[!is.finite(w) | w <= 0] <- .Machine$double.eps
  }
  w_dist <- if (!is.null(w)) 1 / pmax(w, .Machine$double.eps) else NULL

  if ("MCC" %in% metrics) g <- compute_MCC(g)
  if ("MNC" %in% metrics) g <- compute_MNC(g)
  if ("DMNC" %in% metrics) g <- compute_DMNC(g)
  if ("BN" %in% metrics) g <- compute_BN(g)
  if ("radiality" %in% metrics) g <- compute_radiality(g)
  if ("Stress" %in% metrics) g <- compute_Stress(g)
  if ("EPC" %in% metrics) g <- compute_EPC(g, n_iter = 100)

  for (metric in metrics) {
    out[[metric]] <- switch(
      metric,
      ASPL = .node_aspl(g, weights = w_dist),
      AD = igraph::degree(g, mode = "all"),
      DC = .node_degree_centrality(g),
      CC = igraph::closeness(g, mode = "all", normalized = TRUE),
      degree = igraph::degree(g, mode = "all"),
      strength = if (!is.null(w)) igraph::strength(g, mode = "all", weights = w) else rep(NA_real_, length(node_names)),
      betweenness = igraph::betweenness(g, directed = FALSE, normalized = FALSE),
      betweenness_w = if (!is.null(w_dist)) {
        igraph::betweenness(g, directed = FALSE, weights = w_dist, normalized = FALSE)
      } else {
        rep(NA_real_, length(node_names))
      },
      closeness = igraph::closeness(g, mode = "all", normalized = FALSE),
      closeness_w = if (!is.null(w_dist)) {
        igraph::closeness(g, mode = "all", weights = w_dist, normalized = FALSE)
      } else {
        rep(NA_real_, length(node_names))
      },
      pagerank = igraph::page_rank(g, weights = w)$vector,
      eigen_centrality = igraph::eigen_centrality(g, weights = w)$vector,
      coreness = igraph::coreness(g, mode = "all"),
      clustering_coef = igraph::transitivity(g, type = "local", isolates = "zero"),
      eccentricity = igraph::eccentricity(g, mode = "all"),
      radiality = igraph::V(g)$radiality,
      Stress = igraph::V(g)$Stress,
      MCC = igraph::V(g)$MCC,
      MNC = igraph::V(g)$MNC,
      DMNC = igraph::V(g)$DMNC,
      BN = igraph::V(g)$BN,
      EPC = igraph::V(g)$EPC
    )
  }
  out
}

.node_aspl <- function(g, weights = NULL) {
  n <- igraph::vcount(g)
  if (n <= 1L) return(rep(0, n))
  dist_mat <- igraph::distances(g, v = igraph::V(g), to = igraph::V(g), weights = weights, mode = "all")
  out <- apply(dist_mat, 1, function(x) {
    x <- x[is.finite(x) & x > 0]
    if (length(x) == 0L) 0 else mean(x)
  })
  as.numeric(out)
}

.node_degree_centrality <- function(g) {
  n <- igraph::vcount(g)
  if (n <= 1L) return(rep(0, n))
  igraph::degree(g, mode = "all") / (n - 1)
}

.empirical_node_impact_p <- function(observed_ri, random_ri, random_mean) {
  out <- matrix(
    NA_real_,
    nrow = nrow(observed_ri),
    ncol = ncol(observed_ri),
    dimnames = dimnames(observed_ri)
  )
  for (i in seq_len(nrow(observed_ri))) {
    for (j in seq_len(ncol(observed_ri))) {
      random <- random_ri[i, j, ]
      random <- random[is.finite(random)]
      if (length(random) == 0L) {
        out[i, j] <- NA_real_
      } else {
        centered_obs <- abs(observed_ri[i, j] - random_mean[i, j])
        centered_random <- abs(random - random_mean[i, j])
        out[i, j] <- (sum(centered_random >= centered_obs) + 1) / (length(random) + 1)
      }
    }
  }
  out
}

.knockout_distance_info <- function(g, remaining, targets) {
  dist_mat <- igraph::distances(g, v = remaining, to = targets, weights = NA, mode = "all")
  min_dist <- apply(dist_mat, 1, function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else min(x)
  })
  data.frame(
    affected_protein = remaining,
    distance_to_knockout = as.numeric(min_dist[remaining]),
    is_direct_neighbor = as.logical(min_dist[remaining] == 1),
    stringsAsFactors = FALSE
  )
}

.node_impact_long_table <- function(targets,
                                    remaining,
                                    metrics,
                                    before,
                                    after,
                                    delta,
                                    relative_change,
                                    impact_change,
                                    metric_sign,
                                    random_mean,
                                    random_sd,
                                    random_n,
                                    z_delta,
                                    p_value,
                                    p_empirical,
                                    dist_info) {
  rows <- vector("list", length(remaining) * length(metrics))
  k <- 0L
  target_label <- paste(targets, collapse = ";")
  for (protein in remaining) {
    for (metric in metrics) {
      k <- k + 1L
      rows[[k]] <- data.frame(
        knocked_targets = target_label,
        affected_protein = protein,
        metric = metric,
        before = before[protein, metric],
        after = after[protein, metric],
        delta = delta[protein, metric],
        relative_change = relative_change[protein, metric],
        impact_change = impact_change[protein, metric],
        formula = .node_impact_formula(metric),
        metric_sign = metric_sign[metric],
        random_mean = random_mean[protein, metric],
        random_sd = random_sd[protein, metric],
        random_n = random_n[protein, metric],
        z_delta = z_delta[protein, metric],
        Pvalue = p_value[protein, metric],
        P_empirical = p_empirical[protein, metric],
        stringsAsFactors = FALSE
      )
    }
  }
  impact <- do.call(rbind, rows)
  merge(impact, dist_info, by = "affected_protein", all.x = TRUE, sort = FALSE)
}

.node_impact_formula <- function(metric) {
  sign <- .node_impact_metric_sign(metric)
  if (identical(as.numeric(sign), 1)) {
    "(after - before) / abs(before)"
  } else {
    "-(after - before) / abs(before)"
  }
}

.node_impact_protein_summary <- function(impact, alpha = 0.05) {
  split_impact <- split(impact, impact$affected_protein)
  out <- lapply(split_impact, function(x) {
    sig <- x$P_adjust < alpha
    sig[is.na(sig)] <- FALSE
    finite_p <- is.finite(x$P_adjust)
    min_p_adjust <- if (any(finite_p)) min(x$P_adjust[finite_p]) else NA_real_
    min_p_metric <- if (any(finite_p)) x$metric[which.min(x$P_adjust)] else NA_character_
    data.frame(
      knocked_targets = x$knocked_targets[[1]],
      affected_protein = x$affected_protein[[1]],
      distance_to_knockout = x$distance_to_knockout[[1]],
      is_direct_neighbor = x$is_direct_neighbor[[1]],
      n_metrics = nrow(x),
      n_significant_metrics = sum(sig),
      min_P_adjust = min_p_adjust,
      min_P_adjust_metric = min_p_metric,
      max_abs_z = max(abs(x$z_delta), na.rm = TRUE),
      mean_abs_z = mean(abs(x$z_delta), na.rm = TRUE),
      affected = any(sig),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out$min_P_adjust[!is.finite(out$min_P_adjust)] <- NA_real_
  out$max_abs_z[!is.finite(out$max_abs_z)] <- NA_real_
  out$mean_abs_z[!is.finite(out$mean_abs_z)] <- NA_real_
  out <- out[order(out$min_P_adjust, -out$max_abs_z, out$affected_protein), , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' Get node names, falling back to numeric IDs
#' @param g An igraph object.
#' @return Character vector of node names.
#' @keywords internal
#' @noRd
.get_node_names <- function(g) {
  nms <- igraph::V(g)$name
  if (is.null(nms)) {
    nms <- as.character(seq_len(igraph::vcount(g)))
    igraph::V(g)$name <- nms
  }
  return(nms)
}
