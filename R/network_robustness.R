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
#'   closeness
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

  n_metrics    <- 4L
  metric_names <- c("ASPL", "AD", "DC", "CC")

  ## Stage 1: Baseline metrics
  m_base <- .network_metrics(g, weight_attr)

  ## Stage 2: Real attack — remove targets, compute RI (relative change)
  g_ko   <- igraph::delete_vertices(g, targets)
  m_ko   <- .network_metrics(g_ko, weight_attr)
  raw_ri <- .safe_ri(m_ko, m_base)

  ## Stage 3: Permutation null distribution
  rand_ri <- vapply(seq_len(n_perm), function(i) {
    g_rand <- igraph::rewire(g, igraph::keeping_degseq(niter = n_edges * rewire_niter))
    g_rand <- igraph::set_edge_attr(g_rand, weight_attr,
                                    value = sample(weights))
    m_rand_base  <- .network_metrics(g_rand, weight_attr)
    rand_targets <- sample(node_names, n_targets)
    g_rand_ko    <- igraph::delete_vertices(g_rand, rand_targets)
    m_rand_ko    <- .network_metrics(g_rand_ko, weight_attr)
    .safe_ri(m_rand_ko, m_rand_base)
  }, numeric(n_metrics))

  mu_rand <- rowMeans(rand_ri)
  sd_rand <- apply(rand_ri, 1, stats::sd)

  ## Stage 4: Z-score normalisation + p-value (normal approximation)
  sd_safe       <- ifelse(sd_rand == 0, 1, sd_rand)
  normalized_ri <- (raw_ri - mu_rand) / sd_safe

  # Two-sided p-value via Z-score (paper's 95% CI approach)
  p_zscore <- 2 * stats::pnorm(-abs(normalized_ri))

  ## Total Score = Z_ASPL - Z_AD - Z_DC - Z_CC
  total_score <- normalized_ri[1] - normalized_ri[2] -
    normalized_ri[3] - normalized_ri[4]

  ## Total Score p-value from permutation null distribution
  ## Combine per-permutation RIs using the same sign convention
  null_combined <- rand_ri[1, ] - rand_ri[2, ] - rand_ri[3, ] - rand_ri[4, ]
  real_combined <- raw_ri[1]   - raw_ri[2]   - raw_ri[3]   - raw_ri[4]
  mu_comb  <- mean(null_combined)
  sd_comb  <- stats::sd(null_combined)
  sd_comb  <- ifelse(sd_comb == 0, 1, sd_comb)
  z_total  <- (real_combined - mu_comb) / sd_comb
  p_total  <- 2 * stats::pnorm(-abs(z_total))

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
.network_metrics <- function(g, weight_attr) {
  n <- igraph::vcount(g)
  if (n == 0) return(c(ASPL = Inf, AD = 0, DC = 0, CC = 0))

  w <- igraph::edge_attr(g, weight_attr)
  dist_weights <- if (!is.null(w) && length(w) > 0) {
    1 / pmax(w, .Machine$double.eps)
  } else {
    NULL
  }

  # 1. ASPL — weighted average shortest path length
  aspl <- igraph::mean_distance(g, weights = dist_weights, directed = FALSE)
  if (is.nan(aspl) || is.na(aspl)) aspl <- Inf

  # 2. AD — average degree  (= 2|E| / |V|)
  ad <- mean(igraph::degree(g, mode = "all"))

  # 3. DC — Freeman's degree centralization (normalised to [0, 1])
  dc <- igraph::centr_degree(g, mode = "all", normalized = TRUE)$centralization

  # 4. CC — mean normalised closeness centrality
  #    (same convention as compute_nodeinfo: closeness(g, normalized = TRUE))
  cc_vals <- igraph::closeness(g, mode = "all", normalized = TRUE)
  cc <- mean(cc_vals, na.rm = TRUE)

  res <- c(ASPL = aspl, AD = ad, DC = dc, CC = cc)
  return(res)
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
