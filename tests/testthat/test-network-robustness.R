test_that("ppi_knock_impact returns node-level perturbation table", {
  data(demo_ppi)
  target <- igraph::V(demo_ppi)$name[[1]]

  res <- ppi_knock_impact(
    demo_ppi,
    targets = target,
    n_perm = 5,
    seed = 1
  )

  expect_s3_class(res, "tcm_ppi_knock_impact")
  expect_true(all(c("impact", "protein_summary", "params") %in% names(res)))
  expect_equal(res$params$targets, target)
  expect_true(all(c(
    "knocked_targets", "affected_protein", "metric",
    "before", "after", "delta", "relative_change",
    "impact_change", "formula", "metric_sign",
    "random_mean", "random_sd", "random_n",
    "z_delta", "Pvalue", "P_empirical", "P_adjust",
    "distance_to_knockout", "is_direct_neighbor"
  ) %in% names(res$impact)))
  expect_false(target %in% res$impact$affected_protein)
  expect_equal(sort(unique(res$impact$metric)), c("AD", "ASPL", "CC", "DC"))
  expect_equal(nrow(res$protein_summary), igraph::vcount(demo_ppi) - 1L)
})

test_that("ppi_knock_impact supports multiple knockout targets", {
  data(demo_ppi)
  targets <- igraph::V(demo_ppi)$name[seq_len(2)]

  res <- ppi_knock_impact(
    demo_ppi,
    targets = targets,
    metrics = "AD",
    n_perm = 5,
    seed = 2
  )

  expect_s3_class(res, "tcm_ppi_knock_impact")
  expect_equal(sort(res$params$targets), sort(targets))
  expect_false(any(targets %in% res$impact$affected_protein))
  expect_equal(nrow(res$protein_summary), igraph::vcount(demo_ppi) - length(targets))
})

test_that("ppi_knock supports exploratory whole-network metrics", {
  data(demo_ppi)
  target <- igraph::V(demo_ppi)$name[[1]]

  res <- ppi_knock(
    demo_ppi,
    targets = target,
    metrics = c("ASPL", "AD", "density", "components"),
    n_perm = 5,
    seed = 3
  )

  expect_true(all(c("ASPL", "AD", "density", "components") %in% res$Summary$Metric))
  expect_true(is.na(res$Total_Score))
  expect_true(is.na(res$Total_Pvalue))
})

test_that("ppi_knock_impact supports extra PPI node metrics", {
  data(demo_ppi)
  target <- igraph::V(demo_ppi)$name[[1]]

  res <- ppi_knock_impact(
    demo_ppi,
    targets = target,
    metrics = c("ASPL", "pagerank", "MCC"),
    n_perm = 5,
    seed = 4
  )

  expect_true(all(c("ASPL", "pagerank", "MCC") %in% res$impact$metric))
  expect_true(all(res$impact$formula[res$impact$metric == "ASPL"] == "(after - before) / abs(before)"))
  expect_true(all(res$impact$formula[res$impact$metric == "pagerank"] == "-(after - before) / abs(before)"))
})

test_that("ppi_knock_impact validates requested metrics", {
  data(demo_ppi)
  target <- igraph::V(demo_ppi)$name[[1]]

  expect_error(
    ppi_knock_impact(demo_ppi, target, metrics = "not_a_metric", n_perm = 5),
    "Unsupported node-level metric"
  )
})
