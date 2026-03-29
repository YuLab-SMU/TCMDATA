test_that("run_MCL works on Zachary graph", {
  skip_if_not_installed("igraph")

  g <- igraph::make_graph("Zachary")
  g_mcl <- run_MCL(g, inflation = 2)
  
  expect_true(igraph::is_igraph(g_mcl))
  expect_true("MCL_cluster" %in% igraph::vertex_attr_names(g_mcl))
  expect_type(igraph::V(g_mcl)$MCL_cluster, "integer")
})

test_that("run_mcode works on demo_ppi", {
  skip_if_not_installed("igraph")
  
  # demo_ppi should be available as it is lazy loaded or in data/
  # We might need to load it explicitly if testing without full package load
  if (!exists("demo_ppi")) {
    data("demo_ppi", package = "TCMDATA", envir = environment())
  }
  
  # If data loading fails (e.g. package not installed yet), we mock a graph
  if (!exists("demo_ppi")) {
    g <- igraph::make_ring(10)
    # Make it have some structure
    g <- igraph::add_edges(g, c(1,3, 1,4, 1,5, 2,4, 2,5))
  } else {
    g <- demo_ppi
  }
  
  res <- run_mcode(g, max_depth = 10)
  
  expect_true(igraph::is_igraph(res))
  expect_true("mcode_cluster" %in% igraph::vertex_attr_names(res))
  expect_true("mcode_score" %in% igraph::vertex_attr_names(res))
})
