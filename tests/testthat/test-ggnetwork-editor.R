test_that("ggppi_network accepts fixed coordinate matrices", {
  graph <- igraph::make_ring(5)
  coordinates <- cbind(
    x = seq(-1, 1, length.out = 5),
    y = c(0, 1, 0, -1, 0)
  )

  plot <- ggppi_network(
    graph,
    layout = coordinates,
    show_text = FALSE
  )

  expect_equal(plot$data$x, unname(coordinates[, "x"]))
  expect_equal(plot$data$y, unname(coordinates[, "y"]))
  expect_silent(ggplot2::ggplot_build(plot))
})

test_that("ggppi_network reorders named layout coordinates by vertex name", {
  graph <- igraph::make_ring(4)
  graph <- igraph::set_vertex_attr(graph, "name", value = LETTERS[1:4])
  coordinates <- data.frame(
    name = c("D", "B", "A", "C"),
    x = c(4, 2, 1, 3),
    y = c(40, 20, 10, 30)
  )

  plot <- ggppi_network(
    graph,
    layout = coordinates,
    show_text = FALSE
  )

  expect_equal(plot$data$x, 1:4)
  expect_equal(plot$data$y, seq(10, 40, by = 10))
})

test_that("saved visNetwork positions become a reusable layout", {
  positions <- list(
    C = list(x = 30, y = 3),
    A = list(x = 10, y = 1),
    B = list(x = 20, y = 2)
  )

  coordinates <- .ggnetwork_positions_to_layout(
    positions,
    vertex_names = c("A", "B", "C")
  )

  expect_s3_class(coordinates, "tcm_network_layout")
  expect_equal(coordinates$name, c("A", "B", "C"))
  expect_equal(coordinates$x, c(10, 20, 30))
  expect_equal(coordinates$y, c(-1, -2, -3))
})

test_that("layout editor widget is draggable and uses short layout aliases", {
  skip_if_not_installed("visNetwork")
  graph <- igraph::make_ring(5)
  graph <- igraph::set_vertex_attr(graph, "name", value = LETTERS[1:5])
  graph <- igraph::set_edge_attr(
    graph,
    "score",
    value = seq(0.4, 0.8, length.out = igraph::ecount(graph))
  )

  widget <- .ggnetwork_editor_widget(
    graph,
    layout = "fr",
    physics = FALSE,
    seed = 1
  )

  expect_s3_class(widget, "visNetwork")
  expect_equal(.ggnetwork_vis_layout_name("fr"), "layout_with_fr")
  expect_true(all(c("x", "y") %in% names(widget$x$nodes)))
  expect_true("width" %in% names(widget$x$edges))
})

test_that("layout editor starts tripartite data in concentric rings", {
  skip_if_not_installed("visNetwork")
  data <- data.frame(
    herb = c("H1", "H1", "H2", "H2"),
    molecule = c("C1", "C2", "C1", "C3"),
    target = c("T1", "T2", "T1", "T3")
  )
  graph <- .ggtcm_as_graph(
    data,
    herb_col = "herb",
    compound_col = "molecule",
    target_col = "target"
  )

  widget <- .ggnetwork_editor_widget(
    graph,
    layout = "concentric",
    physics = FALSE,
    seed = 1
  )
  radius <- sqrt(widget$x$nodes$x^2 + widget$x$nodes$y^2) / 240
  type <- sub("::.*$", "", widget$x$nodes$id)

  expect_s3_class(widget, "visNetwork")
  expect_equal(unique(round(radius[type == "Herb"], 2)), 0.22)
  expect_equal(unique(round(radius[type == "Molecule"], 2)), 1.05)
  expect_equal(unique(round(radius[type == "Target"], 2)), 2)
  expect_true(all(widget$x$nodes$shape == "dot"))
  expect_setequal(
    unique(widget$x$nodes$color),
    c("#D94B2B", "#F07A3E", "#F5B06A")
  )
})

test_that("layout editor validates inputs before launching the gadget", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("visNetwork")
  expect_error(edit_ggnetwork_layout(data.frame()), "igraph object")
  expect_error(
    edit_ggnetwork_layout(igraph::make_ring(3), layout = NA_character_),
    "non-empty character"
  )
})
