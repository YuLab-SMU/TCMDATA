#' Plot a PPI network
#'
#' @param graph An `igraph` object.
#' @param metric,size_metric Numeric vertex attributes or supported centrality
#'   names mapped to node colour and size.
#' @param cluster Optional vertex attribute mapped to node border colour.
#' @param candidate Optional vertex attribute mapped to node shape (up to five
#'   categories).
#' @param score Optional numeric edge attribute mapped to edge width.
#' @param layout An igraph layout name or function, a coordinate matrix, or an
#'   object returned by [edit_ggnetwork_layout()].
#' @param show_text Whether to label high-scoring nodes.
#' @param label_top Number of nodes to label; use `Inf` for all nodes.
#' @param seed Random seed for layouts and labels; `NULL` leaves RNG unchanged.
#' @param size_range Node-size range.
#' @param colour_low,colour_high Node fill gradient.
#' @param edge_color,edge_alpha,edge_width Constant edge aesthetics.
#' @param edge_width_range Edge-width range when `score` is used.
#' @param base_size,base_family Base text size and family.
#' @param title Optional title.
#' @param ... Arguments passed to the igraph layout.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' data(demo_ppi)
#' p <- ggppi_network(
#'   demo_ppi,
#'   metric = "degree",
#'   size_metric = "betweenness",
#'   label_top = 8
#' )
#' p
#'
#' @importFrom ggtangle geom_edge
#' @importFrom ggplot2 aes coord_equal element_text geom_point ggplot guide_legend guides labs margin scale_colour_manual scale_fill_gradient scale_linewidth_continuous scale_shape_manual scale_size_continuous theme
#' @importFrom ggrepel geom_text_repel
#' @importFrom igraph betweenness closeness coreness degree eccentricity edge_attr edge_attr_names eigen_centrality is_directed is_igraph page_rank set_edge_attr set_vertex_attr strength transitivity vcount vertex_attr vertex_attr_names
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @export
ggppi_network <- function(graph,
                      metric = "degree",
                      size_metric = metric,
                      cluster = NULL,
                      candidate = NULL,
                      score = NULL,
                      layout = "nicely",
                      show_text = TRUE,
                      label_top = 10L,
                      seed = 42L,
                      size_range = c(2.5, 8),
                      colour_low = "#FEE8C8",
                      colour_high = "#D7301F",
                      edge_color = "#B8B8B8",
                      edge_alpha = 0.55,
                      edge_width = 0.35,
                      edge_width_range = c(0.15, 1.2),
                      base_size = 7,
                      base_family = "Arial",
                      title = NULL,
                      ...) {
  base_family <- .resolve_tcm_font_family(base_family)
  if (!is_igraph(graph)) {
    stop("'graph' must be an igraph object.", call. = FALSE)
  }
  if (vcount(graph) == 0L) {
    stop("'graph' must contain at least one node.", call. = FALSE)
  }
  if (!is.character(metric) || length(metric) != 1L || is.na(metric) ||
      !nzchar(metric)) {
    stop("'metric' must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(size_metric) || length(size_metric) != 1L ||
      is.na(size_metric) || !nzchar(size_metric)) {
    stop(
      "'size_metric' must be a single non-empty character string.",
      call. = FALSE
    )
  }
  .ggnetwork_validate_optional_name(cluster, "cluster")
  .ggnetwork_validate_optional_name(candidate, "candidate")
  .ggnetwork_validate_optional_name(score, "score")
  if (!is.logical(show_text) || length(show_text) != 1L || is.na(show_text)) {
    stop("'show_text' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(label_top) || length(label_top) != 1L || is.na(label_top) ||
      label_top < 0) {
    stop("'label_top' must be a non-negative number.", call. = FALSE)
  }
  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed))) {
    stop("'seed' must be NULL or one finite number.", call. = FALSE)
  }
  if (!is.numeric(size_range) || length(size_range) != 2L ||
      any(!is.finite(size_range)) || any(size_range <= 0) ||
      size_range[[1]] > size_range[[2]]) {
    stop("'size_range' must contain two positive, increasing values.", call. = FALSE)
  }
  if (!is.character(edge_color) || length(edge_color) != 1L ||
      is.na(edge_color) || !nzchar(edge_color)) {
    stop("'edge_color' must be a single non-empty colour value.", call. = FALSE)
  }
  if (!is.numeric(edge_alpha) || length(edge_alpha) != 1L ||
      !is.finite(edge_alpha) || edge_alpha < 0 || edge_alpha > 1) {
    stop("'edge_alpha' must be between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(edge_width) || length(edge_width) != 1L ||
      !is.finite(edge_width) || edge_width < 0) {
    stop("'edge_width' must be a non-negative number.", call. = FALSE)
  }
  if (!is.numeric(edge_width_range) || length(edge_width_range) != 2L ||
      any(!is.finite(edge_width_range)) || any(edge_width_range < 0) ||
      edge_width_range[[1]] > edge_width_range[[2]]) {
    stop(
      "'edge_width_range' must contain two non-negative, increasing values.",
      call. = FALSE
    )
  }
  layout <- .ggnetwork_resolve_layout(graph, layout)

  colour_values <- .ggnetwork_clean_metric(
    .ggnetwork_metric(graph, metric),
    metric
  )
  size_values <- if (identical(size_metric, metric)) {
    colour_values
  } else {
    .ggnetwork_clean_metric(
      .ggnetwork_metric(graph, size_metric),
      size_metric
    )
  }
  border_values <- .ggnetwork_vertex_group(
    graph,
    cluster,
    "cluster"
  )
  shape_values <- .ggnetwork_vertex_group(
    graph,
    candidate,
    "candidate"
  )
  if (!is.null(shape_values) && length(unique(shape_values)) > 5L) {
    stop("'candidate' supports at most five categories.", call. = FALSE)
  }
  edge_width_values <- if (is.null(score)) {
    NULL
  } else {
    .ggnetwork_clean_metric(
      .ggnetwork_edge_metric(graph, score),
      score
    )
  }

  graph <- set_vertex_attr(
    graph,
    name = ".tcm_colour_metric",
    value = as.numeric(colour_values)
  )
  graph <- set_vertex_attr(
    graph,
    name = ".tcm_size_metric",
    value = as.numeric(size_values)
  )
  if (!is.null(border_values)) {
    graph <- set_vertex_attr(
      graph,
      name = ".tcm_cluster",
      value = border_values
    )
  }
  if (!is.null(shape_values)) {
    graph <- set_vertex_attr(
      graph,
      name = ".tcm_candidate",
      value = shape_values
    )
  }
  if (!is.null(edge_width_values)) {
    graph <- set_edge_attr(
      graph,
      name = ".tcm_score",
      value = as.numeric(edge_width_values)
    )
  }
  metric_title <- .ggnetwork_metric_title(metric)
  size_metric_title <- .ggnetwork_metric_title(size_metric)
  cluster_title <- if (is.null(cluster)) {
    NULL
  } else {
    .ggnetwork_metric_title(cluster)
  }
  candidate_title <- if (is.null(candidate)) {
    NULL
  } else {
    .ggnetwork_metric_title(candidate)
  }
  score_title <- if (is.null(score)) {
    NULL
  } else {
    .ggnetwork_metric_title(score)
  }

  edge_layer <- if (is.null(score)) {
    geom_edge(
      colour = edge_color,
      alpha = edge_alpha,
      linewidth = edge_width,
      lineend = "round"
    )
  } else {
    geom_edge(
      aes(linewidth = .data[[".tcm_score"]]),
      colour = edge_color,
      alpha = edge_alpha,
      lineend = "round"
    )
  }

  node_layer <- .ggnetwork_node_layer(
    cluster = cluster,
    candidate = candidate
  )

  p <- .ggnetwork_with_seed(
    seed,
    ggplot(graph, layout = layout, ...)
  ) +
    edge_layer +
    node_layer +
    scale_fill_gradient(
      low = colour_low,
      high = colour_high,
      name = metric_title
    ) +
    scale_size_continuous(
      range = size_range,
      name = size_metric_title
    ) +
    guides(
      size = if (identical(size_metric, metric)) {
        "none"
      } else {
        guide_legend(
          override.aes = list(
            fill = colour_high,
            colour = "#4D4D4D",
            shape = 21,
            alpha = 1
          )
        )
      }
    ) +
    coord_equal(clip = "off") +
    labs(title = title) +
    .theme_tcm_void(base_size = base_size, base_family = base_family) +
    theme(
      legend.position = "right",
      legend.title = element_text(face = "plain"),
      plot.margin = margin(5.5, 8, 5.5, 5.5)
    )

  if (!is.null(cluster)) {
    border_levels <- sort(unique(border_values))
    border_palette <- setNames(
      .pal_tcm_methods(length(border_levels)),
      border_levels
    )
    p <- p + scale_colour_manual(
      values = border_palette,
      name = cluster_title,
      guide = guide_legend(
        override.aes = list(fill = "white", size = 4, shape = 21)
      )
    )
  }
  if (!is.null(candidate)) {
    shape_levels <- sort(unique(shape_values))
    shape_palette <- setNames(
      c(21, 24, 22, 23, 25)[seq_along(shape_levels)],
      shape_levels
    )
    p <- p + scale_shape_manual(
      values = shape_palette,
      name = candidate_title,
      guide = guide_legend(
        override.aes = list(
          fill = colour_high,
          colour = "#4D4D4D",
          size = 4
        )
      )
    )
  }
  if (!is.null(score)) {
    p <- p + scale_linewidth_continuous(
      range = edge_width_range,
      name = score_title
    )
  }

  n_labels <- if (is.infinite(label_top)) {
    nrow(p$data)
  } else {
    min(as.integer(label_top), nrow(p$data))
  }
  if (isTRUE(show_text) && n_labels > 0L) {
    label_index <- order(
      p$data[[".tcm_colour_metric"]],
      decreasing = TRUE,
      na.last = NA
    )[seq_len(n_labels)]
    label_data <- p$data
    label_data[[".tcm_label"]] <- ""
    label_data[[".tcm_label"]][label_index] <-
      label_data[["label"]][label_index]
    label_data[[".tcm_point_size"]] <- .ggnetwork_rescale_size(
      label_data[[".tcm_size_metric"]],
      size_range
    )
    p <- p + geom_text_repel(
      data = label_data,
      aes(
        label = .data[[".tcm_label"]],
        point.size = .data[[".tcm_point_size"]]
      ),
      size = max(2, base_size / 3),
      family = base_family,
      fontface = "plain",
      colour = "#272727",
      box.padding = 0.4,
      point.padding = 0.25,
      min.segment.length = 0,
      segment.colour = "#9A9A9A",
      segment.size = 0.25,
      force = 2,
      max.time = 2,
      max.overlaps = Inf,
      seed = seed,
      show.legend = FALSE
    )
  }

  p
}

.ggnetwork_validate_optional_name <- function(value, argument) {
  if (is.null(value)) {
    return(invisible(NULL))
  }
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value)) {
    stop(
      sprintf("'%s' must be NULL or one non-empty character string.", argument),
      call. = FALSE
    )
  }

  invisible(NULL)
}

.ggnetwork_resolve_layout <- function(graph, layout) {
  if (!is.matrix(layout) && !is.data.frame(layout)) {
    return(layout)
  }

  if (is.data.frame(layout)) {
    if (!all(c("x", "y") %in% names(layout))) {
      stop("A layout data frame must contain 'x' and 'y' columns.", call. = FALSE)
    }
    coordinates <- layout[, c("x", "y"), drop = FALSE]
    layout_names <- if ("name" %in% names(layout)) {
      as.character(layout$name)
    } else {
      rownames(layout)
    }
  } else {
    if (ncol(layout) < 2L) {
      stop("A layout matrix must contain at least two columns.", call. = FALSE)
    }
    coordinates <- layout[, 1:2, drop = FALSE]
    layout_names <- rownames(layout)
  }

  if (nrow(coordinates) != vcount(graph)) {
    stop("The layout must contain one row for each graph node.", call. = FALSE)
  }
  graph_names <- if ("name" %in% vertex_attr_names(graph)) {
    as.character(vertex_attr(graph, "name"))
  } else {
    NULL
  }
  if (!is.null(graph_names) && !is.null(layout_names) &&
      length(layout_names) == nrow(coordinates) &&
      !anyNA(layout_names) && !anyDuplicated(layout_names)) {
    match_index <- match(graph_names, layout_names)
    if (anyNA(match_index)) {
      stop("The layout does not contain every graph node name.", call. = FALSE)
    }
    coordinates <- coordinates[match_index, , drop = FALSE]
  }

  coordinates <- as.matrix(coordinates)
  storage.mode(coordinates) <- "double"
  if (any(!is.finite(coordinates))) {
    stop("The layout contains non-finite coordinates.", call. = FALSE)
  }
  fixed_coordinates <- unname(coordinates)
  function(graph, ...) fixed_coordinates
}

.ggnetwork_vertex_group <- function(graph, metric, argument) {
  if (is.null(metric)) {
    return(NULL)
  }
  if (!metric %in% vertex_attr_names(graph)) {
    stop(
      sprintf(
        "Vertex attribute '%s' requested by '%s' was not found.",
        metric,
        argument
      ),
      call. = FALSE
    )
  }

  values <- vertex_attr(graph, metric)
  if (is.list(values) || length(values) != vcount(graph)) {
    stop(
      sprintf("Vertex attribute '%s' must be an atomic vector.", metric),
      call. = FALSE
    )
  }
  values <- as.character(values)
  values[is.na(values) | !nzchar(values)] <- "Unassigned"
  values
}

.ggnetwork_edge_metric <- function(graph, metric) {
  if (!metric %in% edge_attr_names(graph)) {
    stop(
      sprintf("Edge attribute '%s' was not found.", metric),
      call. = FALSE
    )
  }

  values <- edge_attr(graph, metric)
  if (!is.numeric(values)) {
    stop(
      sprintf("Edge attribute '%s' must be numeric.", metric),
      call. = FALSE
    )
  }
  values
}

.ggnetwork_node_layer <- function(cluster, candidate) {
  if (!is.null(cluster) && !is.null(candidate)) {
    return(geom_point(
      aes(
        fill = .data[[".tcm_colour_metric"]],
        size = .data[[".tcm_size_metric"]],
        colour = .data[[".tcm_cluster"]],
        shape = .data[[".tcm_candidate"]]
      ),
      stroke = 0.65,
      alpha = 0.96
    ))
  }
  if (!is.null(cluster)) {
    return(geom_point(
      aes(
        fill = .data[[".tcm_colour_metric"]],
        size = .data[[".tcm_size_metric"]],
        colour = .data[[".tcm_cluster"]]
      ),
      shape = 21,
      stroke = 0.65,
      alpha = 0.96
    ))
  }
  if (!is.null(candidate)) {
    return(geom_point(
      aes(
        fill = .data[[".tcm_colour_metric"]],
        size = .data[[".tcm_size_metric"]],
        shape = .data[[".tcm_candidate"]]
      ),
      colour = "white",
      stroke = 0.3,
      alpha = 0.96
    ))
  }

  geom_point(
    aes(
      fill = .data[[".tcm_colour_metric"]],
      size = .data[[".tcm_size_metric"]]
    ),
    shape = 21,
    colour = "white",
    stroke = 0.3,
    alpha = 0.96
  )
}

.ggnetwork_metric <- function(graph, metric) {
  vertex_attributes <- vertex_attr_names(graph)
  if (metric %in% vertex_attributes) {
    values <- vertex_attr(graph, metric)
    if (!is.numeric(values)) {
      stop(
        sprintf("Vertex attribute '%s' must be numeric.", metric),
        call. = FALSE
      )
    }
    return(values)
  }

  edge_weights <- if ("weight" %in% edge_attr_names(graph)) {
    edge_attr(graph, "weight")
  } else {
    NULL
  }

  values <- switch(
    metric,
    degree = degree(graph, mode = "all"),
    strength = strength(graph, mode = "all", weights = edge_weights),
    betweenness = betweenness(
      graph,
      directed = is_directed(graph),
      normalized = TRUE
    ),
    closeness = closeness(graph, mode = "all", normalized = TRUE),
    eigen_centrality = eigen_centrality(
      graph,
      directed = is_directed(graph),
      weights = edge_weights
    )$vector,
    pagerank = page_rank(
      graph,
      directed = is_directed(graph),
      weights = edge_weights
    )$vector,
    coreness = coreness(graph, mode = "all"),
    clustering_coef = transitivity(
      graph,
      type = "local",
      isolates = "zero"
    ),
    eccentricity = eccentricity(graph, mode = "all"),
    NULL
  )

  if (is.null(values)) {
    available <- sort(unique(c(
      vertex_attributes,
      "degree", "strength", "betweenness", "closeness",
      "eigen_centrality", "pagerank", "coreness",
      "clustering_coef", "eccentricity"
    )))
    stop(
      sprintf(
        "Metric '%s' was not found. Available metrics: %s.",
        metric,
        paste(available, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  values
}

.ggnetwork_clean_metric <- function(values, metric) {
  finite_metric <- is.finite(values)
  if (!any(finite_metric)) {
    stop(sprintf("Metric '%s' has no finite values.", metric), call. = FALSE)
  }
  if (any(!finite_metric)) {
    warning(
      sprintf(
        "Metric '%s' contains non-finite values; they are shown at the minimum finite value.",
        metric
      ),
      call. = FALSE
    )
    values[!finite_metric] <- min(values[finite_metric])
  }

  values
}

.ggnetwork_metric_title <- function(metric) {
  label <- gsub("_", " ", metric, fixed = TRUE)
  paste0(toupper(substr(label, 1L, 1L)), substr(label, 2L, nchar(label)))
}

.ggnetwork_rescale_size <- function(values, size_range) {
  value_range <- range(values, finite = TRUE)
  if (diff(value_range) == 0) {
    return(rep(mean(size_range), length(values)))
  }

  size_range[[1]] +
    (values - value_range[[1]]) / diff(value_range) * diff(size_range)
}

.ggnetwork_with_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(force(expr))
  }

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(as.integer(seed))
  force(expr)
}

#' Edit a network layout
#'
#' @param graph An `igraph` object or a herb-compound-target data frame.
#' @param layout Initial layout; `"concentric"` and short igraph names such as
#'   `"fr"` are supported.
#' @param physics Whether to enable layout physics initially.
#' @param seed Initial-layout seed.
#' @param width,height Window dimensions in pixels.
#' @param viewer Optional Shiny viewer function.
#' @param ... Arguments passed to the igraph layout.
#'
#' @return A `tcm_network_layout` data frame, or `NULL` when cancelled.
#'
#' @examples
#' \dontrun{
#' data(demo_ppi)
#' ppi <- run_mcode(demo_ppi)
#' coordinates <- edit_ggnetwork_layout(ppi, layout = "fr")
#' p <- ggppi_network(
#'   ppi,
#'   layout = coordinates,
#'   cluster = "cluster",
#'   candidate = "candidate",
#'   score = "score"
#' )
#' ggplot2::ggsave("ppi_network.pdf", p, width = 7.2, height = 7.2)
#'
#' hct <- search_herb("lingzhi", "Herb_pinyin_name")
#' coordinates <- edit_ggnetwork_layout(hct, layout = "concentric")
#' ggtcm_network(hct, layout = coordinates, label_top = Inf)
#' }
#'
#' @export
edit_ggnetwork_layout <- function(graph,
                                  layout = "fr",
                                  physics = FALSE,
                                  seed = 42L,
                                  width = 1000,
                                  height = 750,
                                  viewer = NULL,
                                  ...) {
  for (package_name in c("shiny", "visNetwork")) {
    if (!requireNamespace(package_name, quietly = TRUE)) {
      stop(
        sprintf("Package '%s' is required for edit_ggnetwork_layout().", package_name),
        call. = FALSE
      )
    }
  }
  if (is.data.frame(graph)) {
    if (!all(c("herb", "molecule", "target") %in% names(graph))) {
      stop(
        "'graph' must be an igraph object or a herb-compound-target data frame.",
        call. = FALSE
      )
    }
    graph <- .ggtcm_as_graph(
      graph,
      herb_col = "herb",
      compound_col = "molecule",
      target_col = "target"
    )
  }
  if (!is_igraph(graph)) {
    stop(
      "'graph' must be an igraph object or a herb-compound-target data frame.",
      call. = FALSE
    )
  }
  if (vcount(graph) == 0L) {
    stop("'graph' must contain at least one node.", call. = FALSE)
  }
  if (!is.character(layout) || length(layout) != 1L || is.na(layout) ||
      !nzchar(layout)) {
    stop("'layout' must be one non-empty character string.", call. = FALSE)
  }
  if (!is.logical(physics) || length(physics) != 1L || is.na(physics)) {
    stop("'physics' must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
    stop("'seed' must be one finite number.", call. = FALSE)
  }
  if (!is.numeric(width) || length(width) != 1L || !is.finite(width) ||
      width <= 0 || !is.numeric(height) || length(height) != 1L ||
      !is.finite(height) || height <= 0) {
    stop("'width' and 'height' must be positive numbers.", call. = FALSE)
  }

  graph <- .ggnetwork_ensure_vertex_names(graph)
  vertex_names <- vertex_attr(graph, "name")
  widget <- .ggnetwork_editor_widget(
    graph = graph,
    layout = layout,
    physics = physics,
    seed = seed,
    ...
  )
  initial_positions <- widget$x$nodes[, c("id", "x", "y"), drop = FALSE]

  ui <- shiny::fluidPage(
    shiny::fluidRow(
      shiny::column(
        width = 9,
        visNetwork::visNetworkOutput("tcm_network_editor", height = "680px")
      ),
      shiny::column(
        width = 3,
        shiny::h4("Network layout editor"),
        shiny::helpText(
          "Drag nodes to arrange the network. Zoom and pan are enabled. ",
          "Save layout returns the coordinates to R."
        ),
        shiny::checkboxInput(
          "tcm_layout_physics",
          "Enable physics",
          value = physics
        ),
        shiny::actionButton("tcm_layout_fit", "Fit to window"),
        shiny::actionButton("tcm_layout_reset", "Reset layout"),
        shiny::actionButton("tcm_layout_save", "Save layout"),
        shiny::actionButton("cancel", "Cancel"),
        shiny::textOutput("tcm_layout_status")
      )
    )
  )

  server <- function(input, output, session) {
    output$tcm_network_editor <- visNetwork::renderVisNetwork(widget)
    output$tcm_layout_status <- shiny::renderText(
      "The saved layout can be passed to ggppi_network() or ggtcm_network()."
    )

    shiny::observeEvent(input$tcm_layout_physics, {
      visNetwork::visPhysics(
        visNetwork::visNetworkProxy("tcm_network_editor", session = session),
        enabled = isTRUE(input$tcm_layout_physics)
      )
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$tcm_layout_fit, {
      visNetwork::visFit(visNetwork::visNetworkProxy("tcm_network_editor", session = session))
    })

    shiny::observeEvent(input$tcm_layout_reset, {
      proxy <- visNetwork::visNetworkProxy("tcm_network_editor", session = session)
      visNetwork::visPhysics(proxy, enabled = FALSE)
      visNetwork::visUpdateNodes(proxy, initial_positions)
      visNetwork::visFit(proxy)
    })

    shiny::observeEvent(input$tcm_layout_save, {
      visNetwork::visGetPositions(
        visNetwork::visNetworkProxy("tcm_network_editor", session = session),
        input = "tcm_saved_positions"
      )
    })

    shiny::observeEvent(input$tcm_saved_positions, {
      coordinates <- .ggnetwork_positions_to_layout(
        input$tcm_saved_positions,
        vertex_names = vertex_names
      )
      attr(coordinates, "layout_type") <- layout
      shiny::stopApp(coordinates)
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$cancel, {
      shiny::stopApp(NULL)
    })
  }

  if (is.null(viewer)) {
    viewer <- shiny::dialogViewer(
      "TCMDATA network layout editor",
      width = as.integer(width),
      height = as.integer(height)
    )
  }

  shiny::runGadget(
    ui,
    server,
    viewer = viewer,
    stopOnCancel = FALSE
  )
}

.ggnetwork_editor_widget <- function(graph,
                                     layout,
                                     physics,
                                     seed,
                                     ...) {
  network_data <- visNetwork::toVisNetworkData(graph)
  is_concentric <- identical(layout, "concentric")
  if (is_concentric) {
    if (!"type" %in% vertex_attr_names(graph)) {
      stop(
        "The 'concentric' editor layout requires a 'type' vertex attribute.",
        call. = FALSE
      )
    }
    editor_graph <- .ggtcm_prepare_attributes(graph, node_size = "degree")
    coordinates <- .ggtcm_concentric_layout(
      editor_graph,
      ring_radii = c(Herb = 0.22, Molecule = 1.05, Target = 2)
    )
    node_names <- as.character(vertex_attr(editor_graph, "name"))
    node_index <- match(as.character(network_data$nodes$id), node_names)
    if (anyNA(node_index)) {
      stop("Could not match editor nodes to concentric coordinates.", call. = FALSE)
    }
    node_type <- vertex_attr(editor_graph, ".tcm_type")[node_index]
    node_palette <- c(
      Herb = "#D94B2B",
      Molecule = "#F07A3E",
      Target = "#F5B06A"
    )
    node_sizes <- .ggtcm_concentric_sizes(
      editor_graph,
      type_sizes = c(Herb = 11, Molecule = 6, Target = 4),
      node_size = "degree"
    )[node_index]
    network_data$nodes$x <- coordinates[node_index, 1L] * 240
    network_data$nodes$y <- -coordinates[node_index, 2L] * 240
    network_data$nodes$label <- vertex_attr(
      editor_graph,
      ".tcm_label"
    )[node_index]
    network_data$nodes$size <- node_sizes * 2.5
    network_data$nodes$shape <- "dot"
    network_data$nodes$color <- unname(node_palette[node_type])
  } else {
    node_degree <- degree(graph, mode = "all")
    network_data$nodes$size <- .ggnetwork_rescale_size(
      as.numeric(node_degree),
      c(12, 30)
    )
  }
  if ("candidate" %in% vertex_attr_names(graph)) {
    candidate <- as.logical(vertex_attr(graph, "candidate"))
    candidate[is.na(candidate)] <- FALSE
    network_data$nodes$shape <- ifelse(candidate, "triangle", "dot")
  }
  if ("score" %in% edge_attr_names(graph)) {
    score <- edge_attr(graph, "score")
    if (is.numeric(score) && any(is.finite(score))) {
      score[!is.finite(score)] <- min(score[is.finite(score)])
      network_data$edges$width <- .ggnetwork_rescale_size(
        score,
        c(0.5, 4)
      )
    }
  }

  widget <- visNetwork::visNetwork(
    network_data$nodes,
    network_data$edges,
    width = "100%",
    height = "680px",
    background = "#FFFFFF"
  )
  if (!is_concentric) {
    widget <- visNetwork::visIgraphLayout(
      widget,
      layout = .ggnetwork_vis_layout_name(layout),
      physics = physics,
      smooth = FALSE,
      randomSeed = as.integer(seed),
      ...
    )
  }
  widget <- visNetwork::visNodes(
    widget,
    color = list(
      background = "#E95C3A",
      border = "#B73522",
      highlight = list(background = "#D7301F", border = "#7F1D12")
    ),
    borderWidth = 1.5,
    font = list(color = "#272727", face = "Arial", size = 14)
  )
  widget <- visNetwork::visEdges(
    widget,
    color = list(color = "#B8B8B8", opacity = 0.55),
    smooth = FALSE
  )
  widget <- visNetwork::visInteraction(
    widget,
    dragNodes = TRUE,
    dragView = TRUE,
    hover = TRUE,
    keyboard = TRUE,
    navigationButtons = TRUE,
    zoomView = TRUE
  )
  widget <- visNetwork::visOptions(
    widget,
    highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
    nodesIdSelection = TRUE
  )
  visNetwork::visPhysics(widget, enabled = physics, stabilization = FALSE)
}

.ggnetwork_vis_layout_name <- function(layout) {
  aliases <- c(
    nicely = "layout_nicely",
    fr = "layout_with_fr",
    kk = "layout_with_kk",
    circle = "layout_in_circle",
    star = "layout_as_star",
    grid = "layout_on_grid",
    drl = "layout_with_drl",
    lgl = "layout_with_lgl",
    mds = "layout_with_mds",
    graphopt = "layout_with_graphopt"
  )
  if (layout %in% names(aliases)) {
    return(unname(aliases[[layout]]))
  }
  layout
}

.ggnetwork_positions_to_layout <- function(positions, vertex_names) {
  if (is.data.frame(positions)) {
    if (!all(c("id", "x", "y") %in% names(positions))) {
      stop("Saved positions must contain 'id', 'x', and 'y'.", call. = FALSE)
    }
    coordinates <- positions[, c("id", "x", "y"), drop = FALSE]
    names(coordinates)[[1L]] <- "name"
  } else if (is.list(positions) && length(positions) > 0L) {
    position_ids <- names(positions)
    if (is.null(position_ids) || any(!nzchar(position_ids))) {
      stop("Saved node positions do not contain node identifiers.", call. = FALSE)
    }
    coordinates <- data.frame(
      name = position_ids,
      x = vapply(positions, function(value) as.numeric(value$x), numeric(1)),
      y = vapply(positions, function(value) as.numeric(value$y), numeric(1)),
      stringsAsFactors = FALSE
    )
  } else {
    stop("No node positions were returned by the layout editor.", call. = FALSE)
  }

  match_index <- match(vertex_names, as.character(coordinates$name))
  if (anyNA(match_index)) {
    stop("The saved layout does not contain every graph node.", call. = FALSE)
  }
  coordinates <- coordinates[match_index, , drop = FALSE]
  coordinates$name <- as.character(coordinates$name)
  coordinates$x <- as.numeric(coordinates$x)
  # vis.js uses screen coordinates, where y increases downwards. Reverse y so
  # that the final ggplot orientation matches the editor.
  coordinates$y <- -as.numeric(coordinates$y)
  if (any(!is.finite(coordinates$x)) || any(!is.finite(coordinates$y))) {
    stop("The saved layout contains non-finite coordinates.", call. = FALSE)
  }
  rownames(coordinates) <- NULL
  class(coordinates) <- c("tcm_network_layout", "data.frame")
  coordinates
}

.ggnetwork_ensure_vertex_names <- function(graph) {
  vertex_names <- if ("name" %in% vertex_attr_names(graph)) {
    as.character(vertex_attr(graph, "name"))
  } else {
    character(0)
  }
  if (length(vertex_names) != vcount(graph) || anyNA(vertex_names) ||
      any(!nzchar(vertex_names)) || anyDuplicated(vertex_names)) {
    vertex_names <- as.character(seq_len(vcount(graph)))
    graph <- set_vertex_attr(graph, "name", value = vertex_names)
  }
  graph
}
