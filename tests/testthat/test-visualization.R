# Tests for visualization functions in TCMDATA package

test_that("gglollipop works with basic input", {
  # Create a simple mock enrichment result
  mock_enrich <- data.frame(
    ID = paste0("GO:", 1:10),
    Description = paste("Pathway", 1:10),
    GeneRatio = paste0(1:10, "/100"),
    BgRatio = paste0(5:14, "/200"),
    pvalue = runif(10, 0, 0.05),
    p.adjust = runif(10, 0, 0.05),
    qvalue = runif(10, 0, 0.05),
    geneID = paste0("Gene", 1:10),
    Count = 1:10,
    RichFactor = (1:10)/100
  )
  
  # Construct S4 object if class is available, otherwise skip
  if (requireNamespace("clusterProfiler", quietly = TRUE) && methods::isClass("enrichResult", where = asNamespace("clusterProfiler"))) {
     # Try to load class definition to ensure new() works correctly
     suppressPackageStartupMessages(requireNamespace("clusterProfiler"))
     
     mock_enrich <- new("enrichResult", 
                        result = mock_enrich, 
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.05, 
                        organism = "hsa", 
                        ontology = "BP",
                        gene = character(),
                        keytype = "ENTREZID")
     
     # Verify it is S4, otherwise skip (some environments might return S3 if setOldClass is involved or other issues)
     if (!isS4(mock_enrich)) {
       skip("Created enrichResult object is not S4")
     }
  } else {
     skip("enrichResult class not available")
  }
  
  # Test that function runs without error
  expect_silent(p <- gglollipop(mock_enrich, top_n = 5))
  expect_s3_class(p, "ggplot")
})

test_that("ggvenn_plot works with basic input", {
  # Create a data frame for ggvenn (1st col ID, others logical)
  all_genes <- unique(c(paste0("Gene", 1:20), paste0("Gene", 10:30), paste0("Gene", 15:35)))
  venn_df <- data.frame(Element = all_genes)
  venn_df$Set1 <- venn_df$Element %in% paste0("Gene", 1:20)
  venn_df$Set2 <- venn_df$Element %in% paste0("Gene", 10:30)
  venn_df$Set3 <- venn_df$Element %in% paste0("Gene", 15:35)
  
  # Test that function runs without error
  expect_silent(p <- ggvenn_plot(venn_df))
  expect_s3_class(p, "ggplot")
})

test_that("ggdock works with basic input", {
  # Create mock docking data (wide format)
  mock_docking <- matrix(rnorm(15, -7, 1), nrow = 5, ncol = 3)
  rownames(mock_docking) <- paste0("Target", 1:5)
  colnames(mock_docking) <- paste0("Mol", 1:3)
  mock_docking <- as.data.frame(mock_docking)
  
  # Test that function runs without error
  expect_silent(p <- ggdock(mock_docking))
  expect_s3_class(p, "ggplot")
  expect_silent(ggplot2::ggplot_build(ggdock(mock_docking, order = "median")))
})

test_that("TCM_sankey works with basic input", {
  # Create mock data for sankey plot
  mock_data <- data.frame(
    herb = rep(paste0("Herb", 1:3), each = 4),
    molecule = rep(paste0("Mol", 1:4), times = 3),
    target = rep(paste0("Target", 1:6), times = 2)
  )
  
  # Test that function runs without error
  expect_silent(p <- tcm_sankey(mock_data))
  expect_s3_class(p, "ggplot")
})

test_that("PlotNodeHeatmap requires ComplexHeatmap", {
  # Skip if ComplexHeatmap not installed
  testthat::skip_if_not_installed("ComplexHeatmap")
  
  # Create mock node metrics data
  mock_metrics <- data.frame(
    name = paste0("Node", 1:10),
    centrality = runif(10),
    betweenness = runif(10),
    degree = sample(1:10, 10)
  )
  
  # Test that function runs without error when ComplexHeatmap is available
  expect_silent(p <- plot_node_heatmap(mock_metrics))
  expect_s4_class(p, "Heatmap")
})

test_that("gocircle_plot works with basic input", {
  # Create mock data
  mock_data <- data.frame(
    ID = paste0("GO:", 1:10),
    Description = paste("Pathway", 1:10),
    GeneRatio = paste0(1:10, "/100"),
    BgRatio = paste0(1:10, "/100"),
    geneID = paste0("Gene", 1:10),
    ONTOLOGY = "BP",
    Count = 1:10,
    RichFactor = runif(10, 0, 1),
    p.adjust = runif(10, 0, 0.05)
  )
  
  # Test that function runs without error
  # gocircle_plot uses circlize and returns TRUE invisibly
  expect_true(gocircle_plot(mock_data, top = 5))
})

test_that("radar_plot works with basic input", {
  # Create mock data
  mock_data <- data.frame(
    Category = rep(LETTERS[1:5], each = 3),
    Group = rep(paste0("Group", 1:3), times = 5),
    Value = runif(15, 0, 1)
  )
  
  # Test that function runs without error
  expect_silent(p <- radar_plot(mock_data, category = "Category", value = "Value"))
  expect_s3_class(p, "ggplot")
})

test_that("visualization functions handle invalid input gracefully", {
  # Test gglollipop with non-enrichment input
  expect_error(gglollipop(data.frame()))
  
  # Test ggvenn_plot with empty list
  expect_error(ggvenn_plot(list()), "ggvenn_plot only supports 2 to 4 sets")
  
  # Test ggdock with missing columns
  expect_error(ggdock(data.frame()))
})
