test_that("search_herb works", {
  # Skip if tcm_data is not available
  skip_if_not(exists("tcm_data", where = asNamespace("TCMDATA")) || exists("tcm_data"))
  
  # Test with a known herb (if possible) or just structure
  # Assuming "Ginseng" or similar exists in English or "灵芝" in Chinese
  # But safer to just test error handling if we are not sure about data content,
  # or check if we can mock it.
  # For now, let's try a simple case if we can.
  
  # If we cannot guarantee data presence, we can at least test argument validation
  expect_error(search_herb("Something", type = "InvalidType"))
})

test_that("search_herb handles non-existent herbs", {
  # This relies on tcm_data being present to trigger the warning logic
  skip_if_not(exists("tcm_data", where = asNamespace("TCMDATA")) || exists("tcm_data"))
  
  expect_warning(search_herb("NonExistentHerbXYZZ", "Herb_en_name"), "doesn't/don't exist")
})

test_that("search_herb uses the same result schema for each name type", {
  data <- get("tcm_data", envir = asNamespace("TCMDATA"))
  for (type in c("Herb_cn_name", "Herb_pinyin_name", "Herb_en_name")) {
    herb <- data[[type]][which(!is.na(data[[type]]))[1L]]
    result <- search_herb(herb, type)
    expect_named(result, c("herb", "molecule", "target"))
  }
})

test_that("search_target handles non-existent genes", {
   skip_if_not(exists("tcm_data", where = asNamespace("TCMDATA")) || exists("tcm_data"))
   
   expect_warning(search_target("NonExistentGeneXYZ"), "doesn't/don't exist")
})

test_that("legacy DisGeNET interfaces report that the data are unavailable", {
  expect_error(
    search_disease("sepsis"),
    "commercial resource.*licensing and copyright restrictions"
  )
  expect_error(
    search_gene_disease("IL6"),
    "commercial resource.*licensing and copyright restrictions"
  )
})
