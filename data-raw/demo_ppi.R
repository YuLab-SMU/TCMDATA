#' add a demo igraph object for PPI analysis
dir <- "/Users/hinna/Desktop/yulab/casestudy/demo_ppi.rds"
demo_ppi <- readRDS(dir)
usethis::use_data(demo_ppi, overwrite = TRUE)
