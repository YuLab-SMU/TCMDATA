#' add a demo data from GSE142025 for documentation
dir <- "/Users/hinna/Desktop/yulab/casestudy/deg_all.rda"
load(dir)
deg_earlydn <- deg_all[[1]]
usethis::use_data(deg_earlydn, overwrite = TRUE)
