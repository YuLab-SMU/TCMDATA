## load ML_preprocessed data
setwd("/Users/hinna/TCMDATA/data-raw")
covid19 <- readRDS("/Users/hinna/Desktop/yulab/casestudy/ml_covid_bulk_demo.rds")
usethis::use_data(covid19, overwrite = TRUE)
