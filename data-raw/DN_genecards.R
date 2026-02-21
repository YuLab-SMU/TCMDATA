## load DN targets from Genecards
library(dplyr)
dir <- "/Users/hinna/Desktop/yulab/casestudy/DN_GeneCards.csv"
dn_gcds <- read.csv(dir) %>%
  select(Gene.Symbol) %>%
  pull() %>%         
  as.character() %>%  
  unique()           

usethis::use_data(dn_gcds, overwrite = TRUE)
