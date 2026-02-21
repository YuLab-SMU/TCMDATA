#' download target genes associated with DN from Open Targets Platform
dir <- "/Users/hinna/Desktop/yulab/casestudy/OT-EFO_0000401-associated-targets-2_21_2026-v25_12.tsv"
dn_otp <- read.delim(dir)
dn_otp <- dn_otp$symbol
usethis::use_data(dn_otp, overwrite = TRUE)
