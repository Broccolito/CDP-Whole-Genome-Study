library(dplyr)
source("regression_engine.R")

run_regressions("variant_list/andeans_iva_cms_score.csv") %>%
  write.csv("results/regression_iva_exon.csv",quote = FALSE,row.names = FALSE)