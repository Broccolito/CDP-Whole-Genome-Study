library(dplyr)
source("regression_engine.R")

run_regressions("variant_list/cms_iva_minimal_filter_sub5.csv") %>%
  write.csv("results/regression_iva_minimal_filter_sub5.csv",quote = FALSE,row.names = FALSE)