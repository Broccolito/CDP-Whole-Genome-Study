library(dplyr)
source("regression_engine.R")

run_regressions("variant_list/cms_gt6_hypoxia_lit_search_overlap_IVAexport.csv") %>%
  write.csv("results/regression_gt6_hypoxia_lit_search_overlap.csv",quote = FALSE,row.names = FALSE)

