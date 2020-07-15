library(dplyr)
source("regression_engine.R")

run_regressions("variant_list/cms_gt8.1_IVAexport.csv") %>%
  write.csv("results/regression_cms_gt8.1.csv",quote = FALSE,row.names = FALSE)