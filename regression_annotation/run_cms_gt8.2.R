library(dplyr)
source("regression_engine.R")

run_regressions("variant_list/cms_gt8.2_IVAexport.csv") %>%
  write.csv("results/regression_cms_gt8.2.csv",quote = FALSE,row.names = FALSE)