library(dplyr)
source("regression_engine.R")

run_regressions("variant_list/james_list.csv") %>%
  write.csv("results/regression_apriori.csv",quote = FALSE,row.names = FALSE)