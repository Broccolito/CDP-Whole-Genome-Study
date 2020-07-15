library(dplyr)

rm(list = ls())
gc()

cms_gt8 = rbind.data.frame(
  read.csv("results_cms_gt8.1.csv"),
  read.csv("results_cms_gt8.2.csv")
)

cms_gt6_lit_search = read.csv("results_cms_gt6_lit_search.csv")

ipa = read.csv("ipa_biofilter1.csv")

ipa_overlap_lit_search = inner_join(cms_gt6_lit_search,ipa,by = "gene_symbol") %>%
  arrange(pm)

write.csv(ipa_overlap_lit_search,file = "ipa_cms_gt6_lit_search_overlap.csv",quote = FALSE,row.names = FALSE)

ipa_overlap_cms_gt8 = inner_join(cms_gt8,ipa,by = "gene_symbol") %>%
  arrange(pm)

write.csv(ipa_overlap_cms_gt8,file = "ipa_ipa_overlap_cms_gt8_overlap.csv",quote = FALSE,row.names = FALSE)