library(dplyr)

rm(list=ls())
load("results/results_apriori.RData")
load("results/results_exonic.RData")
load("results/results_cms_gt6_lit_search.RData")
load("results/results_cms_gt8.1.RData")
load("results/results_cms_gt8.2.RData")

variants_apriori = read.csv("interested_variants/james_list.csv")
variants_exonic = read.csv("interested_variants/andeans_iva_cms_score.csv")
variants_cms_gt6_lit_search = read.csv("interested_variants/cms_gt6_hypoxia_lit_search_overlap_IVAexport.csv")
variants_cms_gt8.1 = read.csv("interested_variants/cms_gt8.1_IVAexport.csv")
variants_cms_gt8.2 = read.csv("interested_variants/cms_gt8.2_IVAexport.csv")
gene_cms = read.csv("unique_gene_cms.csv") %>%
  mutate(gene_symbol = gene_names) %>%
  select(gene_symbol,cms_score)

suppressWarnings({
  results_apriori = left_join(results_apriori,variants_apriori,by = "variant") %>%
    mutate(gene_symbol = gene) %>%
    select(-gene) %>%
    left_join(gene_cms,by="gene_symbol")  %>%
    arrange(pm) %>%
    select(phenotype:nf2,gene_symbol,region,cms_score)
  
  results_exonic = left_join(results_exonic,variants_exonic,by = "variant") %>%
    left_join(gene_cms,by="gene_symbol",suffix = c("","x")) %>%
    arrange(pm) %>%
    select(phenotype:nf2,gene_symbol,gene_region,cms_score)
  
  results_cms_gt6_lit_search = left_join(results_cms_gt6_lit_search,variants_cms_gt6_lit_search,by = "variant") %>%
    left_join(gene_cms,by="gene_symbol")  %>%
    arrange(pm) %>%
    select(phenotype:nf2,gene_symbol,gene_region,cms_score)
  
  results_cms_gt8.1 = left_join(results_cms_gt8.1,variants_cms_gt8.1,by = "variant") %>%
    left_join(gene_cms,by="gene_symbol")  %>%
    arrange(pm) %>%
    select(phenotype:nf2,gene_symbol,gene_region,cms_score)
  
  results_cms_gt8.2 = left_join(results_cms_gt8.2,variants_cms_gt8.2,by = "variant") %>%
    left_join(gene_cms,by="gene_symbol")  %>%
    arrange(pm) %>%
    select(phenotype:nf2,gene_symbol,gene_region,cms_score)
})

write.csv(results_apriori,file = "annotated_results/results_apriori.csv",quote = FALSE,row.names = FALSE)
write.csv(results_exonic,file = "annotated_results/results_exonic.csv",quote = FALSE,row.names = FALSE)
write.csv(results_cms_gt6_lit_search,file = "annotated_results/results_cms_gt6_lit_search.csv",quote = FALSE,row.names = FALSE)
write.csv(results_cms_gt8.1,file = "annotated_results/results_cms_gt8.1.csv",quote = FALSE,row.names = FALSE)
write.csv(results_cms_gt8.2,file = "annotated_results/results_cms_gt8.2.csv",quote = FALSE,row.names = FALSE)
