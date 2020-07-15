library(dplyr)

rm(list = ls())

ipa = read.csv("ipa_biofilter1.csv")
lit_search_genes = read.csv("hypoxia_literature_search_genes.csv")
exonic_variants = read.csv("results_exonic.csv")
lit_search_variants = read.csv("results_cms_gt6_lit_search.csv")

occurrence = vector()
for(i in 1:dim(lit_search_genes)[1]){
  occurrence = c(occurrence,
                 sum(lit_search_genes$gene_symbol==lit_search_genes$gene_symbol[i]))
}

lit_search_genes = cbind.data.frame(lit_search_genes,occurrence) %>%
  arrange(desc(occurrence))

unique_gene_list = unique(lit_search_genes$gene_symbol)
mentioned_studies_list = vector()
for(gene in unique_gene_list){
  mentioned_studies = paste(filter(lit_search_genes,gene_symbol == gene)$study,
                            collapse = ";")
  mentioned_studies_list = c(mentioned_studies_list,mentioned_studies)
}
lit_search_genes = data.frame(gene_symbol = unique_gene_list, studies = mentioned_studies_list) %>%
  inner_join(lit_search_genes,by="gene_symbol") %>%
  select(gene_symbol,occurrence,studies) %>%
  distinct()

# write.csv(lit_search_genes,file = "unique_lit_search_genes.csv",quote = FALSE,row.names = FALSE)

annotated_exonic_variants = exonic_variants %>%
  left_join(ipa,by="gene_symbol") %>%
  left_join(lit_search_genes,by = "gene_symbol") %>%
  mutate(occurrence = ifelse(is.na(occurrence),0,occurrence)) %>%
  mutate(IPA = ifelse(is.na(IPA),0,1))

annotated_lit_search_variants = lit_search_variants %>%
  left_join(ipa,by="gene_symbol") %>%
  left_join(lit_search_genes,by = "gene_symbol") %>%
  mutate(occurrence = ifelse(is.na(occurrence),0,occurrence)) %>%
  mutate(IPA = ifelse(is.na(IPA),0,1))

write.csv(annotated_exonic_variants,file = "annotated_exonic_variants.csv",quote = FALSE,row.names = FALSE)
write.csv(annotated_lit_search_variants,file = "annotated_lit_search_variants.csv",quote = FALSE,row.names = FALSE)