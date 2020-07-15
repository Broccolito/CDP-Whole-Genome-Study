rm(list = ls())
gc()

library(dplyr)

split_df = function(dff,n){
  d = dim(dff)[1]
  x = 1:d
  chunks = split(x, factor(sort(rank(x)%%n)))
  chunk_df = lapply(chunks, function(i){
    return(dff[i,])
  })
  names(chunk_df) = NULL
  return(chunk_df)
}

load("all_variant_cms_gt6.RData")
lit_search_genes = read.csv("lit_search_gene_list.csv")
unique_gene_list = data.frame(gene.symbol = unique(lit_search_genes$genes))

lit_search_overlap = all_variants %>%
inner_join(all_variants, unique_gene_list, by = "gene.symbol")
