rm(list = ls())

library(dplyr)

get_directory = function(){
  args <- commandArgs(trailingOnly = FALSE)
  file <- "--file="
  rstudio <- "RStudio"
  
  match <- grep(rstudio, args)
  if(length(match) > 0){
    return(dirname(rstudioapi::getSourceEditorContext()$path))
  }else{
    match <- grep(file, args)
    if (length(match) > 0) {
      return(dirname(normalizePath(sub(file, "", args[match]))))
    }else{
      return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
    }
  }
}

wd = get_directory()
setwd(wd)

cms_list = read.csv("CDP_cms_list.csv")

cms_list = subset(cms_list,cms_score>=6) %>%
  mutate_if(is.factor,as.character)

gene_cms_score_list = vector()
for(i in 1:dim(cms_list)[1]){
  gene_cms_score = unlist(strsplit(cms_list$related_genes[i],",")) %>%
    paste(cms_list$cms_score[i],sep = "&") %>%
    paste(collapse = ",")
  gene_cms_score_list = c(gene_cms_score_list,gene_cms_score)
}

gene_cms_score_list = gene_cms_score_list %>% 
  strsplit(",") %>%
  unlist()

gene_cms_score_list = data.frame(gene_cms_score = gene_cms_score_list)

gene_cms_score_list = mutate(gene_cms_score_list,
                             cms_score = sapply(gene_cms_score_list$gene_cms_score,
                                                function(x){
                                                  x = as.character(x)
                                                  unlist(strsplit(x,"&"))[2]
                                                })) %>% 
  mutate(gene_symbol = sapply(gene_cms_score_list$gene_cms_score, 
                              function(x){
                                x = as.character(x)
                                unlist(strsplit(x,"&"))[1]
                              })) %>%
  select(gene_symbol,cms_score) %>%
  mutate(cms_score = as.numeric(cms_score)) %>%
  mutate(gene_symbol = as.character(gene_symbol)) %>%
  arrange(desc(cms_score)) %>%
  filter(gene_symbol != "NA")

# Remove duplicates from the list
for(i in 1:dim(gene_cms_score_list)[1]){
  if(i <= dim(gene_cms_score_list)[1]){
    if(length(which(gene_cms_score_list$gene_symbol == gene_cms_score_list$gene_symbol[i]))>=2){
      gene_cms_score_list = gene_cms_score_list[-which(gene_cms_score_list$gene_symbol==gene_cms_score_list$gene_symbol[i])[-1],]
    }
  }
}

# head(gene_cms_score_list)

literature_list = read.csv("literature_genes.csv")
literature_list = literature_list %>%
  select(Bigham2009:Wuren2014) %>%
  as.list()

literature_list = lapply(literature_list, function(x){
  x = x[x!=""]
  x = x[!is.na(x)]
  x = as.character(x)
  x = unlist(strsplit(x," "))
  x = unique(x)
})

lit_list = vector()
for(i in 1:length(literature_list)){
  lit = paste(literature_list[[i]],names(literature_list)[i],
        sep = "&&")
  lit_list = c(lit_list,lit)
}

lit_list = data.frame(
  gene_symbol = sapply(lit_list, function(x){
    unlist(strsplit(x,"&&"))[1]
  }),
  study = sapply(lit_list,function(x){
    unlist(strsplit(x,"&&"))[2]
  })
)

rownames(lit_list) = NULL

lit_list = lit_list %>%
  mutate(gene_symbol = gsub(" ",replacement = "",gene_symbol)) %>%
  mutate(gene_symbol = gsub(":",replacement = "",gene_symbol)) %>%
  filter(gene_symbol!="") %>%
  filter(!is.na(gene_symbol))

overlapping_genes = inner_join(gene_cms_score_list,lit_list,
                               by = "gene_symbol") %>%
  mutate(occurance = sapply(gene_symbol, function(x){
    length(which(x==gene_symbol))
  })) %>%
  arrange(desc(occurance),desc(cms_score))

cms_gt8 = filter(gene_cms_score_list,cms_score>=8)

write.csv(lit_list,"gene_lists/hypoxia_literature_search_genes.csv",
          quote = FALSE,row.names = FALSE)

write.csv(overlapping_genes,
          "gene_lists/cms_gt6_hypoxia_literature_search_overlap_gene_list.csv",
          quote = FALSE,row.names = FALSE)

write.table(unique(overlapping_genes$gene_symbol),
            file = "gene_lists/cms_gt6_hypoxia_literature_search_overlap_genes.txt",
            quote = FALSE,row.names = FALSE)

write.table(cms_gt8[1:2693,]$gene_symbol,
            file = "gene_lists/cms_gt8.1.txt",quote = FALSE,row.names = FALSE)
write.table(cms_gt8[2694:5386,]$gene_symbol,
            file = "gene_lists/cms_gt8.2.txt",quote = FALSE,row.names = FALSE)
