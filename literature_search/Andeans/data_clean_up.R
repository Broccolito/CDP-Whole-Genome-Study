library(dplyr)
library(readxl)

seperate_genes = function(x){
  x = gsub(pattern = ", ",replacement = " ",x)
  x = gsub(pattern = ",",replacement = " ",x)
  x = gsub(pattern = ";",replacement = " ",x)
  x = gsub(pattern = "; ",replacement = " ",x)
  x = unlist(strsplit(x," "))
  x = gsub(pattern = " ",replacement = "",x)
  return(x)
}

seperate_gene_lists = function(x){
  x = unlist(sapply(x,seperate_genes))
  names(x) = NULL
  x = x[!is.na(x)]
  x = unique(x)
  return(x)
}

make_data_frame = function(l,relevance = ""){
  total_lit_df = vector()
  for(i in 1:length(l)){
    study = names(l[i])
    gene_list = l[[i]]
    lit_df = data.frame(genes = gene_list,
                        study = rep(study,length(gene_list)),
                        relevance = relevance
    )
    total_lit_df = rbind.data.frame(total_lit_df,lit_df)
  }
  return(total_lit_df)
}

# Read Data
s1 = read_excel(path = "overlap.xlsx",sheet = "Selec_And") %>% as.list()
s2 = read_excel(path = "overlap.xlsx",sheet = "Selec_Tib") %>% as.list()
s3 = read_excel(path = "overlap.xlsx",sheet = "Selec_Ethiop") %>% as.list()
s4 = read_excel(path = "overlap.xlsx",sheet = "Selec_mixcomp") %>% as.list()
s5 = read_excel(path = "overlap.xlsx",sheet = "IPA") %>% as.list()
s6 = read_excel(path = "overlap.xlsx",sheet = "GO") %>% as.list()
s7 = read_excel(path = "overlap.xlsx",sheet = "RBC_lit") %>% as.list()
s8 = read_excel(path = "overlap.xlsx",sheet = "Lung_fxn_lit") %>% as.list()
s9 = read_excel(path = "overlap.xlsx",sheet = "Cardiovascular_lit") %>% as.list()
s10 = read_excel(path = "overlap.xlsx",sheet = "Metabolism_lit") %>% as.list()
s11 = read_excel(path = "overlap.xlsx",sheet = "apriorilists") %>% as.list()

# Load as processed data frames
s1 = lapply(s1,FUN = seperate_gene_lists) %>% make_data_frame("Selec_And")
s2 = lapply(s2,FUN = seperate_gene_lists) %>% make_data_frame("Selec_Tib")
s3 = lapply(s3,FUN = seperate_gene_lists) %>% make_data_frame("Selec_Ethiop")
s4 = lapply(s4,FUN = seperate_gene_lists) %>% make_data_frame("Selec_mixcomp")
s5 = lapply(s5,FUN = seperate_gene_lists) %>% make_data_frame("IPA")
s6 = lapply(s6,FUN = seperate_gene_lists) %>% make_data_frame("GO")
s7 = lapply(s7,FUN = seperate_gene_lists) %>% make_data_frame("RBC_lit")
s8 = lapply(s8,FUN = seperate_gene_lists) %>% make_data_frame("Lung_fxn_lit")
s9 = lapply(s9,FUN = seperate_gene_lists) %>% make_data_frame("Cardiovascular_lit")
s10 = lapply(s10,FUN = seperate_gene_lists) %>% make_data_frame("Metabolism_lit")
s11 = lapply(s11,FUN = seperate_gene_lists) %>% make_data_frame("apriorilists")

lit_search_gene_list = rbind.data.frame(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11)
write.csv(lit_search_gene_list,file = "lit_search_gene_list.csv",
          quote = FALSE,row.names = FALSE)