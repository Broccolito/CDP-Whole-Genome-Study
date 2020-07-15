rm(list = ls())
gc()

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

source("boxplot_rendering_engine.R")

apriori_boxplots = read.csv("apriori_boxplots.csv") %>%
  mutate(variant = as.character(variant)) %>%
  mutate(phenotype = as.character(phenotype))

exonic_boxplots = read.csv("exonic_boxplots.csv") %>%
  mutate(variant = as.character(variant)) %>%
  mutate(phenotype = as.character(phenotype))

for(i in 1:dim(apriori_boxplots)[1]){
  plot_pheno_geno(phenotype = apriori_boxplots$phenotype[i],
                  variant = apriori_boxplots$variant[i],
                  phenotype_name = apriori_boxplots$phenotype_name[i],
                  gene_name = apriori_boxplots$gene[i])
}
# try({
#   for(i in 1:dim(exonic_boxplots)[1]){
#     plot_pheno_geno(phenotype = exonic_boxplots$phenotype[i],
#                     variant = exonic_boxplots$variant[i],
#                     phenotype_name = exonic_boxplots$phenotype_name[i],
#                     gene_name = exonic_boxplots$gene[i])
#   }
# },silent = TRUE)




  
  
  
  
  
  

