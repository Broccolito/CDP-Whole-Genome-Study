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

gene_list_other_study = read.csv("genes_from_other_studies.csv") %>%
  mutate(gene_symbol = gene_name)
geno_physio = read.csv("andeans_iva_geno_physio.csv") %>%
  left_join(gene_list_other_study,by = "gene_symbol")

# write.csv(geno_physio, file = "overlaps.csv",quote = FALSE,row.names = FALSE)
