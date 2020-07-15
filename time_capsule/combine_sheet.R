library(dplyr)

rm(list = ls())

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

geno_physio_overall = read.csv("phyno_geno_corr_overall.csv")
geno_physio_male = read.csv("phyno_geno_corr_male.csv")
geno_physio_female = read.csv("phyno_geno_corr_female.csv")
iva_variant_list = read.csv("andeans_iva_cms_score.csv")

generate_geno_physio_list = function(geno_physio){
  geno_physio_list = data.frame(
    loc = unique(geno_physio$variant),
    related_physio = sapply(unique(geno_physio$variant),function(x){
      filter(geno_physio, variant == x) %>%
        select(phynotype) %>%
        unlist() %>%
        as.character() %>%
        paste(collapse = " ")
    })
  )
  return(geno_physio_list)
}

geno_physio_overall = generate_geno_physio_list(geno_physio_overall)
geno_physio_male = generate_geno_physio_list(geno_physio_male)
geno_physio_female = generate_geno_physio_list(geno_physio_female)

iva_variant_list = iva_variant_list %>%
  full_join(geno_physio_overall,by="loc",suffix = c("","_overall")) %>%
  full_join(geno_physio_male,by="loc",suffix = c("","_male")) %>%
  full_join(geno_physio_female,by="loc",suffix = c("","_female"))

write.csv(iva_variant_list, file = "andeans_iva_geno_physio.csv",
          quote = FALSE,row.names = FALSE)
