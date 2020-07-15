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

load("cdp.rds")
james_list = read.csv("pheno_geno_james.csv")
iva_list = read.csv("pheno_geno_master.csv")
cdp2015 = read.csv("2015_aug_dec_2016_dec.csv")
wgs2018 = read.csv("2018_wgs.csv")

selected_variant = c(as.character(iva_list$variant),
                     as.character(james_list$variant)) %>%
  unique()

cdp_selected_impute = cdp[,unlist(sapply(selected_variant, function(x){
  tryCatch({
    cat(paste0("Locating ", x, " in CDP genome...\n"))
    pos = which(names(cdp)==x)
  },error = function(e){
    cat(paste0("Error occurred trying to locate", x, "\n"))
    pos = 1
  })
  return(pos)
}))]

wgs = cdp[,1,drop = FALSE] %>%
  left_join(wgs2018, by = "wgs_id") %>%
  select(wgs_id,id) %>%
  left_join(cdp2015, by = "id") %>%
  cbind.data.frame(cdp_selected_impute) %>%
  mutate_if(is.numeric, as.numeric)

get_imputed_p_list = function(snp_list){
  imputed_p_list = vector()
  geno_physio_list = select(snp_list,phenotype,variant)
  for(i in 1:dim(geno_physio_list)[1]){
    try({
      
      p = geno_physio_list$phenotype[i]
      g = geno_physio_list$variant[i]
      
      geno_physio = select(wgs,sex) %>%
        mutate(geno = cdp[[as.character(g)]]) %>%
        mutate(geno = ifelse(geno==0,0,1)) %>%
        mutate(pheno = unlist(select(wgs,all_of(p)))) %>%
        na.omit() 
      
      geno_physio_male = filter(geno_physio,sex=="M")
      geno_physio_female = filter(geno_physio,sex=="F")
      
      gp_test = data.frame(
        phenotype = p,
        variant = g
      )
      
      gp_test_m = data.frame(
        Nm = dim(geno_physio_male)[1],
        Bm = NA,
        Pm = NA,
        Naltm = sum(geno_physio_male$geno),
        Nrefm = dim(geno_physio_male)[1] - sum(geno_physio_male$geno)
      )
      try({
        mt = t.test(pheno~geno,data = geno_physio_male)
        gp_test_m = data.frame(
          Nm = dim(geno_physio_male)[1],
          Bm = diff(mt$estimate),
          Pm = mt$p.value,
          Naltm = sum(geno_physio_male$geno),
          Nrefm = dim(geno_physio_male)[1] - sum(geno_physio_male$geno)
        )
      },silent = TRUE)
      
      gp_test_f = data.frame(
        Nf = dim(geno_physio_female)[1],
        Bf = NA,
        Pf = NA,
        Naltf = sum(geno_physio_female$geno!=0),
        Nreff = dim(geno_physio_female)[1] - sum(geno_physio_female$geno!=0)
      )
      try({
        ft = t.test(pheno~geno,data = geno_physio_female)
        gp_test_f = data.frame(
          Nf = dim(geno_physio_female)[1],
          Bf = diff(ft$estimate),
          Pf = ft$p.value,
          Naltf = sum(geno_physio_female$geno!=0),
          Nreff = dim(geno_physio_female)[1] - sum(geno_physio_female$geno!=0)
        )
      },silent = TRUE)
      
      gp_test = cbind.data.frame(gp_test,gp_test_m,gp_test_f) %>%
        mutate(
          adj_pi = (-log(ifelse(is.na(Pm),1,Pm),base = 10)) +
            (-log(ifelse(is.na(Pf),1,Pf),base = 10))
        )
      
      imputed_p_list = rbind.data.frame(imputed_p_list,gp_test)
      row.names(imputed_p_list) = NULL
      cat(paste0("Testing ",p," by ",g,"\n"))
      
    },silent = TRUE)
  }
  
  return(imputed_p_list)
  
}

iva_list = get_imputed_p_list(iva_list)
james_list = get_imputed_p_list(james_list)

iva_list = iva_list %>%
  mutate(geno_physio = paste(phenotype,variant))

james_list = james_list %>%
  mutate(geno_physio = paste(phenotype,variant))

pheno_geno_master_imputed = iva_list %>%
  left_join(read.csv("pheno_geno_master.csv") %>%
              mutate(geno_physio = paste(phenotype,variant)) %>%
              select(geno_physio,variation_type:cms_score),
            by = "geno_physio"
            ) %>%
  arrange(desc(adj_pi)) %>%
  select(-geno_physio)

pheno_geno_james_imputed = james_list %>%
  left_join(read.csv("pheno_geno_james.csv") %>%
              mutate(geno_physio = paste(phenotype,variant)) %>%
              select(geno_physio,gene,base_change),
            by = "geno_physio"
            ) %>%
  arrange(desc(adj_pi)) %>%
  select(-geno_physio)

write.csv(pheno_geno_master_imputed,
          file = "pheno_geno_master_imputed.csv",
          quote = FALSE,row.names = FALSE)
write.csv(pheno_geno_james_imputed,
          file = "pheno_geno_james_imputed.csv",
          quote = FALSE,row.names = FALSE)