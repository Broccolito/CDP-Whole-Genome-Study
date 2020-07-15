library(dplyr)
library(compiler)

rm(list = ls())
gc()

load("cdp.rds")
pheno = read.csv("2015_aug_dec_2016_dec.csv")
geno = read.csv("2018_wgs.csv")

wgs = cdp[,1,drop = FALSE] %>%
  left_join(geno,by = "wgs_id") %>%
  select(wgs_id,id) %>%
  left_join(pheno, by = "id") %>%
  select(wgs_id,id,sex,hct_vena_mean) %>%
  cbind.data.frame(cdp[,-1])

rm(pheno,geno,cdp)

get_regression_stats = function(variant,phenotype = "hct_vena_mean"){
  regression_stats = data.frame(phenotype=NULL,variant=NULL,pm=NULL,pf=NULL,
                                bm=NULL,bf=NULL,nm0=NULL,nm1=NULL,nm2=NULL,
                                nf0=NULL,nf1=NULL,nf2=NULL)
  try({
    gp = data.frame(pheno = wgs[[phenotype]],
                    geno = wgs[[variant]],sex = wgs[["sex"]]) %>%
      na.omit() %>%
      mutate(pheno = as.numeric(pheno)) %>%
      mutate(geno = as.numeric(geno)) %>%
      mutate(sex = as.factor(sex))
    gpm = filter(gp,sex == "M")
    gpf = filter(gp,sex == "F")
    fm = pheno ~ geno
    lrm = summary(lm(data = gpm,formula = fm))
    lrf = summary(lm(data = gpf,formula = fm))
    pm = lrm$coefficients[2,4]
    pf = lrf$coefficients[2,4]
    if(any(pm<=0.05,pf<=0.05)){
      bm = lrm$coefficients[2,1]
      bf = lrf$coefficients[2,1]
      nm0 = sum(gpm$geno==0)
      nm1 = sum(gpm$geno==1)
      nm2 = sum(gpm$geno==2)
      nf0 = sum(gpf$geno==0)
      nf1 = sum(gpf$geno==1)
      nf2 = sum(gpf$geno==2)
      regression_stats = data.frame(phenotype=phenotype,variant=variant,pm=pm,pf=pf,
                                    bm=bm,bf=bf,nm0=nm0,nm1=nm1,nm2=nm2,
                                    nf0=nf0,nf1=nf1,nf2=nf2)
    }else{
      stop("Not significant...")
    }
  },silent = TRUE)
  return(regression_stats)
}

variants_apriori = read.csv("interested_variants/james_list.csv") %>%
  select(variant) %>%
  unlist() %>%
  as.character()

variants_exonic = read.csv("interested_variants/andeans_iva_cms_score.csv") %>%
  select(variant) %>%
  unlist() %>%
  as.character()

variants_cms_gt6_lit_search = read.csv("interested_variants/cms_gt6_hypoxia_lit_search_overlap_IVAexport.csv") %>%
  select(variant) %>%
  unlist() %>%
  as.character()

variants_cms_gt8.1 = read.csv("interested_variants/cms_gt8.1_IVAexport.csv") %>%
  select(variant) %>%
  unlist() %>%
  as.character()

variants_cms_gt8.2 = read.csv("interested_variants/cms_gt8.2_IVAexport.csv") %>%
  select(variant) %>%
  unlist() %>%
  as.character()

regress_variants_list = function(variants_list){
  n = 1
  reg_stats_mat = vector()
  for(v in variants_list){
    try({
      reg_stats_mat = rbind.data.frame(reg_stats_mat,get_regression_stats(v))
      n = n + 1
      if(n%%100==0){
        cat(paste0(n," Regressions Evaluated...\n"))
      }
    },silent = TRUE)
  }
  return(reg_stats_mat)
}

results_apriori = regress_variants_list(variants_apriori)
results_exonic = regress_variants_list(variants_exonic)
results_cms_gt6_lit_search = regress_variants_list(variants_cms_gt6_lit_search)
results_cms_gt8.1 = regress_variants_list(variants_cms_gt8.1)
results_cms_gt8.2 = regress_variants_list(variants_cms_gt8.2)

save(results_apriori,results_exonic,
     results_cms_gt6_lit_search,
     results_cms_gt8.1,results_cms_gt8.2,file = "results_apriori.RData")

rm(list = ls())
gc()