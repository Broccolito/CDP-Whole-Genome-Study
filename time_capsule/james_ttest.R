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

run_regression = function(wgs_subset){
  df_physio = select(wgs_subset, wgs_id:abpm_conventional_dbp)
  df_genetics = select(wgs_subset, -(wgs_id:abpm_conventional_dbp))
  result_df = vector()
  for(i in names(df_physio)){
    for(g in names(df_genetics)){
      if(i!="wgs_id"&i!="initials"&i!="id"){
        try({
          xy = data.frame(
            y = df_physio[[i]],
            x = df_genetics[[g]]
          ) %>%
            na.omit() %>%
            mutate(x = ifelse(x!=0,1,0))
          ttest = t.test(y~x,data = xy)
          n = dim(xy)[1]
          pvalue = ttest$p.value
          effsize = diff(ttest$estimate)
          if(pvalue<=0.05){
            result = data.frame(phynotype = i,
                                variant = g,
                                sample_size = n,
                                p_value = pvalue,
                                effect_size = effsize)
            row.names(result) = NULL
            result_df = rbind.data.frame(result_df,result)
            cat(paste0("Find Significant Correlations between ",
                       i, " and ", g,"\n"))
          }
        },silent = TRUE)
      }
    }
  }
  occurrence = as.data.frame(table(result_df$variant))
  snp_occurrence_df = data.frame(variant = occurrence[,1],
                                 snp_occurrence = occurrence[,2])
  result_df = left_join(result_df, snp_occurrence_df, by = "variant")
  return(result_df)
}

wd = get_directory()
setwd(wd)

load("cdp.rds")
cdp2015 = read.csv("2015_aug_dec_2016_dec.csv")
wgs2018 = read.csv("2018_wgs.csv")
cms_snp_james = read.csv("james_list.csv")

wgs_physio = cdp[,1,drop = FALSE] %>%
  left_join(wgs2018, by = "wgs_id") %>%
  select(wgs_id,id) %>%
  left_join(cdp2015, by = "id")

cdp_selected_james = cdp[,unlist(sapply(cms_snp_james$variant, function(x){
  tryCatch({
    cat(paste0("Locating ", x, " in CDP genome...\n"))
    pos = which(names(cdp)==x)
  },error = function(e){
    cat(paste0("Error occurred trying to locate", x, "\n"))
    pos = 1
  })
  return(pos)
}))]

wgs = wgs_physio %>%
  select(wgs_id,id,initials,sex:abpm_conventional_dbp) %>%
  mutate_if(is.numeric, as.numeric) %>%
  cbind.data.frame(cdp_selected_james)

wgsm = subset(wgs, sex == "M")
wgsf = subset(wgs, sex == "F")

james_phyno_geno_corr_overall = run_regression(wgs)
james_phyno_geno_corr_male = run_regression(wgsm)
james_phyno_geno_corr_female = run_regression(wgsf)

pheno_geno_male = mutate(james_phyno_geno_corr_male,Nm = sample_size) %>%
  mutate(Pm = p_value) %>%
  mutate(Bm = effect_size) %>%
  mutate(pheno_geno = paste(phynotype,variant)) %>%
  select(-sample_size,-p_value,-effect_size,-phynotype,-variant,-snp_occurrence)

pheno_geno_female = mutate(james_phyno_geno_corr_female,Nf = sample_size) %>%
  mutate(Pf = p_value) %>%
  mutate(Bf = effect_size) %>%
  mutate(pheno_geno = paste(phynotype,variant)) %>%
  select(-sample_size,-p_value,-effect_size,-phynotype,-variant,-snp_occurrence)

pheno_geno_james = full_join(pheno_geno_male,pheno_geno_female,by = "pheno_geno") %>%
  mutate(phenotype = sapply(pheno_geno, function(x){
    unlist(strsplit(x," "))[1]
  })) %>%
  mutate(variant = sapply(pheno_geno, function(x){
    unlist(strsplit(x," "))[2]
  })) %>%
  select(-pheno_geno) %>%
  select(phenotype,variant,Nm,Pm,Bm,Nf,Pf,Bf) %>%
  left_join(cms_snp_james,by = "variant") %>%
  # select(-chromosome,-position,-gene_names,-(X1000_genomes_frequency:clinvar_id),
  #        -reference_allele,-sample_allele,-transcript_id,-variant_findings,
  #        -classification,-transcriptional_regulator,-regulator,-loc,
  #        -(tfbs_enhancer_binding_score:fused_to)) %>%
  mutate(adj_pi = (-log(ifelse(is.na(Pm),1,Pm))) + (-log(ifelse(is.na(Pf),1,Pf)))) %>%
  arrange(desc(adj_pi)) %>%
  select(phenotype,variant,adj_pi,everything())

write.csv(pheno_geno_james,file = "pheno_geno_james.csv",quote = FALSE,row.names = FALSE)
