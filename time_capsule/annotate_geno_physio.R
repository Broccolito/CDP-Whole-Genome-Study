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

pheno_geno_overall = read.csv("phyno_geno_corr_overall.csv")
pheno_geno_male = read.csv("phyno_geno_corr_male.csv")
pheno_geno_female = read.csv("phyno_geno_corr_female.csv")
iva_variant = read.csv("andeans_iva_cms_score.csv")

pheno_geno_male = mutate(pheno_geno_male,Nm = sample_size) %>%
  mutate(Pm = p_value) %>%
  mutate(Bm = effect_size) %>%
  mutate(pheno_geno = paste(phynotype,variant)) %>%
  select(-sample_size,-p_value,-effect_size,-phynotype,-variant,-snp_occurrence)

pheno_geno_female = mutate(pheno_geno_female,Nf = sample_size) %>%
  mutate(Pf = p_value) %>%
  mutate(Bf = effect_size) %>%
  mutate(pheno_geno = paste(phynotype,variant)) %>%
  select(-sample_size,-p_value,-effect_size,-phynotype,-variant,-snp_occurrence)

pheno_geno = full_join(pheno_geno_male,pheno_geno_female,by = "pheno_geno") %>%
  mutate(phenotype = sapply(pheno_geno, function(x){
    unlist(strsplit(x," "))[1]
  })) %>%
  mutate(variant = sapply(pheno_geno, function(x){
    unlist(strsplit(x," "))[2]
  })) %>%
  select(-pheno_geno) %>%
  select(phenotype,variant,Nm,Pm,Bm,Nf,Pf,Bf) %>%
  left_join(mutate(iva_variant,variant = loc),by = "variant") %>%
  select(-chromosome,-position,-gene_names,-(X1000_genomes_frequency:clinvar_id),
         -reference_allele,-sample_allele,-transcript_id,-variant_findings,
         -classification,-transcriptional_regulator,-regulator,-loc,
         -(tfbs_enhancer_binding_score:fused_to)) %>%
  mutate(adj_pi = (-log(ifelse(is.na(Pm),1,Pm))) + (-log(ifelse(is.na(Pf),1,Pf)))) %>%
  arrange(desc(adj_pi)) %>%
  select(phenotype,variant,adj_pi,everything())

write.csv(pheno_geno,file = "pheno_geno_master.csv",quote = FALSE,row.names = FALSE)

