suppressPackageStartupMessages({
  library(dplyr)
})

load("wgs_concise_physio.RData")

formulas =   paste0(c('hct_vena_mean',
                      'cms_score',
                      'spo2_mean',
                      'hr_mean',
                      'sbp_mean',
                      'dbp_mean',
                      'co_ppm',
                      'glucose',
                      'insulin',
                      'cholesterol',
                      'hdl',
                      'ldl',
                      'triglycerides',
                      'ferritin',
                      'iron',
                      'transferrin',
                      'total_testosterone',
                      'free_testosterone',
                      'height',
                      'hvr_corrected ',
                      'h_hvr_corrected_fixed',
                      'h_hvr_corrected',
                      'hcvr_corrected',
                      'hrr_hypox',
                      'hrr_co2',
                      'tv_ra',
                      'f_ra',
                      'etco2_ra',
                      'sao2_ra',
                      'hr_ra',
                      'peco2_ra',
                      'tv_21c',
                      'f_21c',
                      'etco2_21c',
                      'fio2_21c',
                      'sao2_21c',
                      'hr_21c',
                      'vibtps_kg_21c',
                      'peco2_21c',
                      'peco2_fix_21c',
                      'prdi_b_sleep',
                      'pahi_b_sleep',
                      'odi_b_sleep',
                      'meansatuation_b_sleep',
                      'minsaturation_b_sleep',
                      'maxsaturation_b_sleep',
                      'totalnumberofdesaturations_b_sleep',
                      'satbelow85_b_sleep',
                      'satbelow80_b_sleep',
                      'sleepefficiency_b_sleep',
                      'numberofwakes_b_sleep',
                      'vt_sleep',
                      'frequency_sleep',
                      'sao2_sleep',
                      'hr_bpm_sleep',
                      'vi_lmin_sleep',
                      'vibtps_sleep',
                      'vibtpskg_sleep',
                      'peco2_mmhg_sleep',
                      'hvr_sleep',
                      'hypercapnic_hvr_sleep',
                      'hcvr_sleep',
                      'hvr_corrected__sleep',
                      'h_hvr_corrected_sleep',
                      'hcvr_corrected_sleep',
                      'hrr_hypox_sleep',
                      'hrr_co2_sleep',
                      'sleepefficiency_sleep',
                      'totalahi_sleep',
                      'totalari_sleep',
                      'nadir_desat_sleep_sleep',
                      'spo2below85p_sleep_prct_sleep',
                      'spo2below80p_sleep_prct_sleep',
                      'spo2mean_sleep_sleep',
                      'odi_sleep_sleep',
                      'abpm_sbp_24h',
                      'abpm_dbp_24h',
                      'abpm_map_24h',
                      'abpm_sbp_vig',
                      'abpm_dbp_vig',
                      'abpm_map_vig',
                      'abpm_sbp_sleep',
                      'abpm_dbp_sleep',
                      'abpm_map_sleep',
                      'abpm_conventional_sbp',
                      'abpm_conventional_dbp'
),
" ~ age + ","geno")

llm = function(formula,data_male,data_female,variant){
  
  nm0 = sum(data_male$geno==0,na.rm = TRUE)
  nm1 = sum(data_male$geno==1,na.rm = TRUE)
  nm2 = sum(data_male$geno==2,na.rm = TRUE)
  p.value_male = 1
  effsize_male = 0
  
  nf0 = sum(data_female$geno==0,na.rm = TRUE)
  nf1 = sum(data_female$geno==1,na.rm = TRUE)
  nf2 = sum(data_female$geno==2,na.rm = TRUE)
  p.value_female = 1
  effsize_female = 0
  
  lmr = paste(nm0,nm1,nm2,"1___0",sep = "___")
  lfr = paste(nf0,nf1,nf2,"1___0",sep = "___")
  
  try({
    l_male = lm(formula = formula,data = data_male,na.action = "na.omit") %>%
      summary()
    p.value_male = l_male$coefficients[3,4]
    effsize_male = l_male$coefficients[3,1]
    lmr = paste(nm0,nm1,nm2,p.value_male,effsize_male,sep = "___")
  },silent = TRUE)
  
  try({
    l_female = lm(formula = formula,data = data_female,na.action = "na.omit") %>%
      summary()
    p.value_female = l_female$coefficients[3,4]
    effsize_female = l_female$coefficients[3,1]
    lfr = paste(nf0,nf1,nf2,p.value_female,effsize_female,sep = "___")
  },silent = TRUE)
  
  pheno = unlist(strsplit(formula," ~ "))[1]
  lr = paste(pheno,variant,lmr,lfr,sep = "___")
  
  return(lr)
}

linear_regression = function(variant){
  
  sub_wgs_male = data.frame(geno = wgsm_geno[[variant]]) %>%
    cbind.data.frame(wgsm_physio)
  sub_wgs_female = data.frame(geno = wgsf_geno[[variant]]) %>%
    cbind.data.frame(wgsf_physio)
  
  cat(paste0("Running Regressions on ", variant,"\n"))
  return(sapply(formulas,llm,sub_wgs_male,sub_wgs_female,variant))
  
}

run_regression = function(variant_list){
  omitted_variants = ""
  if(length(which(!variant_list%in%names(wgsm_geno)))!=0){
    omitted_variants = variant_list[which(!variant_list%in%names(wgsm_geno))]
    variant_list = variant_list[-which(!variant_list%in%names(wgsm_geno))]
  }
  regression_stats = sapply(variant_list,linear_regression) %>%
    sapply(FUN = function(x){unlist(strsplit(x,"___"))}) %>%
    t() %>%
    as.data.frame()
  rownames(regression_stats) = NULL
  names(regression_stats) = c("phenotype","variant","nm0","nm1","nm2",
                              "pm","bm","nf0","nf1","nf2","pf","bf")
  cat("Omitted Vriants: \n")
  cat(omitted_variants)
  return(regression_stats)
}

# res = run_regression(variant_list)
