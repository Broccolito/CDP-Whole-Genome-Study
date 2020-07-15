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

load("cdp_selected.rds")
load("wgs_physio.rds")

wgs = wgs_physio %>%
  select(wgs_id,id,initials,sex:abpm_conventional_dbp) %>%
  mutate_if(is.numeric, as.numeric) %>%
  cbind.data.frame(cdp_selected)

wgsm = subset(wgs, sex == "M")
wgsf = subset(wgs, sex == "F")

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

phyno_geno_corr_overall = run_regression(wgs)
phyno_geno_corr_male = run_regression(wgsm)
phyno_geno_corr_female = run_regression(wgsf)

write.csv(phyno_geno_corr_overall, file = "phyno_geno_corr_overall.csv",
          quote = FALSE, row.names = FALSE)
write.csv(phyno_geno_corr_male, file = "phyno_geno_corr_male.csv",
          quote = FALSE, row.names = FALSE)
write.csv(phyno_geno_corr_female, file = "phyno_geno_corr_female.csv",
          quote = FALSE, row.names = FALSE)

