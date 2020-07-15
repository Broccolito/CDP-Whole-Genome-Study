library(dplyr)
library(compiler)

rm(list = ls())
load("intermediate/cms_gt6_hypoxia_lit_search_overlap_wgs.RData")
load("intermediate/cms_gt8.1_wgs.RData")
load("intermediate/cms_gt8.2_wgs.RData")

run_regression = function(wgs){
  pheno = wgs[,1:198] 
  geno = wgs[,199:dim(wgs)[2]]
  
  suppressWarnings({
    pheno_male = pheno[which(pheno$sex=="M"),] %>% 
      select(-wgs_id,-sex,-birthplace,-ancestry_gen,-ee,-bloodletting,-(smokes:address)) 
    pheno_names = names(pheno_male)
    pheno_male = pheno_male %>%
      as.matrix() %>% apply(MARGIN=2,FUN = as.numeric)
    
    geno_male = geno[which(pheno$sex=="M"),] 
    geno_names = names(geno_male)
    geno_male = geno_male %>%
      as.matrix()
    geno_male_t = geno_male
    geno_male_t[geno_male_t==2] = 1
    
    pheno_female = pheno[which(pheno$sex=="F"),] %>% 
      select(-wgs_id,-sex,-birthplace,-ancestry_gen,-ee,-bloodletting,-(smokes:address)) %>%
      as.matrix() %>% apply(MARGIN=2,FUN = as.numeric)
    geno_female = geno[which(pheno$sex=="F"),] %>%
      as.matrix()
    geno_female_t = geno_female
    geno_female_t[geno_female_t==2] = 1
  })
  
  n = 1
  cor_mat = matrix(NA,nrow = (length(pheno_names)*length(geno_names))/2, ncol = 12)
  for(i in 1:length(pheno_names)){
    cat(paste0("Finding variants correlated with ", pheno_names[i],"...\n"))
    for(j in 1:length(geno_names)){
      
      phm = pheno_male[,i]
      getm = geno_male_t[,j]
      gem = geno_male[,j]
      
      pmt = NA
      pm = NA
      bm = NA
      nm0 = NA
      nm1 = NA
      nm2 = NA
      try({
        ttest = t.test(phm~getm,na.action="na.omit")
        if(ttest$p.value<=0.05){
          pmt = ttest$p.value
          l = lm(phm~gem,na.action = "na.omit") %>%
            summary()
          pm = l$coefficients[2,4]
          bm = l$coefficients[2,1]
          nm0 = sum(gem==0)
          nm1 = sum(gem==1)
          nm2 = sum(gem==2)
        }
      },silent =TRUE)
      
      phf = pheno_female[,i]
      getf = geno_female_t[,j]
      gef = geno_female[,j]
      
      pft = NA
      pf = NA
      bf = NA
      nf0 = NA
      nf1 = NA
      nf2 = NA
      try({
        ttest = t.test(phf~getf,na.action="na.omit")
        if(ttest$p.value<=0.05){
          pft = ttest$p.value
          l = lm(phf~gef,na.action = "na.omit") %>%
            summary()
          pf = l$coefficients[2,4]
          bf = l$coefficients[2,1]
          nf0 = sum(gef==0)
          nf1 = sum(gef==1)
          nf2 = sum(gef==2)
        }
      },silent = TRUE)
      
      if(any(!is.na(pmt),!is.na(pft))){
        cor_mat[n,] = c(pheno_names[i],geno_names[j],nm0,nm1,nm2,nf0,nf1,nf2,pm,bm,pmt,pft)
        n = n + 1
        if(n %% 1000 == 0){
          cat(paste0("More than ",n, " significant correlations found...\n"))
        }
      }
      
    }
  }
  colnames(cor_mat) = c("phenotype","variant","nm0","nm1","nm2","nf0","nf1","nf2","pm","bm","pmt","pft")
  return(cor_mat)
}
run_regression =  cmpfun(run_regression)

# regression_cms_gt6_hypoxia_lit_search = run_regression(cms_gt6_hypoxia_lit_search_overlap)
# write.csv(regression_cms_gt6_hypoxia_lit_search,
#           file = "regression_cms_gt6_hypoxia_lit_search.csv",quote = FALSE,row.names = FALSE)
regression_cms_gt8.1 = run_regression(cms_gt8.1)
write.csv(regression_cms_gt8.1,file = "regression_cms_gt8.1.csv",quote = FALSE,row.names = FALSE)
regression_cms_gt8.2 = run_regression(cms_gt8.2)
write.csv(regression_cms_gt8.2,file = "regression_cms_gt8.2.csv",quote = FALSE,row.names = FALSE)



