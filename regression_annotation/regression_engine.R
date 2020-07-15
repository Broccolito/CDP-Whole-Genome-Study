run_regressions = function(gene_list){
  
  library(dplyr)
  
  run_tests = function(wgs){
    
    cat("Starting Regression Engine...\n")
    
    wgs_physio = dplyr::select(wgs,wgs_id:abpm_conventional_dbp)
    wgs_geno = dplyr::select(wgs,-(wgs_id:abpm_conventional_dbp))
    
    wgsm = dplyr::filter(wgs,sex=="M")
    wgsf = dplyr::filter(wgs,sex=="F")
    
    wgsm_physio = dplyr::select(wgsm,wgs_id:abpm_conventional_dbp)
    wgsm_geno = dplyr::select(wgsm,-(wgs_id:abpm_conventional_dbp))
    wgsf_physio = dplyr::select(wgsf,wgs_id:abpm_conventional_dbp)
    wgsf_geno = dplyr::select(wgsf,-(wgs_id:abpm_conventional_dbp))
    
    correlation_list = vector()
    for(i in names(wgs_physio)){
      for(j in names(wgs_geno)){
        if(i!="wgs_id"){
          try({

            xym = data.frame(physio = wgsm_physio[[i]],
                             geno = wgsm_geno[[j]]) %>%
              na.omit() %>%
              mutate(physio = as.numeric(physio)) %>%
              mutate(geno = as.numeric(geno))
            xyf = data.frame(physio = wgsf_physio[[i]],
                             geno = wgsf_geno[[j]]) %>%
              na.omit() %>%
              mutate(physio = as.numeric(physio)) %>%
              mutate(geno = as.numeric(geno))
            
            Nm = dim(xym)[1]
            Nf = dim(xyf)[1]
            Nalt0m = length(which(xym$geno==0))
            Nalt1m = length(which(xym$geno==1))
            Nalt2m = length(which(xym$geno==2))
            Nalt0f = length(which(xyf$geno==0))
            Nalt1f = length(which(xyf$geno==1))
            Nalt2f = length(which(xyf$geno==2))
            
            Pm = NA
            Bm = NA
            try({
              mlm = summary(lm(formula = physio ~ geno, data = xym))
              Pm = mlm$coefficients[2,4]
              Bm = mlm$coefficients[2,1]
            },silent = TRUE)
            
            Pf = NA
            Bf = NA
            try({
              flm = summary(lm(formula = physio ~ geno, data = xyf))
              Pf = flm$coefficients[2,4]
              Bf = flm$coefficients[2,1]
            },silent = TRUE)
            
            Pm_ttest = 1
            Bm_ttest = NA
            try({
              xymt = mutate(xym,geno = ifelse(geno==0,0,1))
              ttestm = t.test(physio ~ geno, xymt)
              Pm_ttest = ttestm$p.value
              Bm_ttest = diff(ttestm$estimate)
            },silent = TRUE)
            
            Pf_ttest = 1
            Bf_ttest = NA
            try({
              xyft = mutate(xyf,geno = ifelse(geno==0,0,1))
              ttestf = t.test(physio ~ geno, xyft)
              Pf_ttest = ttestf$p.value
              Bf_ttest = diff(ttestf$estimate)
            },silent = TRUE)
            
            if(any(Pm_ttest<=0.05,Pf_ttest<=0.05)){
              sig_corr = data.frame(phenotype = i,
                                    variant = j,
                                    Nm,Nf,Nalt0m,Nalt1m,Nalt2m,
                                    Nalt0f,Nalt1f,Nalt2f,Pm,Bm,Pf,Bf,
                                    Pm_ttest,Pf_ttest) %>%
                mutate(adj_pi = -log(ifelse(is.na(Pm_ttest),1,Pm_ttest),base = 10) +
                         -log(ifelse(is.na(Pf_ttest),1,Pf_ttest),base = 10))
              correlation_list = rbind.data.frame(correlation_list,sig_corr)
              cat(paste0("Significnat Correlation found between ",i," and ",j,"\n"))
            }
            
          },silent = TRUE)
        }
      }
    }
    correlation_list = correlation_list %>%
      arrange(desc(adj_pi))
    return(correlation_list)
  }
  
  generate_wgs = function(filepath = "andeans_iva_cms_score.csv"){
    
    cat("Loading Whole Genome Sequencing Data...\n")
    load("cdp.rds")
    cat("Loading Physiology Data...\n")
    cdp2015 = read.csv("2015_aug_dec_2016_dec.csv")
    wgs2018 = read.csv("2018_wgs.csv")
    cat("Loading Interested Variants...\n")
    gene_list = read.csv(filepath)
    
    cat("Selecting Interested Variants from Whole Genome Sequencing...\n")
    # selected_snp1 = cdp[,unlist(sapply(gene_list$variant, function(x){
    #   tryCatch({
    #     cat(paste0("Locating ", x, " in CDP genome...\n"))
    #     pos = which(names(cdp)==x)
    #   },error = function(e){
    #     cat(paste0("Error occurred trying to locate", x, "\n"))
    #     pos = 1
    #   })
    #   return(pos)
    # }))]
    selected_snp = cdp[,(names(cdp) %in% as.character(gene_list$variant))]
    
    suppressWarnings({
      wgs = cdp[,1,drop = FALSE] %>%
        left_join(wgs2018, by = "wgs_id") %>%
        dplyr::select(wgs_id,id) %>%
        left_join(cdp2015, by = "id") %>%
        dplyr::select(wgs_id,id,initials,sex:abpm_conventional_dbp) %>%
        mutate_if(is.numeric, as.numeric) %>%
        cbind.data.frame(selected_snp) %>%
        dplyr::select(-id,-initials,-birthdate)
    })
    
    rm(cdp)
    return(wgs)
    
  }
  
  return(
    generate_wgs(gene_list) %>%
      run_tests() %>%
      left_join(
        read.csv(gene_list),
        by = "variant"
      )
  )
}

