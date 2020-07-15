library(dplyr)
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

cms_gt6_hypoxia_lit_search_overlap = generate_wgs("variant_list/cms_gt6_hypoxia_lit_search_overlap_IVAexport.csv")
save(cms_gt6_hypoxia_lit_search_overlap,file = "intermediate/cms_gt6_hypoxia_lit_search_overlap_wgs.RData")

cms_gt8.1 = generate_wgs("variant_list/cms_gt8.1_IVAexport.csv")
save(cms_gt8.1,file = "intermediate/cms_gt8.1_wgs.RData")

cms_gt8.2 = generate_wgs("variant_list/cms_gt8.2_IVAexport.csv")
save(cms_gt8.2,file = "intermediate/cms_gt8.2_wgs.RData")

# generate_wgs("variant_list/cms_gt6_hypoxia_lit_search_overlap_IVAexport.csv") %>%
#   save(file = "intermediate/cms_gt6_hypoxia_lit_search_overlap_wgs.RData")
# 
# generate_wgs("variant_list/cms_gt8.1_IVAexport.csv") %>%
#   save(file = "intermediate/cms_gt8.1_wgs.RData")
# 
# generate_wgs("variant_list/cms_gt8.2_IVAexport.csv") %>%
#   save(file = "intermediate/cms_gt8.2_wgs.RData")