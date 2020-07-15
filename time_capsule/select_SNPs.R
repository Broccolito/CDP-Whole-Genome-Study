library(dplyr)
library(beepr)

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
cdp2015 = read.csv("2015_aug_dec_2016_dec.csv")
wgs2018 = read.csv("2018_wgs.csv")
cms_snp = read.csv("andeans_iva_cms_score.csv")

wgs_physio = cdp[,1,drop = FALSE] %>%
  left_join(wgs2018, by = "wgs_id") %>%
  select(wgs_id,id) %>%
  left_join(cdp2015, by = "id")

cat("Locating SNPs in CDP genome...\n")

cdp_selected = cdp[,sapply(cms_snp$loc, function(x){
  tryCatch({
    cat(paste0("Locating ", x, " in CDP genome...\n"))
    pos = which(names(cdp)==x)
  },error = function(e){
    cat(paste0("Error occurred trying to locate", x, "\n"))
    pos = 1
  })
  return(pos)
})]

save(cdp_selected, file = "cdp_selected.rds")
save(wgs_physio, file = "wgs_physio.rds")
