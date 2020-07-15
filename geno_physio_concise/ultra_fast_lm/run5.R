source("ultra_fast_lm.R")
file_to_variant_list = function(file_path){
  f = read.delim(file = file_path)
  v= paste(f[,1],f[,2],sep = "_")
  return(v)
}

vlist5 = file_to_variant_list("variant_list/cms_BF6_splitup_5_IVAexport.txt")
run_regression(variant_list = vlist5) %>%
  write.csv(file = "vlist5.csv",quote = FALSE,row.names = FALSE)