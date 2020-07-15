rm(list = ls())
library(dplyr)

load("cdp_cms_gwas_regression_stats.RData")

v1 = read.delim(file = "variant_list/cms_BF6_splitup_1_IVAexport.txt")
v2 = read.delim(file = "variant_list/cms_BF6_splitup_2_IVAexport.txt")
v3 = read.delim(file = "variant_list/cms_BF6_splitup_3_IVAexport.txt")
v4 = read.delim(file = "variant_list/cms_BF6_splitup_4_IVAexport.txt")
v5 = read.delim(file = "variant_list/cms_BF6_splitup_5_IVAexport.txt")
v6 = read.delim(file = "variant_list/cms_BF6_splitup_6_IVAexport.txt")
v7 = read.delim(file = "variant_list/cms_BF6_splitup_7_IVAexport.txt")
v8 = read.delim(file = "variant_list/cms_BF6_splitup_8_IVAexport.txt")

iva_list = rbind.data.frame(v1,v2,v3,v4,v5,v6,v7,v8)
names(iva_list) = tolower(names(iva_list))
iva_list = iva_list %>%
  mutate(variant = paste0(chromosome,"_",position))
rm(v1,v2,v3,v4,v5,v6,v7,v8)

vlist = vlist %>%
  left_join(vlist,iva_list,by = "variant")
