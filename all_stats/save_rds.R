library(dplyr)

v1 = read.csv("vlist1.csv")
v2 = read.csv("vlist2.csv")
v3 = read.csv("vlist3.csv")
v4 = read.csv("vlist4.csv")
v5 = read.csv("vlist5.csv")
v6 = read.csv("vlist6.csv")
v7 = read.csv("vlist7.csv")
v8 = read.csv("vlist8.csv")

vlist = rbind.data.frame(v1,v2,v3,v4,v5,v6,v7,v8)
save(vlist, file = "cdp_cms_gwas_regression_stats.RData")