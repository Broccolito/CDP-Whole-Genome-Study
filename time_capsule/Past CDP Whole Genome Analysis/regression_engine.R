rm(list = ls())

library(here)
library(beepr)

load("cdp_epas1.rds")
cdp_epas1 = cdp
rm(cdp)
matching_file = read.csv("CDP36_matching.csv")
hvr = read.csv("HVR.csv")
met = read.csv("metabolic.csv")

matching_file = matching_file[(matching_file$WGS_ID %in% 
                                 intersect(rownames(cdp_epas1), matching_file$WGS_ID)),]

matching_file = matching_file[order(matching_file$WGS_ID),]
cdp_epas1 = cdp_epas1[order(rownames(cdp_epas1)),]

matching_file = cbind.data.frame(matching_file, cdp_epas1)
rownames(matching_file) = NULL
rm(cdp_epas1)
matching_file = matching_file[order(matching_file$ID),]

hvr = hvr[hvr$id %in% intersect(matching_file$ID, hvr$id),]
hvr_matching_file = matching_file[matching_file$ID %in% intersect(matching_file$ID, hvr$id),]
hvr = hvr[order(hvr$id),]

met = met[met$id %in% intersect(matching_file$ID, met$id),]
met_matching_file = matching_file[matching_file$ID %in% intersect(matching_file$ID, met$id),]
met = met[order(met$id),][-4,] # Remove the duplicate with more NAs

met_significiant_pairs = vector()
for(v in names(met)){
  if(v != "id" & v != "Initials" & v != "sex" & v != "age"){
    for(SNP in names(met_matching_file)){
      if(SNP != "WGS_ID" & SNP != "ID" & SNP != "INITIALS"){
        try({
          xy = data.frame(x = met[[v]], y = met_matching_file[[SNP]])
          xy = na.omit(xy)
          n = dim(xy)[1]
          pvalue = summary(lm(data = xy, y~x))$coefficient[2,4]
          if(pvalue <= 0.05){
            pairing_found = c(v, SNP, n, pvalue)
            met_significiant_pairs = rbind(met_significiant_pairs, pairing_found)
            # beep()
          }
          cat(paste0(v, " Vs. ", SNP, "\n"))
        }, silent = TRUE)
      }
    }
  }
}

hvr_significiant_pairs = vector()
for(v in names(hvr)){
  if(v != "id" & v != "Initials" & v != "sex" & v != "age"){
    for(SNP in names(hvr_matching_file)){
      if(SNP != "WGS_ID" & SNP != "ID" & SNP != "INITIALS"){
        try({
          xy = data.frame(x = hvr[[v]], y = hvr_matching_file[[SNP]])
          xy = na.omit(xy)
          n = dim(xy)[1]
          pvalue = summary(lm(data = xy, y~x))$coefficient[2,4]
          if(pvalue <= 0.05){
            pairing_found = c(v, SNP, n, pvalue)
            hvr_significiant_pairs = rbind(hvr_significiant_pairs, pairing_found)
            # beep()
          }
          cat(paste0(v, " Vs. ", SNP, "\n"))
        }, silent = TRUE)
      }
    }
  }
}

met_significiant_pairs = data.frame(phenotype = met_significiant_pairs[,1],
                                    SNP = met_significiant_pairs[,2],
                                    n = met_significiant_pairs[,3],
                                    pvalue = met_significiant_pairs[,4])
hvr_significiant_pairs = data.frame(phenotype = hvr_significiant_pairs[,1],
                                    SNP = hvr_significiant_pairs[,2],
                                    n = hvr_significiant_pairs[,3],
                                    pvalue = hvr_significiant_pairs[,4])

write.csv(met_significiant_pairs, file = "Metabolic Significant Pairs.csv",
            quote = FALSE, row.names = FALSE)

write.csv(hvr_significiant_pairs, file = "HVR Significant Pairs.csv",
            quote = FALSE, row.names = FALSE)