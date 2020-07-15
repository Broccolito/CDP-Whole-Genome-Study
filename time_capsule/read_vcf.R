library(vcfR)
library(adegenet)
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(ggmap)

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

vcf_wd = list.files(pattern = ".vcf.gz")

vcf = read.vcfR(vcf_wd, verbose = TRUE)

aa.genlight = vcfR2genlight(vcf, n.cores = 6)
locNames(aa.genlight) = paste(vcf@fix[,1],vcf@fix[,2],sep="_") # add
# pop(aa.genlight) <- c(rep("CDP", 36), rep("Tibetan", 27))

cdp = as.data.frame(aa.genlight)

rm(vcf)
gc()

cdp = cbind.data.frame(
  data.frame(wgs_id = sapply(rownames(cdp), function(x){unlist(strsplit(x, "_"))[1]})),
  cdp)
rownames(cdp) = NULL

save(cdp, file = "cdp.rds")

beepr::beep()
