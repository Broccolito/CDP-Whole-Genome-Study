library(vcfR)
library(adegenet)
library(adegraphics)
library(pegas)
library(StAMPP)
library(lattice)
library(gplots)
library(ape)
library(ggmap) 

vcf_wd = list.files(pattern = "vcf.gz")

vcf <- read.vcfR(vcf_wd, verbose = FALSE )

aa.genlight <- vcfR2genlight(vcf, n.cores = 1)
locNames(aa.genlight) <- paste(vcf@fix[,1],vcf@fix[,2],sep="_") # add
pop(aa.genlight) <- c(rep("CDP", 36), rep("Tibetan", 27))

cdp = as.data.frame(aa.genlight)[1:36,]

rownames(cdp) = sapply(rownames(cdp), function(x){unlist(strsplit(x, "_"))[1]})

save(cdp, file = "cdp_epas1.rds")
