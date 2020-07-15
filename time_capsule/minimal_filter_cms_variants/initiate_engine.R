library(dplyr)
library(rstudioapi)

rm(list = ls())
gc()

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

jobRunScript(path = "run_iva_apriori.R", importEnv = TRUE, workingDir = wd, name = "IVA Apriori")
jobRunScript(path = "run_iva_exon.R", importEnv = TRUE, workingDir = wd, name = "IVA Exon")
jobRunScript(path = "run_iva_minimal_sub1.R", importEnv = TRUE, workingDir = wd, name = "IVA Minimal Sub1")
jobRunScript(path = "run_iva_minimal_sub2.R", importEnv = TRUE, workingDir = wd, name = "IVA Minimal Sub2")
jobRunScript(path = "run_iva_minimal_sub3.R", importEnv = TRUE, workingDir = wd, name = "IVA Minimal Sub3")
jobRunScript(path = "run_iva_minimal_sub4.R", importEnv = TRUE, workingDir = wd, name = "IVA Minimal Sub4")
jobRunScript(path = "run_iva_minimal_sub5.R", importEnv = TRUE, workingDir = wd, name = "IVA Minimal Sub5")
jobRunScript(path = "run_iva_minimal_sub6.R", importEnv = TRUE, workingDir = wd, name = "IVA Minimal Sub6")
