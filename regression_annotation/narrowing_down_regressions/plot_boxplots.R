rm(list = ls())
gc()
graphics.off()

source("boxplot_rendering_engine.R")

results_apriori = read.csv("annotated_results/results_apriori.csv")
results_exonic = read.csv("annotated_results/results_exonic.csv")
results_cms_gt6_lit_search = read.csv("annotated_results/results_cms_gt6_lit_search.csv")
results_cms_gt8.1 = read.csv("annotated_results/results_cms_gt8.1.csv")
results_cms_gt8.2 = read.csv("annotated_results/results_cms_gt8.2.csv")

make_boxplots = function(results_list,list_name){
  results_list = results_list %>%
    mutate_if(is.factor,as.character)
  for(i in 1:dim(results_list)[1]){
    try({
      plt = plot_pheno_geno(gene_name = results_list$gene_symbol[i],
                            phenotype = results_list$phenotype[i],
                            variant = results_list$variant[i],
                            phenotype_name = "%Hct")
      ggsave(filename = paste0("boxplots/",list_name," ",
                               results_list$phenotype[i]," Vs. ", results_list$variant[i],".png"),
             device = "png",dpi = 1200,width = 4.5,height = 6,plot = plt)
    },silent = TRUE)
  }
}

cat("Plotting in Progress...\n")
make_boxplots(results_apriori,list_name = "apriori")
cat("Plotting in Progress...\n")
make_boxplots(results_exonic,list_name = "exonic")
cat("Plotting in Progress...\n")
make_boxplots(filter(results_cms_gt6_lit_search,pm<=0.001),list_name = "cms_gt6_lit_search")
cat("Plotting in Progress...\n")
make_boxplots(filter(results_cms_gt8.1,pm<=0.001),list_name = "cms_gt8.1")
cat("Plotting in Progress...\n")
make_boxplots(filter(results_cms_gt8.2,pm<=0.001),list_name = "cms_gt8.2")
cat("Done!")