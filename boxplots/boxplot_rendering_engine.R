library(dplyr)
library(ggplot2)
library(ggpubr)

cat("Loading Genome Files...\n")
load("cdp.rds")
cdp2015 = read.csv("2015_aug_dec_2016_dec.csv")
wgs2018 = read.csv("2018_wgs.csv")

cat("Fusing Genome Files with Phenotypical Files...\n")
suppressWarnings({
  wgs = cdp[,1,drop = FALSE] %>%
    left_join(wgs2018, by = "wgs_id") %>%
    dplyr::select(wgs_id,id) %>%
    left_join(cdp2015, by = "id") %>%
    dplyr::select(wgs_id,id,initials,sex:abpm_conventional_dbp) %>%
    mutate_if(is.numeric, as.numeric) %>%
    select(-wgs_id) %>%
    cbind.data.frame(cdp)
})

rm(cdp)

plot_pheno_geno = function(gene_name = "F5",
                           variant = "1_169519112",
                           phenotype = "hct_vena_mean",
                           phenotype_name = NULL,
                           res = 1200,
                           w = 4,
                           h = 6){
  sex = wgs[["sex"]]
  geno = wgs[[variant]]
  pheno = as.numeric(wgs[[phenotype]])
  
  gp_df = data.frame(sex,geno,pheno) %>%
    mutate(geno = as.factor(geno)) %>%
    mutate(sex = sapply(sex, function(x){
      return(ifelse(x=="M","Male","Female"))
    }))
  
  variant_pos = paste0("chr",paste(unlist(strsplit(variant,"_")),collapse = ":"))
  sub_title = paste0( gene_name, "; ",variant_pos)
  
  if(is.null(phenotype_name)){
    phenotype_name = phenotype
  }
  
  suppressWarnings({
    if(length(unique(gp_df$geno))==2){
      comparisons = list( c("0", "1"))
    }else{
      comparisons = list(c("0","1"),c("1","2"),c("0", "2"))
    }
    
    ggplot(data = gp_df,aes(x = geno,y = pheno)) + 
      geom_boxplot(outlier.shape = 3) +
      facet_grid(rows = .~sex) + 
      geom_point(size = 2, position = position_dodge2(0.3)) + 
      xlab("Copy of Adaptive Allele") + 
      ylab(phenotype_name) + 
      ggtitle(label = NULL,subtitle = sub_title) +
      stat_compare_means(comparisons = comparisons) +
      theme_bw() + 
      ggsave(paste0("boxplots_result/",phenotype," Vs. ",variant,".png"),device = "png",
             dpi = res,width = w,height = h)
  })
  
}

cat("Boxplot Rendering Engine started...\n")

# Example
# plot_pheno_geno(phenotype_name = "%Hct")
