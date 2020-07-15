suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(gridExtra)
})

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
    })) %>%
  na.omit()
  
  pheno_max = max(gp_df$pheno) + 0.25 * diff(range(gp_df$pheno))
  pheno_min = min(gp_df$pheno) - 0.1 * diff(range(gp_df$pheno))
  
  variant_pos = paste0("Chr",paste(unlist(strsplit(variant,"_")),collapse = ":"))
  sub_title = paste0( gene_name, "; ",variant_pos)
  
  if(is.null(phenotype_name)){
    phenotype_name = phenotype
  }
  
  suppressWarnings({
    label_low = max(gp_df$pheno)+0.1*diff(range(gp_df$pheno))
    label_high = max(gp_df$pheno)+0.2*diff(range(gp_df$pheno))
    
    plt = ggplot(data = gp_df,aes(x = geno,y = pheno)) + 
      geom_boxplot(outlier.color = "white") +
      geom_point(size = 2, position = position_dodge2(0.3)) + 
      xlab("Copy of Adaptive Allele") + 
      ylab(phenotype_name) + 
      ylim(c(pheno_min,pheno_max)) + 
      facet_grid(.~sex) + 
      ggtitle(label = NULL,subtitle = sub_title) +
      stat_compare_means(comparisons = list(c("0","1")),method = "t.test",paired = FALSE,label.y = label_low) +
      stat_compare_means(comparisons = list(c("1","2")),method = "t.test",paired = FALSE,label.y = label_low) +
      stat_compare_means(comparisons = list(c("0","2")),method = "t.test",paired = FALSE,label.y = label_high) +
      theme_bw() + 
      theme(text = element_text(size = 15))
  return(plt)
  })
}

cat("Boxplot Rendering Engine started...\n")

# Example
# plot_pheno_geno(phenotype_name = "%Hct")
