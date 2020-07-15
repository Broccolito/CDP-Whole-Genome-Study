library(here)

load("cdp_epas1.rds")
met_sig = read.csv("Metabolic Significant Pairs.csv")
hvr_sig = read.csv("HVR Significant Pairs.csv")

plot_gene = function(sig_mat = met_sig, phenotype = "Hb", color = "red"){
  
}

SNPs = colnames(cdp)
gene_line = rep(0, length(SNPs))

plot_gene = function(significance_mat = met_sig, pheno = "Hb", color = "red",
                     position = FALSE){
  plot(gene_line, type = "l", lwd = 5, xaxt = "n", yaxt = "n",
       xlab = ifelse(position, "Position", ""), ylab = pheno)
  
  find_snp = function(x){
    which(x == SNPs)
  }
  snp_index = unique(unlist(sapply(subset(significance_mat, phenotype == pheno)$SNP,
                                   function(x){which(x == SNPs)})))
  points(snp_index, rep(0, length(snp_index)), col = color, cex = 1.5, pch = 16)
}

windows()
par(mfrow = c(4,1), mai = c(0.1,0.8,0.1,0.1))
plot_gene(pheno = "Hct")
plot_gene(pheno = "CMS_SCORE")
plot_gene(pheno = "Hb")
plot_gene(pheno = "SpO2")

windows()
par(mfrow = c(4,1), mai = c(0.1,0.8,0.1,0.1))
plot_gene(pheno = "BMI")
plot_gene(pheno = "Insulin")
plot_gene(pheno = "Glucose")
plot_gene(pheno = "Iron")

windows()
par(mfrow = c(4,1), mai = c(0.1,0.8,0.1,0.1))
plot_gene(pheno = "bmi", significance_mat = hvr_sig)
plot_gene(pheno = "odi_b", significance_mat = hvr_sig)
plot_gene(pheno = "satbelow80_b", significance_mat = hvr_sig)
plot_gene(pheno = "hvr", significance_mat = hvr_sig)


