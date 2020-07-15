library(ggplot2)
library(dplyr)
library(ggsci)

rm(list = ls())
gc()
graphics.off()

cms_gt8 = rbind.data.frame(
  read.csv(file = "annotated_results/results_cms_gt8.1.csv"),
  read.csv(file = "annotated_results/results_cms_gt8.2.csv")
)

cms_gt8 = cms_gt8 %>%
  mutate(nlogPm = -log(pm,base = 10)) %>%
  mutate_if(is.factor,as.character) %>%
  mutate(chr = sapply(variant,function(x){
    unlist(strsplit(x,"_"))[1] %>%
      as.numeric()
  })) %>%
  mutate(loc = sapply(variant,function(x){
    unlist(strsplit(x,"_"))[2] %>%
      as.numeric()
  })) %>%
  arrange(chr,loc) %>%
  mutate(chr = as.factor(chr)) %>%
  filter(pm<=0.05) %>%
  mutate(ann = ifelse(nlogPm<4,"",gene_symbol))

cms_gt8 = mutate(cms_gt8, index = 1:dim(cms_gt8)[1]) 
  
mp_palette =  c("#9999CC", "gray13", "gold2", "plum4", "darkorange1",
                "lightskyblue2", "firebrick", "burlywood3", "gray51",
                "springgreen4", "lightpink2", "deepskyblue4", "lightsalmon2",
                "mediumpurple4", "orange", "maroon", "yellow3",
                "brown4", "yellow4", "sienna4", "chocolate", "gray19")

plt = ggplot(data = cms_gt8,aes(x = index,y = nlogPm)) + 
  geom_point(aes(y = nlogPm-0.2,color = chr)) + 
  geom_text(aes(label = ann)) + 
  scale_color_manual(values = mp_palette) +
  ylab("-Log(P)") + 
  xlab(NULL) + 
  ggtitle(label = NULL, subtitle = "CMS-GWAS %Hct in Males") + 
  theme_pubclean() + 
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none",
        text = element_text(size = 15)) + 
  ggsave(filename = "manhattan_plots/cms_gwas_hct_males.png",device = "png",width = 12,height = 7,
         dpi = 1200)

# scale_fill_brewer(palette = "Dark2")
# theme(legend.title = element_blank())

