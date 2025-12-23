library(tidyverse)
library(colorspace)
library(phyloseq)
library(ANCOMBC)
setwd("~/WorkForaging/Academia/Nicole/deca/")

#import data
ps <- readRDS('data/16Sdada2/phyloseq.filtered.rds')
meta <- sample_data(ps)

#run with Microhabitat+SampleDate to start
out <- ancombc2(data = ps,
               fix_formula='Microhabitat+SampleDate',
               rand_formula = NULL,
               tax_level="Genus",
               p_adj_method = 'holm',
               prv_cut = 0.1,
               n_cl = 20,
               group='SampleDate',
               pairwise = T,
               struc_zero = T,
               neg_lb = F)

#save
saveRDS('data/16Sancombc/ancombc.out.rds')

#graph results
out$res_pair %>% select(starts_with("diff")) %>% summarize_all(sum)
gdat <- out$res_pair
sigtaxs<-gdat %>% filter(diff_ShrubCreosote|diff_ShrubMariola|diff_ShrubMesquite|diff_ShrubTarbush) %>% pull(taxon)
ggdat <- gdat %>% 
  filter(taxon %in% sigtaxs) %>% 
  group_by(taxon) %>% 
  select(starts_with("lfc")) %>% 
  pivot_longer(starts_with("lfc"),names_to='comparison', values_to = 'LFC') %>% 
  mutate(comparison = sub("^lfc_", "", comparison)) %>% 
  mutate(comparison= sub("*Shrub*","", comparison))%>% 
  mutate(comparison= sub("*Shrub*","", comparison))

glevs <- c("Creosote","Mariola", "Mesquite","Tarbush", "Mariola_Creosote", 
           "Mesquite_Creosote", "Tarbush_Creosote",  "Mesquite_Mariola","Tarbush_Mariola","Tarbush_Mesquite" )

ggdat <- ggdat %>% mutate(comparison=factor(comparison, levels=glevs))

gg <- ggplot(ggdat)+
  geom_tile(mapping=aes(x=comparison, y=taxon, fill=LFC))+
  scale_fill_viridis_c(option='plasma')+
  theme_bw() +
  theme(
    text = element_text(size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(), 
    panel.border = element_blank())+
  xlab("")+
  ylab("Taxon")
gg
ggsave('figs/Genus_lvl_heatmap.png', height=4000, width=4000, dpi=300, units='px')
