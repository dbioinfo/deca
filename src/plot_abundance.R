library(tidyverse)
library(dada2)
library(phyloseq)
library(Biostrings)
library(colorspace)
library(vegan)
setwd('~/WorkForaging/Academia/Nicole/deca')


otu <- read_csv('data/16Sdada2/dada2_otu.csv')
ab <- t(otu[,2:length(otu)])
colnames(ab) <- otu$Sample
rownames(ab) <- colnames(otu)[2:length(otu)]
otus <- colnames(otu)[2:length(otu)]

ps <- readRDS('data/16Sdada2/phyloseq_dada2.rds')
meta <- as_tibble(sample_data(ps)) %>% 
  mutate(Microhabitat = case_when(
    Microhabitat=='ASI'~"SI", 
    Microhabitat=='ASSG'~"SSG", 
    Microhabitat=='BSI'~"SI", 
    Microhabitat=='BSSG'~"SSG", 
    Microhabitat=='ASG'~"SG", 
    Microhabitat=='BSG'~"SG", 
    .default=meta$Microhabitat))
tax <- tax_table(ps)
asv <- rownames(tax)
tax <- as.data.frame(as.matrix(tax))
tax$OTU <- asv


tb <- as_tibble(ab) %>% 
  mutate(OTU=otus) %>% 
  pivot_longer(-OTU, names_to = "SampleID", values_to = "count") %>% 
  filter(count>0) %>% 
  left_join(., tax %>% select(OTU, Phylum), by="OTU") %>% 
  group_by(SampleID, Phylum) %>% 
  summarize(count=sum(count)) %>%
  group_by(SampleID) %>%  
  mutate(tcount=sum(count)) %>% 
  group_by(Phylum) %>% 
  mutate(percent=count*100/tcount) %>% 
  arrange(-percent) 
tb <- tb %>% mutate(Phylum=case_when(percent<1~"<1% Abundance",
                                     .default = Phylum)) %>% 
  ungroup() %>% 
  group_by(SampleID, Phylum) %>% 
  summarize(count=sum(count)) %>%
  ungroup() %>% 
  group_by(SampleID) %>% 
  mutate(tcount=sum(count)) %>% 
  group_by(Phylum) %>%   
  mutate(percent=count*100/tcount) %>% 
  arrange(-percent) %>% 
  left_join(meta %>% select(-c(Notes, ...10, ...11)), by='SampleID') 


plevs <- unique(tb$Phylum)
palette <- qualitative_hcl(length(plevs), palette = "Dynamic") # or "Set3", "Classic Dark"
palette <- setNames(palette, nm = plevs)

for (timepoint in unique(meta$SampleDate)){
  gg <- ggplot(tb %>% filter(SampleDate==timepoint))+
    geom_bar(aes(x=SampleID, y=percent, fill=Phylum), stat='identity', na.rm = T) +
    facet_wrap(~Treatment, nrow=1, scales = "free_x",
               labeller = labeller(label_wrap_gen(width = 50)))+
    scale_fill_manual(values = palette)+
    scale_y_continuous(expand=c(0,0))+
    theme_bw()+
    theme(text = element_text(size=12),
          axis.text.x = element_text(angle =60, vjust = 1, hjust = 1, size=5))+
    xlab("Sample")+
    ylab("Relative Abundance")+
    ggtitle(paste0("Relative abundance across treatment in ", timepoint))
  ggsave(paste0('figs/rel_abundance_treatment_',timepoint,'.png'),gg)
}
