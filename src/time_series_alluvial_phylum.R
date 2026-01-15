library(tidyverse)
library(vegan)
library(phyloseq)
library(ggalluvial)
library(scales)

setwd('~/WorkForaging/Academia/Nicole/deca/')

ps <- readRDS('data/16Sdada2/phyloseq.filtered.rds')
meta <- sample_data(ps) %>% data.frame() %>% rename(Sample=SampleID)
timelevs <- c("Jul_2022","Nov_2022","Feb_2023", "Aug_2023","Nov_2023", "Feb_2024","Aug_2024")


ps_phylum <- ps %>%
  tax_glom(taxrank = "Phylum", NArm = TRUE) %>%
  transform_sample_counts(function(x) x / sum(x))

psdat <- psmelt(ps_phylum) %>%
  as_tibble() %>%
  transmute(
    Sample = as.character(Sample),
    Phylum = as.character(Phylum),
    Abundance = as.numeric(Abundance)
  ) 

thresh <- 0.01 #threshold for 'other' category
psdat <- psdat %>% 
  mutate(Phylum=case_when(Abundance<thresh ~ "<1% Abundance", 
                          .default = psdat$Phylum)) %>% 
  group_by(Sample, Phylum) %>% 
  summarize(Abundance=sum(Abundance)) %>% 
  left_join(., meta, by='Sample')

phy_order <- psdat %>%
  group_by(Phylum) %>%
  summarize(mean_abund = mean(Abundance), .groups = "drop") %>%
  arrange(-mean_abund) %>%
  pull(Phylum)

gdat <- psdat %>% 
  filter(Surface_Subsurface=='S') %>% 
  ungroup() %>% group_by(SampleDate, Phylum) %>% 
  summarize(Abundance=sum(Abundance)) %>% 
  mutate(Abundance = Abundance/sum(Abundance),
         SampleDate = factor(SampleDate, levels=timelevs), 
         Phylum = factor(Phylum, levels=phy_order))

plevs <- unique(gdat$Phylum)
palette <- qualitative_hcl(length(plevs), palette = "Dark 2") # or "Set3", "Classic Dark"
palette <- setNames(palette, nm = rev(plevs) )

ggplot(gdat)+
  geom_alluvium(aes(x=SampleDate, y=Abundance, alluvium=Phylum, fill=Phylum))+
  geom_stratum(aes(x=SampleDate, y=Abundance, stratum=Phylum, fill=Phylum))+
  theme_bw() +
  scale_fill_manual(values=palette)+
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Surface Microbiome Time Stratified")



gdat <- psdat %>% 
  filter(Surface_Subsurface=='SS') %>% 
  ungroup() %>% group_by(SampleDate, Phylum) %>% 
  summarize(Abundance=sum(Abundance)) %>% 
  mutate(Abundance = Abundance/sum(Abundance),
         SampleDate = factor(SampleDate, levels=timelevs), 
         Phylum = factor(Phylum, levels=phy_order))

ggplot(gdat)+
  geom_alluvium(aes(x=SampleDate, y=Abundance, alluvium=Phylum, fill=Phylum))+
  geom_stratum(aes(x=SampleDate, y=Abundance, stratum=Phylum, fill=Phylum))+
  theme_bw() +
  scale_fill_manual(values=palette)+
  scale_y_continuous(expand = c(0,0)) +
  ggtitle("Subsurface Microbiome Time Stratified")


