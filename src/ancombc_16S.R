library(tidyverse)
library(colorspace)
library(phyloseq)
library(ANCOMBC)
setwd("~/WorkForaging/Academia/Nicole/deca/")

#import data
ps <- readRDS('data/16Sdada2/phyloseq.filtered.rds')
meta <- data.frame(sample_data(ps))

timelevs <- c("Jul_2022","Nov_2022","Feb_2023", "Aug_2023","Nov_2023", "Feb_2024","Aug_2024")
treatlevs <- c("Baseline","Control","Water","N", "N+P","N+Water","P","P+Water","N+P+Water")
mhablevs <- c("SI","SG","SSI","SSG")
gilevs <- c("I","G")
sslevs <- c("S","SS")

#fix order
meta <- meta %>% 
  mutate(Microhabitat = factor(Microhabitat, levels=mhablevs) ) %>% 
  mutate(Grass_Interspace = factor(Grass_Interspace, levels=gilevs) ) %>% 
  mutate(Surface_Subsurface = factor(Surface_Subsurface, levels=sslevs) ) %>% 
  mutate(Treatment = factor(Treatment, levels=treatlevs) ) %>% 
  mutate(SampleDate = factor(SampleDate, levels=timelevs) )

sample_data(ps) <- sample_data(meta)

#optional remove baseline samples
ps <- subset_samples(ps, SampleDate!="Jul_2022")
sample_data(ps)$SampleDate <- droplevels(sample_data(ps)$SampleDate)

#optional rename baseline treatment to control (choose this one)
sample_data(ps) <- data.frame(sample_data(ps)) %>%   
  mutate(Treatment=case_when(SampleDate=="Jul_2022"~"Control",
    .default=sample_data(ps)$Treatment))

#run 
out <- ancombc2(data = ps,
               fix_formula='Surface_Subsurface+Grass_Interspace+SampleDate+Treatment',
               rand_formula = NULL,
               tax_level="Phylum",
               p_adj_method = 'holm',
               prv_cut = 0.01,
               n_cl = 10,
               group='Treatment',
               pairwise = T,
               struc_zero = T,
               neg_lb = F)

#save
saveRDS(out, 'data/16Sancombc/ancombc_results_treat_phylum.rds')

