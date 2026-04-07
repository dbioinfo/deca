library(tidyverse)
library(phyloseq)

setwd('~/WorkForaging/Academia/Nicole/deca/')

ps <- readRDS("data/16Sdada2/phyloseq.filtered.tree.rds")
tax <- data.frame(tax_table(ps))
cyanoasv <- tax %>% filter(grepl("Cyanobac",Phylum)) %>% rownames()
ps <- subset_taxa(ps, taxa_names(ps) %in% cyanoasv)
ps <- subset_samples(ps, Surface_Subsurface=="S")
ps <- subset_samples(ps, sample_sums(ps)>0)
saveRDS(ps, "data/cyano/phyloseq.cyano.tree.rds")


