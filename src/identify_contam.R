library(tidyverse)
library(vegan)
library(phyloseq)

setwd('WorkForaging/Academia/Nicole/deca/')
ps <- readRDS('data/16Sdada2/phyloseq_dada2.rds')
otu <- read_csv('data/16Sdada2/dada2_otu.csv')

#identify species-level assignments of each control sample
tax <-as.data.frame(as.matrix(tax_table(ps)))

#first, filter the table to include only control samples
subs <- otu %>% filter(grepl('STANDARD|C6',Sample) )

#get index of ASVs > 0
sidx <- which(colSums(subs[2:length(subs)])>0) +1

#abundance plots of each 
tax <- tax %>% unite(full_name, Phylum, Genus, sep = '_')
tax$ASV <- rownames(tax)
abund <- subs[,c('Sample',names(sidx) )]  %>% pivot_longer(-Sample, names_to = 'ASV', values_to = 'Abundance')
abund <- left_join(abund, tax, by = 'ASV')
abund <- abund %>% group_by(Sample) %>% mutate(Rel_Abundance = Abundance / sum(Abundance))


pal <- colorRampPalette(c("dodgerblue", "darkgreen", "firebrick"))(20)
pal <- sample(pal)

ggplot(abund )+ #%>% filter(grepl('C6',Sample))
  geom_bar(aes(x=Sample, y=Abundance, fill=full_name), stat='identity') +
  theme_bw()+
  scale_fill_manual( values = pal)
