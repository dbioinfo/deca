library(tidyverse)
library(phyloseq)
library(ape)
library(Biostrings)
library(picante)
library(dada2)
setwd("~/WorkForaging/Academia/Nicole/deca/")

##write asv sequences to file
ps <- readRDS('data/16Sdada2/phyloseq.filtered.rds')
refs <- refseq(ps)
writeXStringSet(refs, "data/16Sqiime/asv-seqs.fasta", format = "fasta")

##run gen newick on cluster

##add phylo tree into object
ptree <- read.tree('data/16Sqiime/asv-seqs.newick/tree.nwk')

asvs <- taxa_names(ps)
asvs <- intersect(asvs, ptree$tip.label)

ptree <- keep.tip(ptree, asvs)
ptree <- collapse.singles(ptree)

ps <- merge_phyloseq(ps, phy_tree(ptree))

##annotate taxonomy with cyanoseq, greengenes, and silva
iseqs <- as.character(refs)
tax <- assignTaxonomy(iseqs, "raw_data/cyanoseq/CyanoSeqV1.3_GSRDB_dada2.fastq.gz", multithread=TRUE, tryRC = T)
#tax <- addSpecies(tax, "raw_data/silva_ref/silva_v138.2_assignSpecies.fa.gz")
rownames(tax) <- names(iseqs)

#identify original ASVs labelled Phylum Cyanobacteriota
tt1 <- data.frame(tax_table(ps))
cyanoasvs <- tt1 %>% filter(Phylum=="Cyanobacteriota") %>% rownames()

#figure out taxonomic depth of each annotation
tt2 <- data.frame(tax_table(tax))
tt1$ASV <- rownames(tt1)
tt2$ASV <- rownames(tt2)
tt1 <- tt1 %>% select(Kingdom,Phylum,Class,Order,Family,Genus,ASV)
tt2 <- tt2 %>% select(Kingdom,Phylum,Class,Order,Family,Genus,ASV)
tt1count <- tt1 %>% 
  group_by(ASV) %>% 
  summarise(across(everything(), ~ any(is.na(.)))) %>% 
  rowwise(ASV) %>% 
  summarise(na_tt1 = sum(c_across(everything())))
tt2count <- tt2 %>% 
  group_by(ASV) %>% 
  summarise(across(everything(), ~ any(is.na(.)))) %>% #mark cols as NA (T) or not (F)
  rowwise(ASV) %>% 
  summarise(na_tt2 = sum(c_across(everything()))) #these are the number of NA cols
ttcount <- merge(tt1count, tt2count, by='ASV')

ttcount %>% filter(ASV %in% cyanoasvs, na_tt1<na_tt2) %>% nrow() 


#tax comparison
rseq <- data.frame(refseq(ps))
rseq <- rseq[cyanoasvs, ]
rseq <- data.frame(ASV=cyanoasvs, sequence=rseq)

comp <- merge(tt1[cyanoasvs, ], tt2[cyanoasvs,], by="ASV")
final <- merge(comp, rseq, by="ASV")
write_csv(final, 'data/two_tax_comparison.w.seq.csv')

tax_table(ps) <- tax_table(tax)

saveRDS(ps, 'data/16Sdada2/phyloseq.filtered.tree.rds')

##calculate faith alpha diversity
pd_res <- pd(otu_table(ps), phy_tree(ps), include.root = F)
