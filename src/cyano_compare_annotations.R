library(tidyverse)
library(dada2)
setwd('~/WorkForaging/Academia/Nicole/deca/')
ps <- readRDS('data/cyano/phyloseq.cyano.tree.rds')

seqs <- refseq(ps)
ref_df <- data.frame(ASV = names(seqs), seq = as.character(seqs))

tax <- assignTaxonomy(ref_df$seq, "../drhizo/raw_data/silva_ref/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)
tax <- addSpecies(tax, "../drhizo/raw_data/silva_ref/silva_v138.2_assignSpecies.fa.gz")

ctax <- assignTaxonomy(ref_df$seq, 'raw_data/cyanoseq/CyanoSeqV1.3_GSRDB_dada2.fastq.gz')
ctax <- data.frame(tax_table(ctax)) %>% rownames_to_column("seq")

tmp <- data.frame(tax_table(tax)) %>% rownames_to_column("seq")
tmp <- left_join(tmp, ref_df, by='seq')
tmp <- left_join(tmp, ctax, by='seq')
colnames(tmp) <- sub('\\.y','.CYANOSEQ',sub('\\.x','.SILVA',colnames(tmp)))
tmp <- tmp %>% relocate(ASV, seq)

write_csv(tmp,'data/cyano/taxonomy_comparison.csv')
