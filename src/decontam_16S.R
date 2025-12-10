
library(tidyverse)
library(vegan)
library(phyloseq)
library(decontam)

setwd('WorkForaging/Academia/Nicole/deca/')

ps <- readRDS('data/16Sdada2/phyloseq_dada2.rds')

#check which is a neg control
sample_data(ps)[which(grepl('Neg', sample_data(ps)$Microhabitat)),]
sample_data(ps)$is.neg <- sample_data(ps)$Microhabitat=='Neg_Control'

#run decontam
decdf <- isContaminant(ps, method='prevalence', neg='is.neg')
cindx <- which(decdf$contaminant)

#check abundances
gdat <- as.data.frame(otu_table(ps))[,cindx]
gdat[which(gdat>0)]

