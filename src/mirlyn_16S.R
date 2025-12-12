library(tidyverse)
library(dada2)
library(phyloseq)
library(Biostrings)
library(colorspace)
library(vegan)
library(mirlyn)

#import data
ps <- t(readRDS('data/16Sdada2/phyloseq_dada2.rds'))
otu <- otu_table(ps, taxa_are_rows = F)
meta <- data.frame(sample_data(ps))
psdf <- phyloseq_to_df(t(ps))

#check out abundance plot to be sure
#ibar <- bartax(ps, "Sample", taxrank="Phylum")
#ibar

#rare curve -- cpu and ram are determined by *both* rep and mc.cores kwargs, watch RAM and kill if close to maxing out
whole_rep <- rarefy_whole_rep(ps, rep = 20, mc.cores=80)
png('figs/mir_rarecurve.png', width=1500,height=1500, units='px')
rarecurve(whole_rep, sample = "Sample")
dev.off()

#repeated subsampling rarefaction with mirl
mrare <- mirl(ps, libsize=10000, rep=50, replace=F, mc.cores=80)
saveRDS(mrare, 'data/mirlyn/mirl.rds')
#q()
#mrare <- readRDS('data/mirlyn/mirl.rds')

#lets check the alpha diversity metrics first
alphadiv_df <- alphadivDF(mrare,diversity = "shannon")
galpha <- alphawichVis(alphadiv_df, xvar = "SampleID", colorvar = "Microhabitat") + 
  theme(axis.text.x = element_text(angle=75, vjust = 0.7, hjust = .6))
png('figs/mirl_alpha_10k.png',height=1500,width=1500,units='px')
galpha
dev.off()

acone <- alphacone(ps, rep=20, diversity = "shannon", replace=F, mc.cores=80)
gcone <- alphaconeVis(acone, cols = "Microhabitat")
png('figs/mirl_alphacone_10k.png',height=1500,width=1500,units='px')
gcone
dev.off()
q()

#mirl beta diversity
betamatPCA_object <- betamatPCA(mrare, dsim = "bray", transformation = "hellinger")
betaPCA <- betamatPCA_object$x[,1:5]
betaPCA <- data.frame(PC1=betaPCA[,"PC1"], PC2=betaPCA[,"PC2"], PC3=betaPCA[,"PC3"], Sample=rownames(betaPCA))
betaPCA <- betaPCA %>% separate(Sample, into=c("rep","SampleID"), sep = '-')
betaPCA <- betaPCA %>% left_join(., meta, by="SampleID")
ggplot(betaPCA)+
  geom_point(aes(x=PC2, y=PC3, color=Shrub))+
  ggtitle("Bray-Curtis Beta Diversity PCA Rarefied")+
  theme_bw()

#bespoke beta diversity
#mirl only does PCA of beta diversity, so we'll need to improvise 
#their function just uses vegan under the hood anyways

btrans <- decostand(as.matrix(t(repotu_df(ps))), transformation="hellinger")

