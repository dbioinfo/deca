library(tidyverse)
library(vegan)
library(phyloseq)
library(ggalluvial)
library(scales)

setwd('~/WorkForaging/Academia/Nicole/deca/')

#import data
ps <- readRDS('data/16Sdada2/phyloseq.filtered.rds')
otu <- data.frame(otu_table(ps))
meta <- data.frame(sample_data(ps))
timelevs <- c("Jul_2022","Nov_2022","Feb_2023", "Aug_2023","Nov_2023", "Feb_2024","Aug_2024")
treatlevs <- c("Baseline","Control","Water","N", "N+P","N+Water","P","P+Water","N+P+Water")
treatcols <- c("#E3E571","#E5AE71","#71C4E5","#53DF7B","#C0E5A4","#BCE7D6","#DF53B7","#E7BCCC","#804DEF")
mhablevs <- c("SI","SG","SSI","SSG")
mcols <- c("#B897DD","#BCDD97","#B43FEE","#78EE3F")
mshapes <- c(25,10,17,16)

#run microeco faprotax
data <- phyloseq2meco(ps)
tf <- trans_func$new(data)
tf$cal_spe_func(prok_database = "FAPROTAX")
tf$cal_spe_func_perc(abundance_weighted = T) 
tf$trans_spe_func_perc()
saveRDS(tf, 'data/16Sfaprotax/faprotax_abund.rds')

#process data for alluvial plots
vbles <- unique(tf$res_spe_func_perc_trans$variable)
asv_idx <- which(rowSums(tf$res_spe_func[,vbles])>0)
asv_2_cat <- tf$res_spe_func[,vbles] %>% rownames_to_column("ASV") %>% 
  pivot_longer(!ASV, names_to = 'variable', values_to = 'present') %>% 
  filter(present>0) %>% select(!present)
funcdat <- tf$res_spe_func_perc_trans %>% rename(SampleID=sampname) %>% 
  merge(., asv_2_cat, by='variable')
tax <- tf$tax_table[asv_idx,] %>% rownames_to_column("ASV")
gdat <- otu[,asv_idx] %>% rownames_to_column("SampleID") %>% 
  pivot_longer(!SampleID, names_to = 'ASV', values_to = 'Abundance') %>% 
  group_by(SampleID) %>% mutate(Abundance=Abundance/sum(Abundance)) %>% 
  left_join(., meta, by='SampleID') %>% 
  left_join(., tax, by='ASV') 


#first, the N Surface plot
subasv <- funcdat %>% filter(group=='N-cycle') %>% pull(ASV) %>% unique()
tdat <- gdat %>% filter(ASV %in% subasv) %>% 
  left_join(., funcdat %>% filter(group=="N-cycle"), by=c("SampleID","ASV")) %>% 
  mutate(weight = Abundance * value, tax=paste(Phylum, Genus, sep = '_')) %>% 
  group_by(Microhabitat, group, variable, tax) %>%
  summarize(weight=sum(weight)) %>% 
  filter(weight>1, Microhabitat %in% c("SI","SG"))

ggplot(tdat, aes(y=weight, axis1=Microhabitat, axis2=group, axis3=variable, axis4=tax)) +
  geom_alluvium(aes(fill=Microhabitat)) +
  geom_stratum()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, angle=45) +
  scale_x_discrete(limits=c("Microhabitat", "Cycle","Function", "Tax"))+
  theme_bw()+
  ggtitle("Surface Nitrogen Cycle")
ggsave('figs/16S_Surface_N_Cycle_alluv.png')

#subsurface
tdat <- gdat %>% filter(ASV %in% subasv) %>% 
  left_join(., funcdat %>% filter(group=="N-cycle"), by=c("SampleID","ASV")) %>% 
  mutate(weight = Abundance * value, tax=paste(Phylum, Genus, sep = '_')) %>% 
  group_by(Microhabitat, group, variable, tax) %>%
  summarize(weight=sum(weight)) %>% 
  filter(weight>1, Microhabitat %in% c("SSI","SSG"))

ggplot(tdat, aes(y=weight, axis1=Microhabitat, axis2=group, axis3=variable, axis4=tax)) +
  geom_alluvium(aes(fill=Microhabitat)) +
  geom_stratum()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, angle=45) +
  scale_x_discrete(limits=c("Microhabitat", "Cycle","Function", "Tax"))+
  theme_bw()+
  ggtitle("Subsurface Nitrogen Cycle")
ggsave('figs/16S_Subsurface_N_Cycle_alluv.png')

#first, the C Surface plot
subasv <- funcdat %>% filter(group=='C-cycle') %>% pull(ASV) %>% unique()
tdat <- gdat %>% filter(ASV %in% subasv) %>% 
  left_join(., funcdat %>% filter(group=="C-cycle"), by=c("SampleID","ASV")) %>% 
  mutate(weight = Abundance * value, tax=paste(Phylum, Genus, sep = '_')) %>% 
  group_by(Microhabitat, group, variable, tax) %>%
  summarize(weight=sum(weight)) %>% 
  filter(weight>1, Microhabitat %in% c("SI","SG"))

ggplot(tdat, aes(y=weight, axis1=Microhabitat, axis2=group, axis3=variable, axis4=tax)) +
  geom_alluvium(aes(fill=Microhabitat)) +
  geom_stratum()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, angle=45) +
  scale_x_discrete(limits=c("Microhabitat", "Cycle","Function", "Tax"))+
  theme_bw()+
  ggtitle("Surface Carbon Cycle")
ggsave('figs/16S_Surface_C_Cycle_alluv.png')

#subsurface
tdat <- gdat %>% filter(ASV %in% subasv) %>% 
  left_join(., funcdat %>% filter(group=="C-cycle"), by=c("SampleID","ASV")) %>% 
  mutate(weight = Abundance * value, tax=paste(Phylum, Genus, sep = '_')) %>% 
  group_by(Microhabitat, group, variable, tax) %>%
  summarize(weight=sum(weight)) %>% 
  filter(weight>1e-1, Microhabitat %in% c("SSI","SSG"))

ggplot(tdat, aes(y=weight, axis1=Microhabitat, axis2=group, axis3=variable, axis4=tax)) +
  geom_alluvium(aes(fill=Microhabitat)) +
  geom_stratum()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, angle=45) +
  scale_x_discrete(limits=c("Microhabitat", "Cycle","Function", "Tax"))+
  theme_bw()+
  ggtitle("Subsurface Carbon Cycle")
ggsave('figs/16S_Subsurface_C_Cycle_alluv.png')


#first, the S Surface plot
subasv <- funcdat %>% filter(group=='S-cycle') %>% pull(ASV) %>% unique()
tdat <- gdat %>% filter(ASV %in% subasv) %>% 
  left_join(., funcdat %>% filter(group=="S-cycle"), by=c("SampleID","ASV")) %>% 
  mutate(weight = Abundance * value, tax=paste(Phylum, Genus, sep = '_')) %>% 
  group_by(Microhabitat, group, variable, tax) %>%
  summarize(weight=sum(weight)) %>% 
  filter(Microhabitat %in% c("SI","SG"))

ggplot(tdat, aes(y=weight, axis1=Microhabitat, axis2=group, axis3=variable, axis4=tax)) +
  geom_alluvium(aes(fill=Microhabitat)) +
  geom_stratum()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, angle=45) +
  scale_x_discrete(limits=c("Microhabitat", "Cycle","Function", "Tax"))+
  theme_bw()+
  ggtitle("Surface Sulphur Cycle")
ggsave('figs/16S_Surface_S_Cycle_alluv.png')

#subsurface
tdat <- gdat %>% filter(ASV %in% subasv) %>% 
  left_join(., funcdat %>% filter(group=="S-cycle"), by=c("SampleID","ASV")) %>% 
  mutate(weight = Abundance * value, tax=paste(Phylum, Genus, sep = '_')) %>% 
  group_by(Microhabitat, group, variable, tax) %>%
  summarize(weight=sum(weight)) %>% 
  filter(Microhabitat %in% c("SSI","SSG"))

ggplot(tdat, aes(y=weight, axis1=Microhabitat, axis2=group, axis3=variable, axis4=tax)) +
  geom_alluvium(aes(fill=Microhabitat)) +
  geom_stratum()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, angle=45) +
  scale_x_discrete(limits=c("Microhabitat", "Cycle","Function", "Tax"))+
  theme_bw()+
  ggtitle("Subsurface Sulphur Cycle")
ggsave('figs/16S_Subsurface_S_Cycle_alluv.png')


#surface by treatment N-cycle
subasv <- funcdat %>% filter(group=='N-cycle') %>% pull(ASV) %>% unique()
tdat <- gdat %>% filter(ASV %in% subasv, Surface_Subsurface=='S') %>% 
  left_join(., funcdat %>% filter(group=="N-cycle"), by=c("SampleID","ASV")) %>% 
  mutate(weight = Abundance * value, 
         tax=paste(Phylum, Genus, sep = '_'),
         Treatment=factor(Treatment, levels=treatlevs)) %>% 
  group_by(Treatment, group, variable, tax) %>%
  summarize(weight=sum(weight)) %>% 
  filter(weight>1)
ggplot(tdat, aes(y=weight, axis1=Treatment, axis2=group, axis3=variable, axis4=tax)) +
  geom_alluvium(aes(fill=Treatment)) +
  geom_stratum()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, angle=45) +
  scale_x_discrete(limits=c("Treatment", "Cycle","Function", "Tax"))+
  scale_fill_manual(values=treatcols)+
  theme_bw()+
  ggtitle("Surface Nitrogen Cycle (Treatments)")
ggsave('figs/16S_Surface_N_Cycle_treat_alluv.png')

#subsurface
tdat <- gdat %>% filter(ASV %in% subasv, Surface_Subsurface=='SS') %>% 
  left_join(., funcdat %>% filter(group=="N-cycle"), by=c("SampleID","ASV")) %>% 
  mutate(weight = Abundance * value, 
         tax=paste(Phylum, Genus, sep = '_'),
         Treatment=factor(Treatment, levels=treatlevs)) %>% 
  group_by(Treatment, group, variable, tax) %>%
  summarize(weight=sum(weight)) %>% 
  filter(weight>1)
ggplot(tdat, aes(y=weight, axis1=Treatment, axis2=group, axis3=variable, axis4=tax)) +
  geom_alluvium(aes(fill=Treatment)) +
  geom_stratum()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, angle=45) +
  scale_x_discrete(limits=c("Treatment", "Cycle","Function", "Tax"))+
  scale_fill_manual(values=treatcols)+
  theme_bw()+
  ggtitle("Subsurface Nitrogen Cycle (Treatments)")
ggsave('figs/16S_Subsurface_N_Cycle_treat_alluv.png')

#surface by treatment C-cycle
subasv <- funcdat %>% filter(group=='C-cycle') %>% pull(ASV) %>% unique()
tdat <- gdat %>% filter(ASV %in% subasv, Surface_Subsurface=='S') %>% 
  left_join(., funcdat %>% filter(group=="C-cycle"), by=c("SampleID","ASV")) %>% 
  mutate(weight = Abundance * value, 
         tax=paste(Phylum, Genus, sep = '_'),
         Treatment=factor(Treatment, levels=treatlevs)) %>% 
  group_by(Treatment, group, variable, tax) %>%
  summarize(weight=sum(weight)) %>% 
  filter(weight>.2)
ggplot(tdat, aes(y=weight, axis1=Treatment, axis2=group, axis3=variable, axis4=tax)) +
  geom_alluvium(aes(fill=Treatment)) +
  geom_stratum()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, angle=45) +
  scale_x_discrete(limits=c("Treatment", "Cycle","Function", "Tax"))+
  scale_fill_manual(values=treatcols)+
  theme_bw()+
  ggtitle("Surface Carbon Cycle (Treatments)")
ggsave('figs/16S_Surface_C_Cycle_treat_alluv.png')

#subsurface
tdat <- gdat %>% filter(ASV %in% subasv, Surface_Subsurface=='SS') %>% 
  left_join(., funcdat %>% filter(group=="C-cycle"), by=c("SampleID","ASV")) %>% 
  mutate(weight = Abundance * value, 
         tax=paste(Phylum, Genus, sep = '_'),
         Treatment=factor(Treatment, levels=treatlevs)) %>% 
  group_by(Treatment, group, variable, tax) %>%
  summarize(weight=sum(weight)) %>% 
  filter(weight>1e-2)
ggplot(tdat, aes(y=weight, axis1=Treatment, axis2=group, axis3=variable, axis4=tax)) +
  geom_alluvium(aes(fill=Treatment)) +
  geom_stratum()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, angle=45) +
  scale_x_discrete(limits=c("Treatment", "Cycle","Function", "Tax"))+
  scale_fill_manual(values=treatcols[2:length(treatcols)])+
  theme_bw()+
  ggtitle("Subsurface Carbon Cycle (Treatments)")
ggsave('figs/16S_Subsurface_C_Cycle_treat_alluv.png')


#surface by treatment S-cycle
subasv <- funcdat %>% filter(group=='S-cycle') %>% pull(ASV) %>% unique()
tdat <- gdat %>% filter(ASV %in% subasv, Surface_Subsurface=='S') %>% 
  left_join(., funcdat %>% filter(group=="S-cycle"), by=c("SampleID","ASV")) %>% 
  mutate(weight = Abundance * value, 
         tax=paste(Phylum, Genus, sep = '_'),
         Treatment=factor(Treatment, levels=treatlevs)) %>% 
  group_by(Treatment, group, variable, tax) %>%
  summarize(weight=sum(weight)) %>% 
  filter(weight>1e-5)
ggplot(tdat, aes(y=weight, axis1=Treatment, axis2=group, axis3=variable, axis4=tax)) +
  geom_alluvium(aes(fill=Treatment)) +
  geom_stratum()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, angle=45) +
  scale_x_discrete(limits=c("Treatment", "Cycle","Function", "Tax"))+
  scale_fill_manual(values=c("#71C4E5","#E7BCCC"))+
  theme_bw()+
  ggtitle("Surface Sulphur Cycle (Treatments)")
ggsave('figs/16S_Surface_S_Cycle_treat_alluv.png')

#subsurface
tdat <- gdat %>% filter(ASV %in% subasv, Surface_Subsurface=='SS') %>% 
  left_join(., funcdat %>% filter(group=="S-cycle"), by=c("SampleID","ASV")) %>% 
  mutate(weight = Abundance * value, 
         tax=paste(Phylum, Genus, sep = '_'),
         Treatment=factor(Treatment, levels=treatlevs)) %>% 
  group_by(Treatment, group, variable, tax) %>%
  summarize(weight=sum(weight)) %>% 
  filter(weight>1e-5)
ggplot(tdat, aes(y=weight, axis1=Treatment, axis2=group, axis3=variable, axis4=tax)) +
  geom_alluvium(aes(fill=Treatment)) +
  geom_stratum()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, angle=45) +
  scale_x_discrete(limits=c("Treatment", "Cycle","Function", "Tax"))+
  scale_fill_manual(values=treatcols)+
  theme_bw()+
  ggtitle("Subsurface Sulphur Cycle (Treatments)")
ggsave('figs/16S_Subsurface_S_Cycle_treat_alluv.png')


#surface by treatment Energy source
subasv <- funcdat %>% filter(group=='Energy source') %>% pull(ASV) %>% unique()
tdat <- gdat %>% filter(ASV %in% subasv,
                        SampleDate=="Nov_2023") %>% 
  left_join(., funcdat %>% filter(group=="Energy source"), by=c("SampleID","ASV")) %>% 
  mutate(weight = Abundance * value, 
         tax=paste(Phylum, Genus, sep = '_'),
         Treatment=factor(Treatment, levels=treatlevs)) %>% 
  group_by(Treatment, group, variable, tax) %>%
  summarize(weight=sum(weight)) %>% 
  filter(weight>10)
ggplot(tdat %>% filter(variable %in% c("aerobic_chemoheterotrophy", "anaerobic_chemoheterotrophy")), 
       aes(y=weight, axis1=Treatment, axis2=group, axis3=variable, axis4=tax)) +
  geom_alluvium(aes(fill=Treatment)) +
  geom_stratum()+
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits=c("Treatment", "Cycle","Function", "Tax"))+
  scale_fill_manual(values=treatcols[2:length(treatcols)])+
  theme_bw()+
  ggtitle("Energy source (Treatment)")
ggsave('figs/16S_Surface_N_Cycle_treat_alluv.png')




## FAPROTAX Function-level NMDS
tdat <- gdat %>% 
  filter(ASV %in% names(asv_idx)) %>% 
  left_join(., funcdat , by=c("SampleID","ASV")) %>% 
  group_by(SampleID, variable) %>% 
  summarize(Abundance=sum(Abundance)) %>% 
  left_join(., meta, by="SampleID")

#prep distance mat
tmat <- tdat %>% pivot_wider(id_cols = SampleID, names_from = variable, values_from = Abundance) %>% data.frame()
rownames(tmat) <- tmat[,1]
tmat <- tmat[2:ncol(tmat)]
X <- vegdist(tmat, method='bray')

#run NMDS
out <- metaMDS(X,
               k=3,
               wascores = T,
               weakties=T,
               try=50,
               trymax=100,
               parallel=8,
               maxit=300)
print(paste0("Stress of NMDS model: ", out$stress))

gdat <- data.frame(out$points)
gdat$goodness <- goodness(out)
gdat$Sample <- rownames(gdat)
gdat <- left_join(gdat, meta %>% rename(Sample=SampleID), by='Sample') %>% 
  mutate(SampleDate=factor(SampleDate, levels=timelevs))

ggplot(gdat)+
  geom_jitter(aes(x=SampleDate, y=goodness, color=Treatment, shape=Microhabitat))+
  theme_bw()+
  scale_color_manual(values = treatcols)+
  scale_shape_manual(values = mshapes)+
  xlab("Date")+
  ylab("Goodness of Fit")


stressplot(out)

ggplot(gdat)+
  geom_point(aes(x=MDS1, y=MDS2, color=Treatment, shape=Microhabitat),
             size=3, alpha=0.85)+
  theme_bw()+
  theme(text = element_text(size=15))+
  scale_color_manual(values = treatcols)+
  scale_shape_manual(values = mshapes)+
  facet_wrap(~SampleDate)+
  ggtitle('MDS1 v MDS2 (Function-level)')
