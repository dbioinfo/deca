library(tidyverse)
library(colorspace)
library(phyloseq)
library(ggh4x)

setwd('~/WorkForaging/Academia/Nicole/deca/')

ps <- readRDS('data/16Sdada2/phyloseq.filtered.tree.rds')

otu <- t(data.frame(otu_table(ps)))
meta <- data.frame(sample_data(ps))
tax <- data.frame(tax_table(ps)) %>% rownames_to_column("ASV")
timelevs <- c("Jul_2022","Nov_2022","Feb_2023", "Aug_2023","Nov_2023", "Feb_2024","Aug_2024")
meta <- meta %>% mutate(SampleDate = factor(SampleDate, levels=timelevs), 
                        sample=SampleID) %>% arrange(SampleDate) %>% #specific colnames needed for this genius framework
                mutate(collection_date = case_when(SampleDate=="Jul_2022"~as.Date("7/1/22"),
                                                   SampleDate=="Nov_2022"~as.Date("11/1/22"),
                                                   SampleDate=="Feb_2023"~as.Date("2/1/23"),
                                                   SampleDate=="Aug_2023"~as.Date("8/1/23"),
                                                   SampleDate=="Nov_2023"~as.Date("11/1/23"),
                                                   SampleDate=="Feb_2024"~as.Date("2/1/24"),
                                                   SampleDate=="Aug_2024"~as.Date("8/1/24")))


#let's summarize at the 'phylum' level
iotu <- data.frame(otu) %>% rownames_to_column("ASV") %>% relocate(ASV) #no specific colname here, just first col
gotu <- left_join(iotu, tax, by='ASV') %>% 
  mutate(taxon=paste(Kingdom, Phylum, sep="_")) %>% 
  select(taxon, colnames(iotu)) %>%
  select(-ASV) %>%  
  pivot_longer(-taxon, names_to="sample", values_to = "abundance") %>% 
  group_by(sample, taxon) %>% 
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from='sample', values_from='abundance')
  
itax <- tax %>% mutate(taxon=paste(Kingdom, Phylum, sep="_")) %>% select(ASV, taxon)

gdat <- left_join(gotu %>% pivot_longer(-taxon,names_to = "sample", values_to = "abundance"), meta, by='sample')
gdat <- gdat %>% 
  group_by(SampleDate, taxon, Microhabitat) %>% 
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>% 
  group_by(taxon) %>%
  mutate(    z_score = (abundance - mean(abundance, na.rm = TRUE)) / sd(abundance, na.rm = TRUE) ) %>%
  ungroup() %>% 
  mutate(abs_z = abs(z_score), sign_z = case_when(sign(z_score)>0~"Positive", .default = "Negative" ) )

#complete the missing entries
gdat <- gdat %>% 
  group_by(taxon, Microhabitat) %>%
  complete(SampleDate, sign_z = c("Positive", "Negative"), fill = list(z_score = 0)) 

#rewriting the zscores
iotu <- data.frame(otu) %>% rownames_to_column("ASV") %>% relocate(ASV) #no specific colname here, just first col
gotu <- left_join(iotu, tax, by='ASV') %>% 
  mutate(taxon=paste(Kingdom, Phylum, sep="_")) %>% 
  select(taxon, colnames(iotu)) %>%
  select(-ASV) %>%  
  pivot_longer(-taxon, names_to="sample", values_to = "abundance") %>% 
  group_by(sample, taxon) %>% 
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from='sample', values_from='abundance')
gdat <- left_join( gotu %>% pivot_longer(-taxon, names_to = "sample", values_to = "abundance"),
                    meta,
                    by = "sample") %>%
                    group_by(SampleDate, taxon, Microhabitat) %>%
                    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

all_dates <- sort(unique(meta$SampleDate))
all_microhabs <- sort(unique(meta$Microhabitat))
all_taxa <- sort(unique(gdat$taxon))

gdat <- gdat %>%
  complete(
    taxon = all_taxa,
    Microhabitat = all_microhabs,
    SampleDate = all_dates,
    fill = list(abundance = 0)
  )
gdat <- gdat %>%
  group_by(taxon) %>%
  mutate(
    z_score = {
      mu <- mean(abundance, na.rm = TRUE)
      sig <- sd(abundance, na.rm = TRUE)
      
      if (is.na(sig) || sig == 0) {
        rep(0, dplyr::n())
      } else {
        (abundance - mu) / sig
      }
    }
  ) %>%
  ungroup() %>%
  tidyr::crossing(sign_z = c("Positive", "Negative")) %>%
  mutate(
    z_plot = dplyr::case_when(
      sign_z == "Positive" & z_score > 0 ~ z_score,
      sign_z == "Negative" & z_score < 0 ~ z_score,
      TRUE ~ 0
    ),
    absent = abundance == 0
  )

#poor man's biomehorizon
samp_phy <-gdat %>% group_by(taxon) %>% summarize(z_score=sum(z_score)) %>% filter(z_score!=0) %>% pull(taxon)
samp_phy <- sample(samp_phy, 10) #replace with list of chosen kingdom-phylum pairs c("Bacteria_Cyanobacteriota","Archaea_Thermoplasmatota")
tmp <- gdat %>% filter(taxon %in% samp_phy, !is.na(sign_z))
gg <- ggplot(tmp, 
       aes(x = SampleDate, y = z_plot, group = sign_z)) +
  geom_ribbon(aes(ymin = 0, ymax = z_plot, fill = sign_z), alpha = 0.9) +
  theme_bw() +
  scale_x_discrete(expand=c(0,0))+
  labs(title = "10 Random Phyla Timeseries",
    x = "Date",
    y = "Z-Score (Std Dev from Mean)" ) +
  scale_fill_manual(values = c("#81322C", "#2C812C"))+
  theme(axis.text.x = element_text(size=10,angle=60, vjust=0.3),
        strip.text = element_text(size=15),
        axis.text.y=element_text(size=10),
        title = element_text(size=20))+
  facet_grid2(taxon~Microhabitat, scales = "free_y")
gg
ggsave("figs/microhabitat_biomehorizon.png",gg,height=7000, width=5000, units = 'px')

