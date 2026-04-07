library(tidyverse)
library(phyloseq)
library(treeio)
library(ggtreeExtra)
library(tidytree)
library(ggnewscale)
library(ggiraph)
library(iggtree)
library(ggtree)

setwd("~/WorkForaging/Academia/Nicole/deca/")

## Import all the data
#### Phyloseq
ps <- readRDS('data/16Sdada2/phyloseq.filtered.tree.rds') 

#### ANCOMBC
ancom_out <- data.frame() 
itaxlevs <- c("phylum","genus")
variables <- c("Microhabitat","Treatment","Surface_Subsurface","Grass_Interspace","SampleDate")
for (ivar in variables){
  
  #handle each var 
  loadpath <- 'data/16Sancombc/ancombc_results'
  if (ivar=='Microhabitat') { iref <- "SI"; ipath <- 'mh' } else 
    if (ivar=="Treatment") { iref <- "Control"; ipath <- 'treat'} else 
      if (ivar=="Surface_Subsurface") {iref <- "S"; ipath <- 'ss'} else 
        if (ivar=="Grass_Interspace") {iref <- "I"; ipath <- 'gi'} else
          if (ivar=="SampleDate") {iref <- "Jul_2022"; ipath <- 'date'}
  
  
  #iter through taxlevs
  for (itax in itaxlevs){
    
    out <- readRDS(paste0( paste(loadpath, ipath, itax, sep = '_'), '.rds')  )
    
    
    if (ivar %in% c("Surface_Subsurface","Grass_Interspace")) { #handle vars with 2 levels
      res <- out$res %>% select(starts_with(c('tax',"lfc","diff", "p","W"))) %>%  
        select(contains(c("taxon",ivar)))
    } else {
      res <- out$res_pair %>% select(starts_with(c('tax',"lfc","diff", "p","W"))) 
    }
    
    res_long <- res %>% #rename implicit cols to make reference column explicit
      rename_with( ~{ 
        x <- .x
        term <- ivar #category
        ref <- iref #reference condition
        
        #regex string clean (special handle for dates) 
        if (term == 'SampleDate') { lvl_pat <- "(?![A-Za-z]{3}_\\d{4}_[A-Za-z]{3}_\\d{4}$)([^_]+(?:_[^_]+)?)$" } else {  lvl_pat <- "([A-Za-z0-9\\+]+)"}
        x <- str_replace(x, 
                         pattern = paste0("^(lfc|diff|p|W)_", term, lvl_pat, "$"), 
                         replacement = paste0("\\1_\\2_", ref))
        x <- str_replace_all(x, paste0(term), "")
        x <- str_replace(x, "^p_", "p_val_")
        
        
        x
      } , .cols = -taxon) %>%
      #pivot longer for graphing
      pivot_longer(cols = -taxon,
                   names_to = c(".value", "contrast"),
                   names_pattern = "^(lfc|diff|p_val|W)_(.+)$") %>%
      #clean up, add variable column
      mutate(contrast = str_replace_all(contrast, "^_", ""),  
             contrast = str_replace_all(contrast, "__+", "_"), 
             variable=ivar,
             tax=itax) 
    
    ancom_out <- rbind(ancom_out, res_long)
  }
} #when selecting data from this output, you must choose tax (level), variable, AND contrast

#### Faprotax
fap <-readRDS('data/16Sfaprotax/faprotax_abund.rds')
fap <- fap$res_spe_func


#### ANCOMBC: pick a variable, contrast, and tax level
ivar <- "Surface_Subsurface"
icontrast <- "SS_S"
itaxlev <- "phylum"

sps <- tax_glom(ps, taxrank=str_to_title(itaxlev), NArm=T)
tree <- phy_tree(sps)
tax <- tax_table(sps) %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("label")
gdat <- left_join(as_tibble(tree), tax, by = "label")

#add ancombc result
tmp <- ancom_out %>% filter(variable==ivar,tax==itaxlev,contrast==icontrast)
if (itaxlev=="phylum"){
  gdat <- gdat %>% mutate(taxon=paste(Kingdom, Phylum, sep="_"))
  } else if (itaxlev=="genus") {
  gdat <- gdat %>% mutate(taxon=paste(Kingdom, Phylum, Class, Order, Family, Genus, sep="_"))
  }
gdat <- left_join(gdat, tmp, by='taxon')

#add Faprotax: pick tax level and variables to viz
ivars <- c("phototrophy","nitrogen_respiration","iron_respiration")
ivars <- unique(c(ivars, sample( colnames(fap[which(colSums(fap)>0)]) , 10))) #choose 10 random non-0
tmp <- fap[,ivars] %>% rownames_to_column("label")
itax <- data.frame(tax_table(ps)) %>% rownames_to_column("label")
tmp <- left_join(tmp, itax, by='label')

if (itaxlev=="phylum"){
  tmp <- tmp %>% mutate(taxon=paste(Kingdom, Phylum, sep="_"))
} else if (itaxlev=="genus") {
  tmp <- tmp %>% mutate(taxon=paste(Kingdom, Phylum, Class, Order, Family, Genus, sep="_"))
}
tmp <- tmp %>% group_by(taxon) %>% summarize(across(all_of(ivars), sum, na.rm=T))

gdat <- left_join(gdat, tmp, by='taxon')

#add rel_abundance: pick metadata conditions
meta <- data.frame(sample_data(ps))
isamps <- meta %>% filter(SampleDate=="Jul_2022", Surface_Subsurface=="S") %>% pull(SampleID)

##test env
otu <- data.frame(otu_table(sps))[isamps, ] %>% 
  rownames_to_column('SampleID') %>% 
  left_join(., meta %>% select(SampleID, Surface_Subsurface), by="SampleID") %>% 
  group_by(Surface_Subsurface) %>% 
  summarise(across(where(is.numeric), sum), .groups = "drop") %>% 
  mutate(across(where(is.numeric), ~ .x / sum(.x))) %>% 
  pivot_longer(!Surface_Subsurface, names_to="label", values_to="rel_abundance") %>% 
  left_join(., gdat %>% select(node, label), by='label')


otu <- colSums(data.frame(otu_table(sps))[isamps, ])
otu <- otu / sum(otu)
otu <- data.frame(rel_abundance=otu) %>% rownames_to_column("label")
gdat <- left_join(gdat, otu, by="label")

#format the data for viz
is.waive <- ggplot2::is_waiver
tdat <- gdat %>% as.treedata()


## Main Plot Generation

#first generate the main tree plot
gg<- ggtree(tdat, layout='circular')

#next add on the ancombc results
gg <- gg +
  geom_fruit(mapping=aes(fill=diff), geom = geom_tile, pwidth=0.01, width=0.05)+
  scale_fill_manual(values=c("#E2908D","#89AAE1","#FFFFFF"),name="Significant(SS,S)",guide = guide_legend(order = 1))+
  new_scale_fill()+
  geom_fruit(mapping=aes(fill=lfc), geom = geom_tile, pwidth=0.075, width=0.2)+
  scale_fill_gradient2(low='#A22E2A',mid='#FFFFFF',high='#2A56A2',midpoint=0, name="LFC(SS,S)",guide = guide_colorbar(order = 2))
  
#add on specific faprotax  
gg <- gg +
  new_scale_fill()+
  geom_fruit(mapping=aes(fill=phototrophy), color='black',linewidth=0.1, geom = geom_tile, pwidth = 0.01, width=.1, offset = 25e-3)+ #initial offset v import
  scale_fill_gradient2(low='#FFFFFF', high="#53AC53",name="phototrophy (N ASV)",guide = guide_colorbar(order = 3))+
  new_scale_fill()+
  geom_fruit(mapping=aes(fill=nitrogen_respiration), color='black',linewidth=0.1, geom = geom_tile, pwidth = 0.01, width=.1, offset = 25e-3)+
  scale_fill_gradient2(low='#FFFFFF', high="#E2B78D",name="nitrogen respiration (N ASV)",guide = guide_colorbar(order = 4))+
  new_scale_fill()+
  geom_fruit(mapping=aes(fill=iron_respiration), color='black',linewidth=0.1, geom = geom_tile, pwidth = 0.01, width=.1, offset = 25e-3)+
  scale_fill_gradient2(low='#FFFFFF', high="#D98DE2",name='iron respiration (N ASV)',guide = guide_colorbar(order = 5))

#add in random faprotax
iorder <- 1
for (i in setdiff(ivars, c("phototrophy","nitrogen_respiration","iron_respiration"))){ 
  gg <- gg +
    new_scale_fill()+
    geom_fruit(mapping=aes(fill= !!sym(i),
                           tooltip = paste0(str_to_title(itaxlev),
                                            ": ", taxon, '\n ',
                                            "N ", i," ASV: ", !!sym(i))), 
               color='black',
               geom = geom_tile, linewidth=0.1, pwidth = 0.02, width=0.1, offset=25e-3, show.legend=F)+
    scale_fill_gradient2(low='#FFFFFF', high="#1F7A57",name=paste0(i," (N ASV)"),guide = guide_colorbar(order = 5+iorder))
  iorder <- iorder + 1
}

##experiment
gg <- gg+
  new_scale_fill()+
  geom_fruit(data=otu, mapping=aes(fill=Surface_Subsurface, y=label, x=rel_abundance),
             geom = geom_bar,orientation = "y", stat='identity', 
             position = position_dodge(preserve = "total"),
             pwidth=0.1, offset=.5, show.legend=F)

#add in rel abundance
gg <- gg +
  geom_fruit(mapping=aes(size=rel_abundance), geom=geom_point, offset = 0.17, pwidth=.05)+
  scale_size(name='Jul_2022 Surface Relative Abundance')

#add the labels
gg <- gg+
  geom_tiplab(mapping=aes(label=taxon),size = 2, pwidth=.5, offset = 4.5, align = TRUE,linetype = 'blank') 

#fix theme elements
gg <- gg+
  theme(legend.margin = margin(t = 0, r = 0, b = 0, l = 2),
        plot.margin = margin(1, 1, 1, 1, "cm"),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.key = element_rect(fill = "transparent", color = NA),
        text=element_text(size=5))

#gg

#ggsave("figs/16S_example_phylo.png",gg,width = 7000,height=5500,units='px')

girafe(ggobj = gg, 
       options = list(
         width_svg = 10,  # Increase the internal "canvas" width (in inches)
         height_svg = 8, # Increase height to match your 7k:5.5k ratio (~1.27)
         opts_zoom(max = 5), # Allows the user to zoom in up to 5x
         opts_toolbar(position = "topright") # Adds the zoom/reset buttons
       ))
