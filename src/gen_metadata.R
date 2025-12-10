library(tidyverse)
setwd('~/WorkForaging/Academia/Nicole/NOPe/')

#import
manif <- read_tsv('raw_data/manifest.tsv', col_names=c('Bytes','fpath') )

#format
manif <- manif %>%   mutate(fname = sub("^.*/", "", fpath)) %>%   
  separate(fname, into = c("uknwn1", "sid", "uknwn2", "drop4", "read"), 
           sep = "_", remove = FALSE, fill = "right", extra = "drop") %>%
  select(-drop4)  %>% 
  mutate(
    Microenvironment = str_extract(sid, "[A-Za-z]+"),
    uknwn3 = str_extract(sid, "[0-9]+")
  ) %>% 
  mutate(sample_name = paste(uknwn1,sid,sep='_'),
         Microenvironment = case_when(Microenvironment=='S'~'C', .default=Microenvironment))

#check some stats for sanity
loners <- manif %>% group_by(sample_name) %>% summarise(nreads = n()) %>% filter(nreads<2) %>% pull(sample_name)
manif <-manif  %>% mutate(paired = case_when(sample_name %in% loners ~ F, .default=T)) #add a column indicating unpaired

ampd <- manif %>% group_by(sample_name) %>% summarise(nreads = n()) %>% filter(nreads>2) %>% pull(sample_name)
manif <- manif %>% mutate(amplified = case_when(sample_name %in% ampd ~T, .default=F) )

manif %>% group_by(Microenvironment) %>% summarize(n=n()) #how many samples of each microenvironment do we have

manif %>% filter(Bytes < 1000) %>% nrow() #how many files are probably too small to contain any usable information?

#out
write_csv(manif, 'data/metadata.csv')
