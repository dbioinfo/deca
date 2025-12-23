library(tidyverse)
library(dbplyr)
library(DBI)
library(RPostgres)

#set working dir
setwd('~/WorkForaging/Academia/Nicole/deca')

#load data
ps <- readRDS('data/16Sdada2/phyloseq_dada2.rds')
otu <- read_csv('data/16Sdada2/dada2_otu.csv')

# Connect to the deca_db
con <- dbConnect(
  RPostgres::Postgres(),
  dbname = "deca_db",
  host   = "127.0.0.1",
  port   = 1738,
  user   = 'postgres' )

#store otu in viable format (PSQL needs <1.6k columns)
otu <- otu %>% pivot_longer(-Sample, names_to = 'ASV', values_to = 'Abundance')
dbWriteTable(con, 'otu', otu, overwrite=T)

#metadata and tax are fine as they are
meta <- data.frame(sample_data(ps)) %>% rename(Sample=SampleID)
dbWriteTable(con, 'meta', meta, overwrite=T)

tax <- data.frame(tax_table(ps))
tax$ASV <- rownames(tax)
dbWriteTable(con, 'tax', tax, overwrite=T)

#add index for much faster searching
dbExecute(con, 'CREATE INDEX IF NOT EXISTS idx_otu ON otu ("Sample", "ASV")')

#quick query to check speed
tmp <- tbl(con, 'otu')
ssg <- tbl(con, 'meta') %>% filter(Microhabitat=='SSG') %>% pull(Sample)
taxs <- tbl(con, 'tax') %>% filter(Phylum=='Thermoproteota') %>% pull(ASV)
#how many thermoproteota species are found in all ssg samples? (did that run quickly? if so, proceed, if not, db needs optimization)
tmp %>% filter(Sample %in% ssg, ASV %in% taxs, Abundance>0) %>% tally()

#now write all the mirlyn reps to the database
mirl <- readRDS('data/mirlyn/mirl.rds')
dbExecute(con, 'DROP TABLE IF EXISTS mirl;')
for (i in 1:length(mirl)){
  tmp <- data.frame(otu_table(mirl[[i]]))
  tmp$rep <- i
  dbWriteTable(con, 'mirl', tmp, append=T)
  print(i)
}


