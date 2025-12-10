library(tidyverse)
library(dada2)
library(phyloseq)
library(Biostrings)
library(colorspace)
library(vegan)

#import metadata
meta <- read_csv('raw_data/Raw_Seq_Files_and_Metadata_V2.csv') 
meta <- meta %>% mutate(SampleID=str_remove(R1, "\\_L001_R1_001.fastq.gz$")) %>% as.data.frame()
rownames(meta) <- meta$SampleID

#load forward and reverse reads
fnFs <- paste0('raw_data/250204_Pietrakiak_demultiplxed/',meta$R1) 
fnRs <- paste0('raw_data/250204_Pietrakiak_demultiplxed/',meta$R2)


#plot qc
#for (i in 1:5){
#    png(paste0('figs/dada2_fwd_qc',i,'.png'))
#    print(plotQualityProfile(fnFs[(i*100):(i*100 +3)]))
#    dev.off()
#    png(paste0('figs/dada2_rev_qc',i,'.png'))
#    print(plotQualityProfile(fnRs[(i*100):(i*100+3)]))
#    dev.off()
#}

#set params according to qc
trunc_f <- 250
trunc_r <- 250
maxEE_f <- 2
maxEE_r <- 2


#filter by params
#filtFs <- file.path("raw_data/filt_fastq/", meta$R1)
#filtRs <- file.path("raw_data/filt_fastq/", meta$R2)
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen = c(trunc_f, trunc_r),
#                     maxN=0, maxEE=c(maxEE_f,maxEE_r), rm.phix=TRUE,
#                     compress=TRUE, multithread=TRUE)
#write.csv(out, 'data/16Sdada2/dada2_reads_filt.csv', row.names = T)
#3 samples removed A24_53SSI_S920_L001_R2_001 A24_52SSI_S916_L001_R2_001 A24_53SSG_S918_L001_R2_001

out <- read.csv('data/16Sdada2/dada2_reads_filt.csv')
out <- out %>% mutate(SampleID=str_remove(X,"\\_L001_R1_001.fastq.gz$"))

#prep matched files 
filt_fs <- c()
filt_rs <- c()
for (sname in out$SampleID){
    if ( (file.exists(paste0('raw_data/filt_fastq/',sname,'_L001_R1_001.fastq.gz'))) && (file.exists(paste0('raw_data/filt_fastq/',sname,'_L001_R2_001.fastq.gz')))) {
        filt_fs <- c(filt_fs, paste0('raw_data/filt_fastq/',sname,'_L001_R1_001.fastq.gz'))
        filt_rs <- c(filt_rs, paste0('raw_data/filt_fastq/',sname,'_L001_R2_001.fastq.gz'))
    }
}
#learn error rates (~10 min run time w/32 cores)
#err_f <- learnErrors(filt_fs, multithread = T)
#err_r <- learnErrors(filt_rs, multithread = T)

#run dada2
#dadaFs <- dada(filt_fs, err_f, multithread = T)
#dadaRs <- dada(filt_rs, err_r, multithread = T)

#intout <- list(err_f=err_f, err_r=err_r, dadaFs=dadaFs, dadaRs=dadaRs)
#saveRDS(intout, 'data/16Sdada2/dada2_errs.rds')

#intout <- readRDS('data/16Sdada2/dada2_errs.rds')
#list2env(intout, .GlobalEnv)

#mergers <- mergePairs(dadaFs, filt_fs, dadaRs, filt_rs, verbose=TRUE)
#save
#saveRDS(mergers,'data/16Sdada2/merged_reads_dada2.rds')
mergers <- readRDS('data/16Sdada2/merged_reads_dada2.rds')

seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
1 - sum(seqtab.nochim)/sum(seqtab) #what percent of reads were chimera/bimera

#assign taxa from SILVA
tax <- assignTaxonomy(seqtab.nochim, "raw_data/silva_ref/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)
tax <- addSpecies(tax, "raw_data/silva_ref/silva_v138.2_assignSpecies.fa.gz")

#create phyloseq obj
rownames(seqtab.nochim) <- gsub('_L001_R1_001.fastq.gz','',rownames(seqtab.nochim))
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows = F), sample_data(meta), tax_table(tax))

#fix ASV names, exact seqs are annoying
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

#save phyloseq
otu <- data.frame(otu_table(ps))
otu$Sample <- rownames(otu)
otu <- otu %>% relocate(Sample)
write_csv(otu,'data/16Sdada2/dada2_otu.csv')
saveRDS(ps, 'data/16Sdada2/phyloseq_dada2.rds')





