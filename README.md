# Desert Ecology Community Assembly

src/deca_16Sdada2.Rmd
- This script performs qc, error correction, abundance quantification and taxonomic classification. The resulting data product is a phyloseq object containing an otu table, taxonomic classifications, sample metadata and reference sequences.  

src/decontam_16S.Rmd
- We then use the decontam package and our negative control to remove any identified ASVs. We also filter out any samples with <2500 Abundance as well as ASVs unique to a sample.

src/deca_16Smirlyn.Rmd
- We then have to rarefy the samples due to varying library sizes before calculating alpha diversities and differential abundance.  

src/deca_16Sordination.Rmd
- We then compute BC distance and perform an NMDS ordination analysis. 

src/deca_16Spermanova.Rmd
- From this BC distance matrix, we can order the variables in our metadata by how well they separate the data using adonis2. 

src/deca_16Sancombc.Rmd
- Then we generate several runs of ancombc2 to test for differential expression using these ordered variables. 

src/deca_16S_faprotax_eda.Rmd
- Then we run faprotax and generate an interactive UI to aid in data exploration. Here, each ASV is grouped by inclusion in functional cycles and specific processes, as well as taxonomic rank. Alluvial plots are generated to visualize relative abundance in functional space. 

src/deca_16S_ancombc_eda.Rmd
- We also generate an interactive UI for the ancombc2 results using bubble plots that allow the user to explore the data visually. 

