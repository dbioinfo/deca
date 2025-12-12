# Desert Ecology Community Assembly

src/deca_16Sdada2.Rmd
- This script performs qc, error correction, abundance quantification and taxonomic classification. The resulting data product is a phyloseq object containing an otu table, taxonomic classifications, sample metadata and reference sequences.  

src/decontam_16S.Rmd
- We then use the decontam package and our negative control to remove any identified ASVs. We also filter out any samples with <5000 Abundance as well as ASVs unique to a sample.

src/deca_16Smirlyn.Rmd
- We then have to rarefy the samples due to varying library sizes before calculating alpha diversities and differential abundance.  

src/generate_phytree.sh
- Then we create a phylogenetic tree using RAxML.

