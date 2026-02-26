set -o pipefail

export TMPDIR=/scratch/han_lab/dbarth/nicole/deca/tmp
export TMP=/scratch/han_lab/dbarth/nicole/deca/tmp
export TEMP=/scratch/han_lab/dbarth/nicole/deca/tmp
export MPLCONFIGDIR=/scratch/han_lab/dbarth/nicole/deca/tmp
export APPTAINER_TMPDIR=/scratch/han_lab/dbarth/nicole/deca/tmp
export APPTAINER_CACHEDIR=/scratch/han_lab/dbarth/nicole/deca/tmp


# Import FASTA as QIIME2 artifact
#qiime tools import \
#  --type 'FeatureData[Sequence]' \
#  --input-path data/16Sqiime/asv-seqs.fasta \
#  --output-path data/16Sqiime/rep-seqs.qza

# Run SEPP fragment insertion using a reference database artifact
qiime fragment-insertion sepp \
  --i-representative-sequences data/16Sqiime/rep-seqs.qza \
  --i-reference-database raw_data/silva_ref/sepp-refs-silva-128.qza \
  --p-threads 80 \
  --o-tree data/16Sqiime/insertion-tree.qza \
  --o-placements data/16Sqiime/insertion-placements.qza

# Export the resulting tree to Newick
qiime tools export \
  --input-path data/16Sqiime/insertion-tree.qza \
  --output-path data/16Sqiime/asv-seqs.newick


