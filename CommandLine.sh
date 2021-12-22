# Title: Evolution of pesticide tolerance and associated changes in the microbiome in Daphnia magna
# Authors: Lizanne Janssens, Marlies Van de Maele, Vienna Delnat, Charlotte Theys, Shinjini Mukherjee, Luc De Meester, Robby Stoks
# Command line used on the VSC (Vlaams Supercomputer Centrum) using MobaXterm
# Based on the Physalia course ‘16 S/ITS Metabarcoding Of Microbial Communities’


# Add fastq.gz files to DATA server

###### Install Miniconda on the VSC to create QIIME2 environment ######

# Download the most recent version of Miniconda (Linux) from link and install on DATA server
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

# Download QIIME2 using Miniconda 
wget https://data.qiime2.org/distro/core/qiime2-2020.6-py36-linux-conda.yml
conda env create -n qiime2-2020.6 --file qiime2-2020.6-py36-linux-conda.yml

# Set suitable UTF-8 supporting locales
export LC_ALL=en_GB.utf8
export LANG=en_GB.utf8

# Activate conda environment qiime2-2020.6
conda activate qiime2-2020.6

# Work on the SCRATCH server from now on using batch files (.pbs)


###### Importing and visualizing reads in QIIME2 ######

# https://docs.qiime2.org/2020.6/tutorials/importing/ (July, 2020)
# Demultiplexed data, “Fastq manifest” formats: PairedEndFastqManifestPhred33V2
# Add the manifest.txt to your working directory

source activate qiime2-2020.6

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path /manifest.txt \
  --output-path /demux-paired-end-20200723.qza

qiime demux summarize \
  --i-data /demux-paired-end-20200723.qza \
  --o-visualization /demux-paired-end-20200723.qzv


###### Denoising using DADA2 in QIIME2 ######

# https://docs.qiime2.org/2020.6/plugins/available/dada2/denoise-paired/ (July, 2020)

source activate qiime2-2020.6

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs /demux-paired-end-20200723.qza \
  --p-n-threads 12 \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 240 \
  --p-trim-left-f 14 \
  --p-trim-left-r 14 \
  --o-table /table_16S_trim14trunc240.qza \
  --o-representative-sequences /rep-seqs_16S_trim14trunc240.qza \
  --o-denoising-stats /denoising-stats_16S_trim14trunc240.qza

qiime metadata tabulate \
  --m-input-file /denoising-stats_16S_trim14trunc240.qza \
  --o-visualization /denoising-stats_16S_trim14trunc240.qzv


###### Summarizing Feature Table and Feature Data ######

# Add the sample-metadata.txt to your working directory
# Details on QIIME 2 metadata requirements: https://docs.qiime2.org/2020.6/tutorials/metadata/
# Validate sample/feature metadata files with Keemei: https://keemei.qiime2.org

source activate qiime2-2020.6

qiime feature-table summarize \
  --i-table /table_16S_trim14trunc240.qza \
  --o-visualization /table_16S_trim14trunc240.qzv \
  --m-sample-metadata-file /sample-metadata.txt

qiime feature-table tabulate-seqs \
  --i-data /rep-seqs_16S_trim14trunc240.qza \
  --o-visualization /rep-seqs_16S_trim14trunc240.qzv


###### Taxonomy assignment ######

# Download full OTU sequences and taxonomy of the Silva 138 SSU Ref NR 99 database 
wget \
  -O "silva-138-99-sequences.qza" \
"https://data.qiime2.org/2020.6/common/silva-138-99-seqs.qza"
wget \
  -O "silva-138-99-taxonomy.qza" \
"https://data.qiime2.org/2020.6/common/silva-138-99-tax.qza"

# Create a naive Bayesian classifier containing the extracted sequences of the V3-V4 region with a forward (CCTACGGGNGGCWGCAG) and a reverse (GACTACHVGGGTATCTAATCC) primer of 16S-IllumTS and a truncation length of 460 bp
# Remark: feature-classifier fit-classifier-naive-bayes needs more than the default 5 gb available memory (Genius server with nodes=1:ppn=1 and pmem=50gb)

source activate qiime2-2020.6

qiime feature-classifier extract-reads \
--i-sequences /silva-138-99-sequences.qza \
--p-f-primer CCTACGGGNGGCWGCAG \
--p-r-primer GACTACHVGGGTATCTAATCC \
--p-trunc-len 460 \
--o-reads /silva-138-99-v3-v4-sequences.qza

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads /silva-138-99-v3-v4-sequences.qza \
--i-reference-taxonomy /silva-138-99-taxonomy.qza \
--o-classifier /silva-138-99-v3-v4-nb-classifier.qza

# Taxonomy assignment
# Remark: this needs more than default available memory (Genius server with nodes=1:ppn=4 and pmem=15gb)

source activate qiime2-2020.6

qiime feature-classifier classify-sklearn \
  --i-classifier /silva-138-99-v3-v4-nb-classifier.qza \
  --i-reads /rep-seqs_16S_trim14trunc240.qza \
  --p-n-jobs 4 \
  --p-reads-per-batch 5000 \
  --o-classification /taxonomy_16S_trim14trunc240.qza

qiime metadata tabulate \
  --m-input-file /taxonomy_16S_trim14trunc240.qza \
  --o-visualization /taxonomy_16S_trim14trunc240.qzv

qiime taxa barplot \
  --i-table /table_16S_trim14trunc240.qza \
  --i-taxonomy /taxonomy_16S_trim14trunc240.qza \
  --m-metadata-file /sample-metadata.txt \
  --o-visualization /taxa-bar-plots_16S_trim14trunc240.qzv


###### Generating a phylogenetic tree ######

source activate qiime2-2020.6

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /rep-seqs_16S_trim14trunc240.qza \
  --o-alignment /aligned-rep-seqs_16S_trim14trunc240.qza \
  --o-masked-alignment /masked-aligned-rep-seqs_16S_trim14trunc240.qza \
  --o-tree /unrooted-tree_16S_trim14trunc240.qza \
  --o-rooted-tree /rooted-tree_16S_trim14trunc240.qza
