#!/bin/bash

#conda init --all

#conda create --name gs_utils
#conda activate gs_utils
#conda install gs_utils -c bioconda

cd ./resources

gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/GRCh38.primary_assembly.genome.fa ./
gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/GRCh38_gencode.v27.refFlat.txt ./
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz

gunzip hg38.refGene.gtf.gz

ln -sT ./GRCh38.primary_assembly.genome.fa genome.fasta
ln -sT ./GRCh38_gencode.v27.refFlat.txt genome.refFlat.txt
ln -sT ./hg38.refGene.gtf genome.gtf

#genome_noAlt.gtf is produced in order to run featureCounts generating results compatible with STAR counts output (primary for the purpose of comparing the results in order evaluate their similarity).

#cd ./resources
#sed '/_alt/d' genome.gtf > genome_noAlt.gtf
#sed -i '/_decoy/d' genome_noAlt.gtf 
#sed -i '/_random/d' genome_noAlt.gtf 
#sed -i '/^HLA/d' genome_noAlt.gtf 
#sed -i '/^chrUn_/d' genome_noAlt.gtf 

