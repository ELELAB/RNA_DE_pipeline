#!/bin/bash

###The following lines create local conda installation of the gsutils package which is needed to download reference resources.

#conda init --all
#conda create --name gs_utils
#conda activate gs_utils
#conda install -c conda-forge gsutil

cd ./resources

gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/GRCh38.primary_assembly.genome.fa ./
gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/GRCh38_gencode.v27.refFlat.txt ./
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz

gunzip hg38.refGene.gtf.gz

ln -sT ./GRCh38.primary_assembly.genome.fa genome.fasta
ln -sT ./GRCh38_gencode.v27.refFlat.txt genome.refFlat.txt
ln -sT ./hg38.refGene.gtf genome.gtf
