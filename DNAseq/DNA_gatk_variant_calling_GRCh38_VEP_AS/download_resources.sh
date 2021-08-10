#!/bin/bash

conda init --all

#conda create --name gs_utils
#conda activate gs_utils
#conda install gs_utils -c bioconda
#conda install rename #not tested

RESOURCES=/home/people/niktom/nik_projects/snakemake/resources/
cd $RESOURCES
mkdir -p gatk_bundle
cd gatk_bundle

gsutil -m cp -r gs://genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf \
gs://genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.idx \
gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz \
gs://genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi \
gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
gs://genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi \
gs://genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
gs://genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi \
gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf \
gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx \
gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict \
gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta \
gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.alt \
gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.amb \
gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.ann \
gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.bwt \
gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.pac \
gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.64.sa \
gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta.fai \
gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz \
gs://genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi \
gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
gs://genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi \
gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz \
gs://genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi \
gs://genomics-public-data/resources/broad/hg38/v0/scattered_calling_intervals/ \
gs://genomics-public-data/resources/broad/hg38/v0/wgs_calling_regions.hg38.interval_list .


gsutil -m cp -r  “gs://genomics-public-data/resources/broad/hg38/v0/*.*“
cd ..


ln -sT ./gatk_bundle/Homo_sapiens_assembly38.dict genome.dict
ln -sT ./gatk_bundle/Homo_sapiens_assembly38.fasta genome.fasta
ln -sT ./gatk_bundle/Homo_sapiens_assembly38.fasta.64.alt genome.fasta.alt
ln -sT ./gatk_bundle/Homo_sapiens_assembly38.fasta.64.amb genome.fasta.amb
ln -sT ./gatk_bundle/Homo_sapiens_assembly38.fasta.64.ann genome.fasta.ann
ln -sT ./gatk_bundle/Homo_sapiens_assembly38.fasta.64.bwt genome.fasta.bwt
ln -sT ./gatk_bundle/Homo_sapiens_assembly38.fasta.64.pac genome.fasta.pac
ln -sT ./gatk_bundle/Homo_sapiens_assembly38.fasta.64.sa genome.fasta.sa
ln -sT ./gatk_bundle/Homo_sapiens_assembly38.fasta.fai genome.fasta.fai


wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/

gsutil -m cp -r gs://broad-public-datasets/funcotator/ .


###using renamed version of database is not recommended; 
###CHROMOSOMES ARE RENAMED BY SUBTRACTING "CHR". OTHER THAN CHROMOSOMES 1-23 ARE NOT COMPATIBLE WITH ENSEMBL REFERENCE (NEEDS ANOTHER CURATION)
mkdir -p edited_chr
cd edited_chr
for file in ../gatk_bundle/*.vcf.gz
do
filename=$(basename "$file")
echo "processing $filename"
zcat $file | sed 's/chr//g' > $filename.tmp
rename 's/.gz.tmp//g' ./*
done

for file in ./*.vcf
do
bgzip $file
tabix $file.gz
done

for file in ../gatk_bundle/*.vcf
do
filename=$(basename "$file")
echo "processing $filename"
cat $file | sed 's/chr//g' > $filename
done


