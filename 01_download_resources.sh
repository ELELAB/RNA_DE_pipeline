ln -s /home/projects/dtu_00011/data/icope_analysis/ALL/rna/rnaseq/resources ./

# #!/bin/bash



# ###The following lines create local conda installation of the gsutils package which is needed to download reference resources.

# #conda init --all
# #conda create --name ref_tools
# #conda activate ref_tools
# #conda install -c conda-forge gsutil
# #conda install -c conda-forge picard
# #conda install -c conda-forge samtools

# cd ./resources

# ###Download reference genome and annotations
# gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/GRCh38.primary_assembly.genome.fa ./
# gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/GRCh38_gencode.v27.refFlat.txt ./
# wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
# gunzip gencode.v38.annotation.gtf.gz

# ###Alternatively, GENCODE with alternative contigs or RefSeq annotations might be used.
# ###GENCODE is supposed to be more comprehensive and better for the analysis of regulatory regions; GENCODE is used by GDC
# ##wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.chr_patch_hapl_scaff.annotation.gtf.gz
# ##wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz

# ###Create links to reference genome and annotations
# ln -sT ./GRCh38.primary_assembly.genome.fa genome.fasta
# ln -sT ./GRCh38_gencode.v27.refFlat.txt genome.refFlat.txt
# ln -sT ./gencode.v38.annotation.gtf genome.gtf

# ###Create reference sequence dictionary and index for picard
# picard CreateSequenceDictionary R=genome.fasta O=genome.dict
# samtools faidx genome.fasta

# ###rRNA intervals - These files are necessary for additional QC (picard); interval list was created using picard
# cat gencode.v38.annotation.gtf | grep rRNA | cut -s -f 1,4,5,7,9 > gencode.v38.annotation_rRNA.bed
# picard BedToIntervalList I=gencode.v38.annotation_rRNA.bed O=gencode.v38.annotation_rRNA.interval_list SD=genome.dict
# ###Alternatively, "hg38_primary_ribosomal.interval_list" which seems to be more comprehensive (but cannot find the source - some forum) could be used.

# ###get stratification regions include file with CDS; downloaded from:
# wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/GRCh38/FunctionalRegions/
# gunzip ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/GRCh38/FunctionalRegions/GRCh38_refseq_cds.bed.gz
# picard BedToIntervalList I=ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/GRCh38/FunctionalRegions/GRCh38_refseq_cds.bed O=GRCh38_refseq_cds.interval_list SD=genome.dict

