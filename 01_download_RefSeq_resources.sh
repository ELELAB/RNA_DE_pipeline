#!/bin/bash

ln -sT /home/projects/dtu_00011/data/icope_analysis/ALL/rna/rnaseq/resources_RefSeq ./resources

# ###The following lines create local conda installation of the gsutils package which is needed to download reference resources. 
#conda init --all
#conda create --name ref_tools
#conda activate ref_tools
#conda install -c conda-forge gsutil
#conda install -c conda-forge -c bioconda picard
#conda install -c conda-forge -c bioconda samtools



# ###Download reference genome and annotations
# ###REFERENCE GENOME:
# gsutil -m cp -r gs://genomics-public-data/references/hg38/v0/GRCh38.primary_assembly.genome.fa ./

# ###Alternatively, GENCODE with alternative contigs or RefSeq annotations might be used.
# ###GENCODE is supposed to be more comprehensive and better for the analysis of regulatory regions; GENCODE is used by GDC

# ###REFSEQ:
# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.refGene.gtf.gz  ###	UCSC states it's based on NCBI Refseq
# wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refFlat.txt.gz ./  ###the same data as download using UCSC tables - track (NCBI refseq) and table (refflat)
# gunzip hg38.refGene.gtf.gz
# gunzip refFlat.txt.gz

# ###Create links to reference genome and annotations
# ln -sT ./GRCh38.primary_assembly.genome.fa genome.fasta
# ln -sT ./refFlat.txt genome.refFlat.txt
# ln -sT ./hg38.refGene.gtf genome.gtf


# ###Create reference sequence dictionary and index for picard
# picard CreateSequenceDictionary R=genome.fasta O=genome.dict
# samtools faidx genome.fasta

# ###rRNA intervals - These files are necessary for additional QC (picard); interval list was created using picard
# ###get rRNA.gtf file from UCSC Table Browser:

# #http://genome.ucsc.edu/cgi-bin/hgTables
# #Select "All Tables" from the group drop-down list
# #Select the "rmsk" table from the table drop-down list
# #Choose "GTF" as the output format
# #Type a filename in "output file" so your browser downloads the result
# #Click "create" next to filter
# #Next to "repClass," type rRNA
# #Next to free-form query, select "OR" and type repClass = "tRNA"
# #Click submit on that page, then get output on the main page

# ###RefSeq rRNA copied manually
# ###only chromosomes 1-23 kept
# picard BedToIntervalList I=rRNA_HG38_UCSC_tables1.bed O=rRNA_HG38_UCSC_tables1.interval_list SD=genome.dict
# ln -sT ./rRNA_HG38_UCSC_tables1.interval_list rRNA.interval_list

# # ###get stratification regions include file with CDS; downloaded from:
# wget -r ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/GRCh38/FunctionalRegions/
# gunzip ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/GRCh38/FunctionalRegions/GRCh38_refseq_cds.bed.gz
# picard BedToIntervalList I=ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/GRCh38/FunctionalRegions/GRCh38_refseq_cds.bed O=GRCh38_refseq_cds.interval_list SD=genome.dict

