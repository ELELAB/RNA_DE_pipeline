#!/bin/bash

###conda init --all
###conda activate nik_tools


###Define variables for indexed genome
REF_SEQ=/home/projects/dtu_00011/people/niktom/snakemake/test_gatk_local/dna-seq-gatk-variant-calling/resources/genome.fasta
REF_SEQ_DIR=/home/projects/dtu_00011/people/niktom/snakemake/test_gatk_local/dna-seq-gatk-variant-calling/resources
DICT=/home/projects/dtu_00011/people/niktom/snakemake/test_gatk_local/dna-seq-gatk-variant-calling/resources/genome.dict

###Define input/output dir variables
INPUT_DIR=/home/projects/dtu_00011/people/niktom/snakemake/test_gatk_local/dna-seq-gatk-variant-calling/results/trimmed/
RESULTS_QC=/home/projects/dtu_00011/people/niktom/snakemake/test_gatk_local/dna-seq-gatk-variant-calling/results/QC2

DEDUP=/home/people/niktom/nik_projects/snakemake/test_gatk_local/dna-seq-gatk-variant-calling/results/dedup
CODING_REGIONS_BED=/home/people/niktom/nik_projects/snakemake/resources/giab_stratification/ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v2.0/GRCh38/FunctionalRegions/edited_chr/GRCh38_refseq_cds.bed


##############################################################


FASTQC=$RESULTS_QC/fastqc
MULTIQC=$RESULTS_QC/multiqc 
PICARD_METRICS=$RESULTS_QC/picard
BEDTOOLS=$RESULTS_QC/bedtools
SAMTOOLS=$RESULTS_QC/samtools

###Create results directories
echo "Create results DIRs"
mkdir -p $RESULTS_QC
mkdir -p $FASTQC
mkdir -p $MULTIQC
mkdir -p $PICARD_METRICS
mkdir -p $BEDTOOLS
mkdir -p $SAMTOOLS
echo "Results DIRs created"


###Define suffixes
APPENDIX=".fastq.gz" 
APPENDIX1=".1.fastq.gz" # File suffix for forward reads
APPENDIX2=".2.fastq.gz" # File suffix for reverse reads

for file in $INPUT_DIR/*$APPENDIX1
do
    ###Define sample variables
    FORWARD=$file
    REVERSE=${FORWARD%$APPENDIX1*}$APPENDIX2
    SAMPLE="$(basename ${FORWARD%$APPENDIX1})"


    ###FASTQC on trimmed reads
    echo "Run fastqc on original set of reads"
    fastqc -t 10 $FORWARD $REVERSE -o $FASTQC/
    echo "Run fastqc on original set of reads done"


    ###Coverage analysis using bedtools
    echo "Calculate coverage histogram and genomecoverage for sample $SAMPLE"
    bedtools coverage -hist -abam $DEDUP/$SAMPLE.bam -b $CODING_REGIONS_BED | grep ^all > $BEDTOOLS/$SAMPLE.bam.rmdup.hist.all.txt
    bedtools genomecov -ibam $DEDUP/$SAMPLE.bam -bga > $BEDTOOLS/$SAMPLE.genomecov 
    echo "Calculate coverage histogram and genomecoverage for sample $SAMPLE done"

    echo "calculate picard for deduplicated alignments"
    ###ADD HEADER TO INTERVAL DATA
    picard BedToIntervalList \
    I=$CODING_REGIONS_BED \
    O=$CODING_REGIONS_BED.list.interval_list \
    SD=$DICT
    
    ###COLLECT HS METRICS FOR CODING REGIONS
    picard CollectHsMetrics I=$DEDUP/${SAMPLE}.bam O=$PICARD_METRICS/${SAMPLE}.hs_metrics.txt PER_TARGET_COVERAGE=$PICARD_METRICS/${SAMPLE}.hs_perTargetCov.txt R=$REF_SEQ BAIT_INTERVALS=$CODING_REGIONS_BED.list.interval_list TARGET_INTERVALS=$CODING_REGIONS_BED.list.interval_list

    ###COLLECT WGS METRICS
    picard CollectWgsMetrics \
    I=$DEDUP/$SAMPLE.bam \
    O=$PICARD_METRICS/$SAMPLE.collect_wgs_metrics.txt \
    R=$REF_SEQ
    
    ###SAMTOOLS STATISTICS - STAT AND FLAGSTAT
    samtools stats $DEDUP/${SAMPLE}.bam > $SAMTOOLS/${SAMPLE}.stats
    samtools flagstat $DEDUP/${SAMPLE}.bam > $SAMTOOLS/${SAMPLE}.flagstat
        
    ###DO MULTIQC ON THE QC RESULTS    
    multiqc $RESULTS_QC 


###NEXT, QC ON VEP RESULTS SHOULD BE RUN
 #   echo "Run Variant Effect Predictor on $SAMPLE"
 #   mkdir -p $VARCALL/vep
 #   vep -i $VARCALL/$SAMPLE.vcf -o $VARCALL/vep/$SAMPLE.vep.vcf --fasta $REF_SEQ --force_overwrite -vcf --database --assembly GRCh38 --port 3337
 
 ###FOR maftools, multisample vcf would need to be split per sample, converted to maf and those maf files concatenated.
 #   echo "Run maf2vcf on $SAMPLE"
 #   mkdir -p $VARCALL/maf
 #   vcf2maf.pl --input-vcf $VARCALL/vep/$SAMPLE.vcf --output-maf $VARCALL/maf/$SAMPLE.maf --ref-fasta $REF_SEQ --inhibit-vep 
 #   echo "Run maf2vcf on $SAMPLE done"

done

