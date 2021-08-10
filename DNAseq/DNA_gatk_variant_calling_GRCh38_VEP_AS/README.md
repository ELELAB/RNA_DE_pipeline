# Snakemake workflow: DNA gatk Alelle Specific (AS) variant calling

This is a modified working version of the "standard" workflow from https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling.
The majority of the modifications are introduced in order to run VQSR.
As the reference genome, primary assembly GRCh38 is used.
As known variant sets, modified databases (removed "chr" prefix; not for alternative contigs yet) from gatk bundle are used.
Biological annotations are based on VEP.
To run this workflow, resources from gatk bundle must be downloaded or cached based on existing resource.

This pipeline is using standard group annotations (-G)

This version of gatk AS variant calling is not working properly yet. It's NOT part of gatk best practices now (it's beta version).
Specific problem in this version is: 
Merging snvs and indels after recalibration has different format (not including FILTER info) in comparison to the standard pipeline.

Pipeline needs further testing and debugging plus further should be considered:

This pipeline will be most likely replaced by the pipeline based on:
1) hg38 reference genome with .alt index (reference with alternative contigs - HLA, etc.)
2) Funcotator for variant annotations
3) VEP annotations might be added to the modified ("chr" removed) vcf file


Before starting the analysis, resources must be downloaded using download_resources.sh
Additional quality metrics could be calculated using qc.sh
