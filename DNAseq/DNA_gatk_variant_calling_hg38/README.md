# Snakemake workflow: DNA gatk variant calling based on gatk bundle (hg38) and funcotator

This is a modified working version of the "standard" workflow from https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling.
The majority of the modifications are introduced in order to run VQSR.
As the reference genome, primary assembly hg38 is used.
As known variant sets from gatk bundle are used.
To run this workflow, resources from gatk bundle must be downloaded prior the analysis or cached based on existing resource.



This pipeline is based on:
1) hg38 reference genome with .alt index (reference with alternative contigs - HLA, etc.)
2) Funcotator for variant annotations (config files for the databases needs to be set up; probably other debuging needed - pipeline run finished)
3) VEP annotations might be added to the modified ("chr" removed) vcf file


Before starting the analysis, resources must be downloaded using download_resources.sh
Additional quality metrics could be calculated using qc.sh
