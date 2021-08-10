# Snakemake workflow: rna-seq-star-deseq2

<!-- [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4737358.svg)](https://doi.org/10.5281/zenodo.4737358)
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.1.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/snakemake-workflows/rna-seq-star-deseq2/workflows/Tests/badge.svg?branch=master)](https://github.com/snakemake-workflows/rna-seq-star-deseq2/actions?query=branch%3Amaster+workflow%3ATests) -->

This pipeline is a modified version of the "best practices" snakemake workflow written by Johannes Koester:
https://github.com/snakemake-workflows/rna-seq-star-deseq2
This workflow performs a differential gene expression analysis using Cutadapt, STAR and Deseq2. Alternatively, genes (or other features) are counted also using FeatureCounts. 
Report.html file offers an interactive workflow diagram in the form of Direct Acyclic Graph (DAG); clicking on the node will display code for each step (code visualization might not work well for robust workflows). Report include also run time statistics, some details about pipeline configuration (from confing.yaml) and selected results (PCA, DEG graphs). Example MultiQC output is represented by multiqc_report.html


## Usage

The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/?usage=snakemake-workflows%2Frna-seq-star-deseq2). As our code is a modified version, please, skip Step 2 (deploy workflow) and download the code from this github pages.



