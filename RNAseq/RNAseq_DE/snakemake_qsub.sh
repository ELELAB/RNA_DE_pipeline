#!/bin/sh
### Note: No commands may be executed until after the #PBS lines
### Account information
#PBS -W group_list=dtu_00011 -A dtu_00011
### Job name (comment out the next line to get the name of the script used as the job name)
#PBS -N nik_snakemake_rna
### Output files (comment out the next 2 lines to get the job name used instead)
#PBS -e rnaseq_test.err
#PBS -o rnaseq_test.log
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=38
### Memory
#PBS -l mem=150gb
### Requesting time - format is <days>:<hours>:<minutes>:<seconds> (here, 12 hours)
#PBS -l walltime=120:00:00

# Go to the directory from where the job was submitted (initial directory is $HOME)
MY_dir=/home/projects/dtu_00011/data/icope_analysis/ALL/rna/rnaseq/gene_expression/runs/initial_run_93_samples/code/ngs_pipeline/RNAseq/RNAseq_DE
echo Working directory is $MY_dir
cd $MY_dir

### Here follows the user commands:
# Define number of processors
NPROCS=30
echo This job has allocated $NPROCS nodes

# Load all required modules for the job
module load tools
module load miniconda3/4.10.1
conda init --all
conda activate /home/people/niktom/.conda/envs/snakemake/


# This is where the work is done
# Make sure that this script is not bigger than 64kb ~ 150 lines, otherwise put in seperat script and execute from here
#export SNAKEMAKE_OUTPUT_CACHE=/home/people/niktom/nik_projects/snakemake/test_RNAseq/rna-seq-star-deseq2/resources
#snakemake --cores 30 --use-conda --cache get_genome get_annotation genome_faidx bwa_index star_index star_genome
snakemake --cores 38 --use-conda --rerun-incomplete
